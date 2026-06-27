"""
patch_algos.py  --  Log the optimizer's cost history (convergence helper)
=========================================================================

PURPOSE
-------
The ORR optimizer (`orr_recon` / `_orr_cpu`) computes a data-fidelity cost at
every gradient-descent step but, by default, returns only the final image. For
this study we wanted to *verify that the optimizer converges* on dense-breast
data -- i.e. that the cost drops substantially across iterations. That requires
access to the per-iteration cost history.

This helper applies a small, surgical edit to a LOCAL working copy of
`umbms/recon/algos.py` so that `orr_recon` returns `(img, cost_history)`. It is a
convenience for evaluation; it does not change the reconstruction algorithm.

DESIGN NOTES (so the edit is safe to apply and re-apply)
--------------------------------------------------------
  * It restores a clean copy from git first, so running it repeatedly is safe.
  * It locates `_orr_cpu` and `orr_recon` BY NAME and edits only those, bounded by
    the next function -- so the beamformers `fd_dmas` / `fd_das` (which also end in
    `return img`) are never touched.
  * `orr_recon` has a GPU branch and a CPU branch sharing one return statement; the
    GPU branch does not track cost, so an empty cost list is injected there to keep
    the shared return valid.

Call sites should guard with `isinstance(out, tuple)` -- see
`contrib/reconstruction usage` in the notebook / README.

USAGE
-----
    python contrib/patch_algos.py /path/to/ORR-EPM
"""
import subprocess
import re
import sys
import os


def patch_algos(repo_path):
    algos = os.path.join(repo_path, "umbms/recon/algos.py")

    # Restore a pristine copy first so re-running is always safe.
    subprocess.run(["git", "-C", repo_path, "checkout", "--",
                    "umbms/recon/algos.py"], check=False)

    lines = open(algos).read().split("\n")

    def find_def(name):
        for i, l in enumerate(lines):
            if l.startswith(f"def {name}("):
                return i
        return None

    i_orr_recon = find_def("orr_recon")
    i_orr_cpu = find_def("_orr_cpu")
    i_fd_das = find_def("fd_das")  # boundary: orr code ends before this

    # Detect the cost-history variable name used inside _orr_cpu.
    cost_var = None
    for i in range(i_orr_cpu, i_fd_das):
        m = re.search(r"(\w*cost\w*)\s*=", lines[i])
        if m:
            cost_var = m.group(1)
            break
    if cost_var is None:
        raise RuntimeError("Could not locate cost variable inside _orr_cpu")

    def patch_return(start, end, new_return):
        for i in range(end - 1, start, -1):
            if lines[i].strip() == "return img":
                indent = lines[i][:len(lines[i]) - len(lines[i].lstrip())]
                lines[i] = indent + new_return
                return i
        return None

    # _orr_cpu:  return img  ->  return img, cost_funcs
    patch_return(i_orr_cpu, i_fd_das, f"return img, {cost_var}")

    # orr_recon CPU branch: capture the cost returned by _orr_cpu(...)
    for j in range(i_orr_recon, i_orr_cpu):
        if "_orr_cpu(" in lines[j] and "=" in lines[j].split("_orr_cpu(")[0]:
            lhs = lines[j].split("=")[0].strip()
            indent = lines[j][:len(lines[j]) - len(lines[j].lstrip())]
            rhs = lines[j].split("=", 1)[1].strip()
            lines[j] = f"{indent}{lhs}, _cost_hist = {rhs}"
            break

    # orr_recon GPU branch: define an empty _cost_hist for the shared return.
    for j in range(i_orr_recon, i_orr_cpu):
        if "img = _orr_gpu(" in lines[j]:
            k = j
            while lines[k].strip() != ")":
                k += 1
            indent = lines[j][:len(lines[j]) - len(lines[j].lstrip())]
            lines.insert(k + 1, f"{indent}_cost_hist = []  # GPU path: no cost history")
            break

    # orr_recon shared return:  return img  ->  return img, _cost_hist
    patch_return(i_orr_recon, i_orr_cpu, "return img, _cost_hist")

    open(algos, "w").write("\n".join(lines))
    print("orr_recon now returns (img, cost_history); beamformers untouched.")


if __name__ == "__main__":
    repo = sys.argv[1] if len(sys.argv) > 1 else "."
    patch_algos(repo)
