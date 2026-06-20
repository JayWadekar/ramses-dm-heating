from pathlib import Path

import numpy as np
from scipy.io import FortranFile


def read_hydro(path):
    with FortranFile(path, "r") as f:
        ncpu = f.read_ints(np.int32)[0]
        nvar = f.read_ints(np.int32)[0]
        ndim = f.read_ints(np.int32)[0]
        nlevelmax = f.read_ints(np.int32)[0]
        nboundary = f.read_ints(np.int32)[0]
        _gamma = f.read_reals(np.float64)[0]

        fields = [[] for _ in range(nvar)]
        for _level in range(nlevelmax):
            for _domain in range(nboundary + ncpu):
                _ilevel = f.read_ints(np.int32)[0]
                ncache = f.read_ints(np.int32)[0]
                if ncache == 0:
                    continue
                for _cell in range(2**ndim):
                    for ivar in range(nvar):
                        fields[ivar].append(f.read_reals(np.float64))

    return {
        "density": np.concatenate(fields[0]),
        "pressure": np.concatenate(fields[4]),
    }


def summarize(name, delta):
    return [
        f"{name}_n_cells {delta.size}",
        f"{name}_max_abs_delta {np.max(np.abs(delta)):.12e}",
        f"{name}_mean_delta {np.mean(delta):.12e}",
        f"{name}_median_delta {np.median(delta):.12e}",
        f"{name}_min_delta {np.min(delta):.12e}",
        f"{name}_max_delta {np.max(delta):.12e}",
    ]


if __name__ == "__main__":
    here = Path(__file__).resolve().parent
    stock = read_hydro(here / "run_stock" / "output_00003" / "hydro_00003.out00001")
    dmheat = read_hydro(here / "run_dmheat" / "output_00003" / "hydro_00003.out00001")

    lines = []
    for field in ("density", "pressure"):
        lines.extend(summarize(field, dmheat[field] - stock[field]))

    print("\n".join(lines))
