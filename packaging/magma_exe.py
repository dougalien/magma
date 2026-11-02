"""Entry point of the packaged Magma.exe (see magma.spec)."""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")      # before numpy loads (see magma.app.__main__)

from magma.app.__main__ import main  # noqa: E402

if __name__ == "__main__":
    main()
