# PyInstaller recipe for Magma.exe: one file holding Python, numpy, scipy, the app and its data tables.
# Build with packaging/build_exe.ps1 (or: pyinstaller packaging/magma.spec from the repository root).
from pathlib import Path

ROOT = Path(SPECPATH).parent
SRC = ROOT / "src" / "magma"

a = Analysis(
    [str(ROOT / "packaging" / "magma_exe.py")],
    pathex=[str(ROOT / "src")],
    datas=[(str(SRC / "data" / "*.tsv"), "magma/data"),
           (str(SRC / "app" / "static"), "magma/app/static")],
    excludes=["matplotlib", "pandas", "tkinter", "pytest", "IPython"],
)
pyz = PYZ(a.pure)
exe = EXE(
    pyz, a.scripts, a.binaries, a.datas,
    name="Magma",
    console=True,          # the window shows where the app runs; closing it stops the app
    upx=False,
    icon=str(ROOT / "packaging" / "magma.ico"),            # the Magma M in the Dougalien Systems style (packaging/icon/make_icons.py)
    version=str(ROOT / "packaging" / "version_info.txt"),  # company and product in Properties > Details
)
