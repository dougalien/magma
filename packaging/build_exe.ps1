# Build dist\Magma.exe on Windows.  Run from the Magma folder:  powershell -ExecutionPolicy Bypass -File packaging\build_exe.ps1
$ErrorActionPreference = "Stop"
python -m pip install --upgrade pip
python -m pip install numpy scipy pyinstaller
python -m PyInstaller --noconfirm --clean packaging\magma.spec
Write-Host "Built dist\Magma.exe"
