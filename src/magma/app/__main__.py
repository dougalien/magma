"""``python -m magma.app`` starts the local server and opens the control window."""
from __future__ import annotations

import argparse
import multiprocessing
import os
import sys
import threading
import webbrowser

# One thread for numpy's linear algebra: Magma's matrices are small, and on Windows several threads
# fighting over them made each step many times slower (the app looked frozen).  Set before numpy loads.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from .server import DEFAULT_PORT, make_server  # noqa: E402


def _no_quick_edit() -> None:
    """Windows consoles pause the whole program while text in them is selected ("QuickEdit"), so one
    click in the Magma console froze the app.  Turn that off for this console."""
    try:
        import ctypes
        k = ctypes.windll.kernel32
        h = k.GetStdHandle(-10)                    # standard input
        mode = ctypes.c_uint32()
        if k.GetConsoleMode(h, ctypes.byref(mode)):
            k.SetConsoleMode(h, (mode.value & ~0x0040) | 0x0080)   # clear QUICK_EDIT, set EXTENDED_FLAGS
    except Exception:
        pass


def main() -> None:
    multiprocessing.freeze_support()
    ap = argparse.ArgumentParser(description="Magma: crystallization of silicate melts, by Dougalien Systems")
    ap.add_argument("--port", type=int, default=DEFAULT_PORT)
    ap.add_argument("--no-browser", action="store_true", help="do not open the windows")
    args = ap.parse_args()
    if sys.platform == "win32":               # the console window says Magma, not the path of the exe
        import ctypes
        ctypes.windll.kernel32.SetConsoleTitleW("Magma, by Dougalien Systems")
        _no_quick_edit()
    httpd = make_server(args.port)
    url = f"http://127.0.0.1:{httpd.server_address[1]}"
    print("Magma, by Dougalien Systems (dougalien.com)")
    print(f"Magma is running at {url}/control  (results window: {url}/results)")
    print("Press Ctrl+C to stop.")
    if not args.no_browser:
        threading.Timer(0.5, lambda: webbrowser.open(f"{url}/control")).start()
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        httpd.server_close()


if __name__ == "__main__":
    main()
