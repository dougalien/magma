"""Magma on Streamlit Community Cloud: ``streamlit run streamlit_app.py``.

Two windows, as in the desktop app: Input (this address) and Results (``/results``), sharing one run
through the private ``?w=`` id in the address.
"""
import os
import sys
from pathlib import Path

# One thread for numpy's linear algebra (Magma's matrices are small); set before numpy loads.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

import streamlit as st  # noqa: E402

from magma.web.app import input_page, results_page  # noqa: E402

page = st.navigation([st.Page(input_page, title="Input", url_path="input", default=True),
                      st.Page(results_page, title="Results", url_path="results")], position="hidden")
page.run()
