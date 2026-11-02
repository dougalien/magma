"""The Magma M is on every page, as the browser-tab icon and at /favicon.ico."""
import http.client
import threading

import pytest

from magma.app.server import make_server


@pytest.fixture(scope="module")
def server():
    httpd = make_server(0)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    yield httpd
    httpd.shutdown()
    httpd.server_close()


def get(server, path):
    conn = http.client.HTTPConnection("127.0.0.1", server.server_address[1], timeout=10)
    conn.request("GET", path)
    r = conn.getresponse()
    return r.status, r.read()


def test_pages_use_the_magma_m(server):
    for path in ("/control", "/results"):
        page = get(server, path)[1].decode()
        assert "/static/magma_mark.svg" in page and "/favicon.ico" in page, path


def test_favicon_is_the_magma_ico(server):
    status, data = get(server, "/favicon.ico")
    assert status == 200 and data[:4] == b"\x00\x00\x01\x00"     # a Windows .ico, not a 404
