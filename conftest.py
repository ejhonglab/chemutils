"""
Configures a temporary chemutils cache (via the CHEMUTILS_CACHE environment variable)
for all tests, so they never read or modify the real cache, and blocks network access
in tests not marked with `@pytest.mark.network`.
"""

import os
from pathlib import Path
import shutil
import socket
import sqlite3
import tempfile

import pytest


_session_cache_dir = None

def pytest_configure(config):
    global _session_cache_dir
    # Set before anything (e.g. doctest collection) imports chemutils, so even code run
    # outside of the fixture below can not use the real cache.
    _session_cache_dir = tempfile.mkdtemp(prefix='chemutils_test_cache_')
    os.environ['CHEMUTILS_CACHE'] = str(Path(_session_cache_dir) / 'cache.sqlite')


def pytest_unconfigure(config):
    if _session_cache_dir is not None:
        shutil.rmtree(_session_cache_dir, ignore_errors=True)


class NetworkAccessBlocked(RuntimeError):
    # Not an OSError subclass, so libraries do not catch it as a connection error.
    pass


class NetworkGuard:
    """Records (and raises on) network access while `blocking` is True.
    """
    def __init__(self, monkeypatch, blocking: bool):
        self.blocking = blocking
        self.attempts = []

        orig_getaddrinfo = socket.getaddrinfo
        orig_connect = socket.socket.connect

        def check(address):
            if self.blocking:
                self.attempts.append(address)
                raise NetworkAccessBlocked(f'test tried to access network: {address}')

        def getaddrinfo(host, port, *args, **kwargs):
            check((host, port))
            return orig_getaddrinfo(host, port, *args, **kwargs)

        def connect(sock, address):
            check(address)
            return orig_connect(sock, address)

        monkeypatch.setattr(socket, 'getaddrinfo', getaddrinfo)
        monkeypatch.setattr(socket.socket, 'connect', connect)


@pytest.fixture(autouse=True)
def chemutils_test_cache(tmp_path, monkeypatch):
    """Points chemutils at a new, empty cache for each test. Yields its path.
    """
    import chemutils as cu

    cache_path = tmp_path / 'chemutils_cache.sqlite'
    monkeypatch.setenv('CHEMUTILS_CACHE', str(cache_path))

    real_cache = Path(cu.default_cache_file).resolve()
    orig_connect = sqlite3.connect

    def connect(database, *args, **kwargs):
        assert Path(database).resolve() != real_cache, 'test tried to open real cache'
        return orig_connect(database, *args, **kwargs)

    monkeypatch.setattr(sqlite3, 'connect', connect)

    cu._reset_cache_state()
    yield cache_path
    cu._reset_cache_state()


@pytest.fixture(autouse=True)
def network_guard(request, monkeypatch):
    """Blocks network access, unless test is marked `network`.

    Tests marked `network` can set `network_guard.blocking = True` (e.g. to check a
    second call is served from the cache), and check `network_guard.attempts`.
    """
    blocking = request.node.get_closest_marker('network') is None
    return NetworkGuard(monkeypatch, blocking)
