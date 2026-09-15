"""Boundary experiments never contact external servers."""

import socket

import pytest


@pytest.fixture(autouse=True)
def forbid_live_network(monkeypatch: pytest.MonkeyPatch) -> None:
    def denied(*args: object, **kwargs: object) -> None:
        raise AssertionError(
            "Boundary tests require an injected HTTP response, not live network access"
        )

    monkeypatch.setattr(socket.socket, "connect", denied)
