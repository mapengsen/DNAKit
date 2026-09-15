# Adapted from ammawla/encode-toolkit @ b1d90bc9e7ead40d49b8c89f0b9b9e7b8e592a5b.
# Original file: src/encode_connector/client/__init__.py. See LICENSE and UPSTREAM.json at the package root.
# DNAKit changes: package namespace, isolated cache, lazy imports and native execution.
"""ENCODE Project Python client library."""

from dnakit.encode._backend.client.auth import CredentialManager
from dnakit.encode._backend.client.encode_client import EncodeClient

__all__ = ["EncodeClient", "CredentialManager"]
