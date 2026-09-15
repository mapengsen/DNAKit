"""Agent-facing discovery and execution for stable DNAKit public functions."""

from dnakit._verbosity import register_public as _register_public_verbose

from .agent import AgentResult, run_agent
from .catalog import AgentToolSpec, ToolEffect, build_tool_catalog
from .embeddings import EmbeddingConfig, ToolSearchMode
from .execution import PythonExecutionConfig
from .registry import AgentToolRegistry, default_tool_registry
from .serialization import AgentSerializationLimits, serialize_for_agent
from .server import create_server

__all__ = [
    "AgentResult",
    "AgentSerializationLimits",
    "AgentToolRegistry",
    "AgentToolSpec",
    "EmbeddingConfig",
    "PythonExecutionConfig",
    "ToolEffect",
    "ToolSearchMode",
    "build_tool_catalog",
    "create_server",
    "default_tool_registry",
    "run_agent",
    "serialize_for_agent",
]

_register_public_verbose(globals())
