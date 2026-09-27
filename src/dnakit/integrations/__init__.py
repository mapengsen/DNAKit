"""Reviewed external genomics analyses with optional, caller-installed backends."""

from dnakit._verbosity import register_public as _register_public_verbose

from .catalog import describe_analysis, list_analyses, search_projects
from .execution import plan_analysis, run_analysis
from .mcp import call_mcp_tool, call_mcp_workflow, list_mcp_servers, list_mcp_tools

__all__ = [
    "call_mcp_tool",
    "call_mcp_workflow",
    "describe_analysis",
    "list_analyses",
    "list_mcp_servers",
    "list_mcp_tools",
    "plan_analysis",
    "run_analysis",
    "search_projects",
]

_register_public_verbose(globals(), tuple(__all__))
