"""Local registration protocol for DNAKit's licensed scientific adapters.

No plugin discovery, upstream package import, or tool execution happens here.
"""

from __future__ import annotations

import inspect

_classes: dict[str, type] = {}


def register_tool(tool_type_name=None, config=None):
    if inspect.isclass(tool_type_name):
        _classes[tool_type_name.__name__] = tool_type_name
        return tool_type_name

    def register(cls):
        _classes[tool_type_name or cls.__name__] = cls
        return cls

    return register
