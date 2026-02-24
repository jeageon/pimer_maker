from __future__ import annotations


class UTGError(RuntimeError):
    """Base exception for UTG."""


class ToolError(UTGError):
    """Wrap all transport/API related errors."""


class NoMappingError(UTGError):
    """No mapping could be resolved for the given UniProt accession."""


class SequenceLengthMismatchError(UTGError):
    """Fetched sequence does not match expected coordinate span."""

