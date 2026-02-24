from __future__ import annotations


class PrimerMakerError(RuntimeError):
    """Base exception for Primer Maker."""


class ToolError(PrimerMakerError):
    """Wrap all transport/API related errors."""


class NoMappingError(PrimerMakerError):
    """No mapping could be resolved for the given UniProt accession."""


class SequenceLengthMismatchError(PrimerMakerError):
    """Fetched sequence does not match expected coordinate span."""
