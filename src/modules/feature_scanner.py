from __future__ import annotations

from typing import Any, Optional

from ..config import DEFAULT_FEATURES, FeatureScanOptions
from ..models.data_schemas import GenomicCoordinates, NegativeFeature
from ..utils import seq_utils
from ..utils.coord_utils import ensembl_to_relative
from ..utils.feature_utils import dedupe_features, merge_by_type
from ..utils.seq_utils import scan_ambiguous, scan_extreme_gc_windows, scan_homopolymers


def _first_non_null(*values: Any) -> Any:
    for value in values:
        if value is not None:
            return value
    return None


def _to_float(value: Any) -> Optional[float]:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _to_int(value: Any, default: Optional[int] = None) -> Optional[int]:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


class FeatureScanner:
    def __init__(self) -> None:
        pass

    def scan(
        self,
        coordinates: GenomicCoordinates,
        full_sequence: str,
        requested_features: Optional[list[str]] = None,
        options: Optional[FeatureScanOptions] = None,
    ) -> tuple[list[NegativeFeature], list[str]]:
        requested_features = requested_features or list(DEFAULT_FEATURES)
        options = options or FeatureScanOptions()
        warnings: list[str] = []
        seq_len = len(full_sequence)
        requested = set(requested_features)

        collected: list[NegativeFeature] = []
        if coordinates.coordinate_source == "ncbi":
            collected.extend(
                self._scan_ncbi_annotations(
                    coordinates=coordinates,
                    requested=requested,
                    full_seq_len=seq_len,
                )
            )
        else:
            warnings.append("NCBI-only pipeline expects coordinate_source='ncbi'")

        collected.extend(
            self._scan_internal(
                full_sequence=full_sequence,
                requested=requested,
                options=options,
            )
        )

        deduped = dedupe_features(collected)
        merge_gaps = {
            "extreme_gc": options.gc_step,
            "homopolymer": 0,
            "ambiguous": 0,
            "annotation": 0,
        }
        normalized = merge_by_type(deduped, merge_gaps=merge_gaps)
        return normalized, warnings

    def _scan_ncbi_annotations(
        self,
        coordinates: GenomicCoordinates,
        requested: set[str],
        full_seq_len: int,
    ) -> list[NegativeFeature]:
        requested_annotations = requested.intersection({"annotation", "annotations", "ncbi"})
        if not requested_annotations:
            return []

        results: list[NegativeFeature] = []
        for feature in coordinates.ncbi_annotations:
            start = _to_int(feature.get("start"), 0)
            end = _to_int(feature.get("end"), 0)
            if start <= 0 or end <= 0:
                continue
            rel_start, rel_end = ensembl_to_relative(start, end, coordinates.ext_start_1based, full_seq_len)
            if rel_start >= rel_end:
                continue

            feature_type = _first_non_null(feature.get("feature_type"), "annotation")
            display_name = _first_non_null(feature.get("display_name"), "")
            source_type = _first_non_null(feature.get("source_type"), "ncbi_feature")
            desc = f"NCBI annotation: {feature_type}"
            if display_name:
                desc = f"NCBI annotation: {display_name} ({feature_type})"

            qualifiers = feature.get("qualifiers") or {}
            if not isinstance(qualifiers, dict):
                qualifiers = {"raw": str(qualifiers)}

            attributes = {
                "source_type": source_type,
                "feature_type_raw": feature_type,
            }
            for key, val in qualifiers.items():
                if val:
                    attributes[str(key)] = str(val)

            results.append(
                NegativeFeature(
                    feature_type="annotation",
                    start=rel_start,
                    end=rel_end,
                    description=desc,
                    source="ncbi_genbank",
                    strand=_to_int(feature.get("strand"), 0),
                    attributes=attributes,
                )
            )
        return results

    def _scan_internal(
        self,
        full_sequence: str,
        requested: set[str],
        options: FeatureScanOptions,
    ) -> list[NegativeFeature]:
        results: list[NegativeFeature] = []
        seq_len = len(full_sequence)
        if not requested:
            return results

        if "extreme_gc" in requested:
            windows = scan_extreme_gc_windows(
                full_sequence,
                window_size=options.gc_window,
                step=options.gc_step,
                gc_min=options.gc_min,
                gc_max=options.gc_max,
            )
            merged = seq_utils.merge_intervals_with_gap(windows, gap=options.gc_step)
            for start, end, gc in merged:
                if start >= end or end > seq_len:
                    continue
                results.append(
                    NegativeFeature(
                        feature_type="extreme_gc",
                        start=start,
                        end=min(end, seq_len),
                        description=f"Extreme GC window(s): GC<{options.gc_min}% or GC>{options.gc_max}%",
                        source="internal_gc",
                        score=gc,
                    )
                )

        if "homopolymer" in requested:
            hits = scan_homopolymers(full_sequence, at_run=options.homopolymer_at, gc_run=options.homopolymer_gc)
            for base, start, end in hits:
                if start >= end or end > seq_len:
                    continue
                results.append(
                    NegativeFeature(
                        feature_type="homopolymer",
                        start=start,
                        end=end,
                        description=f"Homopolymer run: {base}x{end-start}",
                        source="internal_regex",
                        score=float(end - start),
                    )
                )

        if "ambiguous" in requested:
            blocks = scan_ambiguous(full_sequence)
            for start, end in blocks:
                if start >= end or end > seq_len:
                    continue
                results.append(
                    NegativeFeature(
                        feature_type="ambiguous",
                        start=start,
                        end=end,
                        description="Ambiguous base(s) present",
                        source="internal_regex",
                    )
                )

        return results


def _scan_extreme_gc(
    sequence: str,
    *,
    window: int,
    step: int,
    gc_min: float,
    gc_max: float,
) -> list[NegativeFeature]:
    """Backward compatible helper used by legacy tests."""
    hits = scan_extreme_gc_windows(
        sequence.upper(),
        window_size=window,
        step=step,
        gc_min=gc_min,
        gc_max=gc_max,
    )
    results: list[NegativeFeature] = []
    for start, end, gc in hits:
        if start >= end:
            continue
        results.append(
            NegativeFeature(
                feature_type="extreme_gc",
                start=start,
                end=end,
                description=f"Extreme GC window(s): GC<{gc_min}% or GC>{gc_max}%",
                source="internal_gc",
                score=gc,
            )
        )
    return results
