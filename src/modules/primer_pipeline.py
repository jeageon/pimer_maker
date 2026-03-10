from __future__ import annotations

from dataclasses import dataclass, asdict
from datetime import datetime
from io import StringIO
import hashlib
import re
from typing import Any, Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import MeltingTemp

from ..config import (
    DEFAULT_GC_CLAMP_MAX,
    DEFAULT_GC_CLAMP_MIN,
    DEFAULT_GC_MAX,
    DEFAULT_GC_MIN,
    DEFAULT_MANUAL_HAIRPIN_MAX_K,
    DEFAULT_MANUAL_HAIRPIN_MIN_K,
    DEFAULT_MANUAL_SELF_DIMER_EXCLUDE_IDENTICAL_WINDOW,
    DEFAULT_MANUAL_OFFTARGET_SEED_LEN,
    DEFAULT_MANUAL_OFFTARGET_SEED_WARNING_LIMIT,
    DEFAULT_MANUAL_PAIR_DIMER_MAX_OVERLAP,
    DEFAULT_MANUAL_PAIR_DIMER_MIN_OVERLAP,
    DEFAULT_MANUAL_REQUIRE_3P_DIMER,
    DEFAULT_MANUAL_SELF_DIMER_MAX_OVERLAP,
    DEFAULT_MANUAL_SELF_DIMER_MIN_OVERLAP,
    DEFAULT_MANUAL_TM_GAP_FAIL,
    DEFAULT_INTERFERENCE_GC_MAX,
    DEFAULT_INTERFERENCE_GC_MIN,
    DEFAULT_INTERFERENCE_HOMOPOLYMER_AT,
    DEFAULT_INTERFERENCE_HOMOPOLYMER_GC,
    DEFAULT_INTERFERENCE_GC_STEP,
    DEFAULT_INTERFERENCE_GC_WINDOW,
    DEFAULT_INTERFERENCE_REPEAT_RUN,
    DEFAULT_IDEAL_PRIMER_TM_MAX,
    DEFAULT_IDEAL_PRIMER_TM_MIN,
    DEFAULT_MAX_PRODUCT_SIZE,
    DEFAULT_PRIMER_TM_GAP_MAX,
    DEFAULT_IDEAL_REPEAT_UNIT_MAX,
    DEFAULT_MIN_PRODUCT_SIZE,
    DEFAULT_PRIMER_LEN_MAX,
    DEFAULT_PRIMER_LEN_MIN,
    DEFAULT_PRIMER_TM_TARGET,
    DEFAULT_PRIMER_TM_TOLERANCE,
    PRIMER_OUTPUT_FILE_SUFFIX,
)
from ..models.data_schemas import NegativeFeature
from ..utils.seq_utils import gc_percent, scan_ambiguous, scan_extreme_gc_windows, scan_homopolymers


class PrimerPipelineError(RuntimeError):
    """Runtime error for primer design pipeline."""


@dataclass(frozen=True)
class PrimerCandidate:
    primer_id: str
    sequence: str
    start: int
    end: int
    length: int
    strand: int
    tm: float
    gc: float
    gc_clamp: int
    score: float
    is_ideal: bool


@dataclass(frozen=True)
class RejectedPrimerFeature:
    reason: str
    sequence: str
    start: int
    end: int
    length: int
    details: dict[str, Any]


@dataclass(frozen=True)
class PrimerPairBinding:
    forward_start: int
    forward_end: int
    forward_strand: int
    reverse_start: int
    reverse_end: int
    reverse_strand: int
    product_size: int
    tm_forward: float
    tm_reverse: float
    tm_gap: float
    seed_warnings: int
    dimer_risk: bool


_VALID_DNA = re.compile(r"^[ACGTNacgtn]+$")


def run_primer_pipeline(
    gb_bytes: bytes | bytearray,
    *,
    source_filename: str | None = None,
    tm_target: float = DEFAULT_PRIMER_TM_TARGET,
    tm_tolerance: float = DEFAULT_PRIMER_TM_TOLERANCE,
    gc_min: float = DEFAULT_GC_MIN,
    gc_max: float = DEFAULT_GC_MAX,
    len_min: int = DEFAULT_PRIMER_LEN_MIN,
    len_max: int = DEFAULT_PRIMER_LEN_MAX,
    repeat_run_limit: int = DEFAULT_INTERFERENCE_REPEAT_RUN,
    max_candidates: Optional[int] = None,
    product_min: int = DEFAULT_MIN_PRODUCT_SIZE,
    product_max: int = DEFAULT_MAX_PRODUCT_SIZE,
    gc_clamp_min: int = DEFAULT_GC_CLAMP_MIN,
    gc_clamp_max: int = DEFAULT_GC_CLAMP_MAX,
    ideal_tm_min: float = DEFAULT_IDEAL_PRIMER_TM_MIN,
    ideal_tm_max: float = DEFAULT_IDEAL_PRIMER_TM_MAX,
    ideal_tm_gap: float = DEFAULT_PRIMER_TM_GAP_MAX,
    ideal_repeat_unit_limit: int = DEFAULT_IDEAL_REPEAT_UNIT_MAX,
    interference_window: int = DEFAULT_INTERFERENCE_GC_WINDOW,
    interference_step: int = DEFAULT_INTERFERENCE_GC_STEP,
    interference_gc_min: float = DEFAULT_INTERFERENCE_GC_MIN,
    interference_gc_max: float = DEFAULT_INTERFERENCE_GC_MAX,
    interference_repeat_at: int = DEFAULT_INTERFERENCE_HOMOPOLYMER_AT,
    interference_repeat_gc: int = DEFAULT_INTERFERENCE_HOMOPOLYMER_GC,
    include_input_features: bool = True,
    self_dimer_exclude_identical_window: bool = DEFAULT_MANUAL_SELF_DIMER_EXCLUDE_IDENTICAL_WINDOW,
) -> dict[str, Any]:
    """Run GB upload -> interference mapping -> primer scan pipeline."""
    parse_candidates = list(_iter_decoded_texts(gb_bytes))
    record = None
    record_index = 1
    record_count = 0
    raw_seq_text = parse_candidates[0] if parse_candidates else ""

    for candidate_text in parse_candidates:
        record_match = _try_parse_record(candidate_text)
        if record_match:
            record, record_index, record_count = record_match
            raw_seq_text = candidate_text
            break

    if record is None:
        raw_sequence: str = ""
        for candidate_text in parse_candidates:
            raw_sequence = _extract_raw_sequence(candidate_text)
            if raw_sequence:
                break
        if not raw_sequence:
            raise PrimerPipelineError(
                "유효한 GenBank(.gb/.gbk), FASTA(.fasta), 또는 염기서열 텍스트를 찾지 못했습니다. "
                "빈 파일이거나 파일 인코딩/헤더(LOCUS, >) 또는 포맷이 깨진 것 같습니다."
            )
        if not _VALID_DNA.fullmatch(raw_sequence):
            raise PrimerPipelineError(
                "업로드 파일에서 A/C/G/T/N 이외의 문자만 발견되었습니다. "
                "GenBank/FASTA 형식이거나 염기서열만 포함된 텍스트인지 확인하세요."
            )
        record = SeqRecord(
            Seq(raw_sequence),
            id="raw_sequence",
            description="converted_from_raw_sequence",
        )
        record_count = 1

    sequence = str(record.seq).upper()

    if not sequence:
        raise PrimerPipelineError("GenBank sequence is empty.")
    if not _VALID_DNA.fullmatch(sequence):
        invalid = sorted(set(re.findall(r"[^ACGTN]", sequence, flags=re.IGNORECASE)))
        raise PrimerPipelineError(f"Invalid DNA bases in sequence: {', '.join(invalid)}")

    plasmid_name = _resolve_plasmid_name(record, source_filename=source_filename)

    interference = scan_interference_regions(
        record=record,
        sequence=sequence,
        window=interference_window,
        step=interference_step,
        gc_min=interference_gc_min,
        gc_max=interference_gc_max,
        repeat_run_at=interference_repeat_at,
        repeat_run_gc=interference_repeat_gc,
        repeat_run_limit=repeat_run_limit,
        include_input_features=include_input_features,
    )

    primer_candidates, rejected_primer_features = find_candidates(
        sequence=sequence,
        interference_regions=interference,
        tm_target=tm_target,
        tm_tolerance=tm_tolerance,
        gc_min=gc_min,
        gc_max=gc_max,
        len_min=len_min,
        len_max=len_max,
        gc_clamp_min=gc_clamp_min,
        gc_clamp_max=gc_clamp_max,
        ideal_tm_min=ideal_tm_min,
        ideal_tm_max=ideal_tm_max,
        ideal_tm_gap=ideal_tm_gap,
        ideal_repeat_unit_limit=ideal_repeat_unit_limit,
        repeat_run_limit=repeat_run_limit,
        max_candidates=max_candidates,
        self_dimer_exclude_identical_window=self_dimer_exclude_identical_window,
    )
    palindrome_regions = _build_palindrome_interference_regions(
        rejected_primer_features,
        seq_len=len(sequence),
    )

    run_id = hashlib.sha1(gb_bytes).hexdigest()[:10]
    annotated_record = build_annotated_record(
        record=record,
        interference_regions=interference,
        primer_candidates=primer_candidates,
        palindrome_regions=palindrome_regions,
        plasmid_name=plasmid_name,
    )
    gb_text = _to_genbank_text(annotated_record)
    safe_name = _safe_name(plasmid_name)

    metadata = {
        "run_id": run_id,
        "generated_at": datetime.utcnow().isoformat() + "Z",
        "input_id": record.id,
        "source_filename": source_filename,
        "plasmid_name": plasmid_name,
        "sequence_length": len(sequence),
        "tm_target": tm_target,
        "tm_tolerance": tm_tolerance,
        "gc_min": gc_min,
        "gc_max": gc_max,
        "len_min": len_min,
        "len_max": len_max,
        "repeat_run_limit": repeat_run_limit,
        "include_input_features": include_input_features,
        "product_size_range": [product_min, product_max],
        "interference_count": len(interference),
        "interference_by_type": _count_by_feature_type(interference),
        "rejected_primer_feature_count": len(rejected_primer_features),
        "palindrome_rejected_primer_feature_count": sum(
            1
            for feature in rejected_primer_features
            if feature.reason in {"hairpin", "self_dimer"}
        ),
        "palindrome_interference_count": len(palindrome_regions),
        "self_dimer_exclude_identical_window": self_dimer_exclude_identical_window,
        "primer_candidate_count": len(primer_candidates),
        "primer_ideal_count": sum(1 for item in primer_candidates if item.is_ideal),
        "input_record_count": record_count,
        "used_record_index": record_index,
    }

    return {
        "record_id": record.id,
        "record_name": plasmid_name,
        "filename": f"{safe_name}{PRIMER_OUTPUT_FILE_SUFFIX}",
        "gb_text": gb_text,
        "sequence": sequence,
        "metadata": metadata,
        "interference_regions": [item.model_dump() for item in interference],
        "palindrome_interference_regions": [item.model_dump() for item in palindrome_regions],
        "primer_candidates": [asdict(item) for item in primer_candidates],
        "rejected_primer_features": [asdict(item) for item in rejected_primer_features],
    }


def _resolve_plasmid_name(record: SeqRecord, source_filename: str | None = None) -> str:
    if source_filename:
        candidate = source_filename.strip()
        if candidate:
            candidate = candidate.rsplit("/", 1)[-1].rsplit("\\", 1)[-1]
            if "." in candidate:
                candidate = ".".join(candidate.split(".")[:-1]) or candidate
            candidate = re.sub(r"^\.+", "", candidate).strip()
            if candidate:
                return candidate

    record_name = (record.name or "").strip() if record.name else ""
    if record_name and record_name != ".":
        return record_name

    record_id = (record.id or "").strip() if record.id else ""
    if record_id and record_id != ".":
        return record_id

    source = (record.annotations.get("source", "") or "").strip()
    if source:
        return source

    return "plasmid"


def scan_interference_regions(
    record: SeqRecord,
    sequence: str,
    *,
    window: int,
    step: int,
    gc_min: float,
    gc_max: float,
    repeat_run_at: int,
    repeat_run_gc: int,
    repeat_run_limit: int,
    include_input_features: bool = True,
) -> list[NegativeFeature]:
    seq_len = len(sequence)
    regions: list[NegativeFeature] = []

    if include_input_features:
        for idx, feature in enumerate(record.features):
            if feature.type.lower() == "source":
                continue
            if not hasattr(feature, "location") or feature.location is None:
                continue
            for start, end, strand in _location_intervals(feature.location, seq_len):
                if start < 0:
                    start = 0
                if end > seq_len:
                    end = seq_len
                if start >= end:
                    continue
                # skip over-broad annotations that would block whole-template search
                if end - start >= seq_len:
                    continue
                regions.append(
                    NegativeFeature(
                        feature_type="existing",
                        start=start,
                        end=end,
                        description=f"existing feature in input: {feature.type}",
                        source="input_gb",
                        strand=strand,
                        attributes={"index": idx, "original_type": feature.type},
                    )
                )

    for start, end, gc in scan_extreme_gc_windows(
        sequence,
        window_size=window,
        step=step,
        gc_min=gc_min,
        gc_max=gc_max,
    ):
        if 0 <= start < end <= seq_len:
            regions.append(
                NegativeFeature(
                    feature_type="extreme_gc",
                    start=start,
                    end=end,
                    description=f"GC out of preferred range ({gc_min} - {gc_max}): {gc:.1f}%",
                    source="pipeline",
                    score=gc,
                    attributes={"gc": gc},
                )
            )

    for base, start, end in scan_homopolymers(
        sequence,
        at_run=repeat_run_at,
        gc_run=repeat_run_gc,
    ):
        if start >= end or end > seq_len:
            continue
        run_len = end - start
        if run_len >= repeat_run_limit:
            regions.append(
                NegativeFeature(
                    feature_type="homopolymer",
                    start=start,
                    end=end,
                    description=f"homopolymer run: {base}x{run_len}",
                    source="pipeline",
                    score=float(run_len),
                    attributes={"base": base, "run_len": run_len},
                )
            )

    for start, end in scan_ambiguous(sequence):
        if start >= end or end > seq_len:
            continue
        regions.append(
            NegativeFeature(
                feature_type="ambiguous",
                start=start,
                end=end,
                description="ambiguous base range",
                source="pipeline",
            )
        )

    return _merge_intervals([item for item in regions if item.end > item.start])


def find_candidates(
    sequence: str,
    interference_regions: list[NegativeFeature],
    *,
    tm_target: float,
    tm_tolerance: float,
    gc_min: float,
    gc_max: float,
    len_min: int,
    len_max: int,
    gc_clamp_min: int,
    gc_clamp_max: int,
    ideal_tm_min: float,
    ideal_tm_max: float,
    ideal_tm_gap: float,
    ideal_repeat_unit_limit: int,
    repeat_run_limit: int,
    max_candidates: Optional[int],
    self_dimer_min_overlap: int = DEFAULT_MANUAL_SELF_DIMER_MIN_OVERLAP,
    self_dimer_max_overlap: int = DEFAULT_MANUAL_SELF_DIMER_MAX_OVERLAP,
    self_dimer_exclude_identical_window: bool = DEFAULT_MANUAL_SELF_DIMER_EXCLUDE_IDENTICAL_WINDOW,
) -> tuple[list[PrimerCandidate], list[RejectedPrimerFeature]]:
    if len_min < 12:
        len_min = 12
    if len_max < len_min:
        len_max = len_min

    seq_len = len(sequence)
    forbidden = _intervals_from_features(interference_regions)
    forward_candidates: list[PrimerCandidate] = []
    rejected_features: list[RejectedPrimerFeature] = []
    seq = sequence.upper()

    if ideal_tm_min > ideal_tm_max:
        ideal_tm_min, ideal_tm_max = ideal_tm_max, ideal_tm_min

    for primer_len in range(len_min, len_max + 1):
        for start in range(0, seq_len - primer_len + 1):
            end = start + primer_len
            template_window = seq[start:end]
            if "N" in template_window:
                rejected_features.append(
                    RejectedPrimerFeature(
                        reason="contains_N",
                        sequence=template_window,
                        start=start,
                        end=end,
                        length=primer_len,
                        details={},
                    )
                )
                continue
            if _overlaps_any(start, end, forbidden):
                rejected_features.append(
                    RejectedPrimerFeature(
                        reason="overlaps_interference",
                        sequence=template_window,
                        start=start,
                        end=end,
                        length=primer_len,
                        details={},
                    )
                )
                continue
            if len(template_window) != primer_len or _has_repeat_run(template_window, repeat_run_limit):
                rejected_features.append(
                    RejectedPrimerFeature(
                        reason="repeat_run",
                        sequence=template_window,
                        start=start,
                        end=end,
                        length=primer_len,
                        details={},
                    )
                )
                continue
            if not _gc_clamp_ok(template_window, gc_clamp_min, gc_clamp_max):
                rejected_features.append(
                    RejectedPrimerFeature(
                        reason="gc_clamp",
                        sequence=template_window,
                        start=start,
                        end=end,
                        length=primer_len,
                        details={"gc_clamp": _gc_clamp_count(template_window)},
                    )
                )
                continue
            if _has_hairpin_like(template_window):
                rejected_features.append(
                    RejectedPrimerFeature(
                        reason="hairpin",
                        sequence=template_window,
                        start=start,
                        end=end,
                        length=primer_len,
                        details=_hairpin_overlap_info(template_window) or {},
                    )
                )
                continue
            self_dimer_info = _pair_dimer_overlap_info(
                template_window,
                template_window,
                min_overlap=self_dimer_min_overlap,
                max_overlap=self_dimer_max_overlap,
                full_sequence_scan=True,
                allow_identical_window_match=not self_dimer_exclude_identical_window,
            )
            if self_dimer_info:
                rejected_features.append(
                    RejectedPrimerFeature(
                        reason="self_dimer",
                        sequence=template_window,
                        start=start,
                        end=end,
                        length=primer_len,
                        details=self_dimer_info,
                    )
                )
                continue

            tm = _calc_tm(template_window)
            if tm is None:
                rejected_features.append(
                    RejectedPrimerFeature(
                        reason="tm_calc_failed",
                        sequence=template_window,
                        start=start,
                        end=end,
                        length=primer_len,
                        details={},
                    )
                )
                continue
            if not (tm_target - tm_tolerance <= tm <= tm_target + tm_tolerance):
                rejected_features.append(
                    RejectedPrimerFeature(
                        reason="tm_out_of_range",
                        sequence=template_window,
                        start=start,
                        end=end,
                        length=primer_len,
                        details={
                            "tm": tm,
                            "tm_target": tm_target,
                            "tm_tolerance": tm_tolerance,
                        },
                    )
                )
                continue

            gc = gc_percent(template_window)
            if not (gc_min <= gc <= gc_max):
                rejected_features.append(
                    RejectedPrimerFeature(
                        reason="gc_out_of_range",
                        sequence=template_window,
                        start=start,
                        end=end,
                        length=primer_len,
                        details={"gc": gc, "gc_min": gc_min, "gc_max": gc_max},
                    )
                )
                continue

            score = _score_candidate(
                tm=tm,
                tm_target=tm_target,
                gc=gc,
                gc_min=gc_min,
                gc_max=gc_max,
                length=primer_len,
                len_min=len_min,
                len_max=len_max,
            )
            is_ideal = _is_ideal_primer_candidate(
                sequence=template_window,
                tm=tm,
                gc=gc,
                tm_target=tm_target,
                gc_min=gc_min,
                gc_max=gc_max,
                len_min=len_min,
                len_max=len_max,
                gc_clamp_min=gc_clamp_min,
                gc_clamp_max=gc_clamp_max,
                repeat_run_limit=repeat_run_limit,
                ideal_tm_min=ideal_tm_min,
                ideal_tm_max=ideal_tm_max,
                ideal_tm_gap=ideal_tm_gap,
                ideal_repeat_unit_limit=ideal_repeat_unit_limit,
            )
            forward_candidates.append(
                PrimerCandidate(
                    primer_id="",
                    sequence=template_window,
                    start=start,
                    end=end,
                    length=primer_len,
                    strand=1,
                    tm=tm,
                    gc=gc,
                    gc_clamp=_gc_clamp_count(template_window),
                    score=score,
                    is_ideal=is_ideal,
                )
            )

    forward_candidates.sort(key=lambda item: item.score)
    if max_candidates is not None and max_candidates > 0:
        top = forward_candidates[:max_candidates]
    else:
        top = forward_candidates

    result: list[PrimerCandidate] = []
    for idx, cand in enumerate(top, start=1):
        result.append(
            PrimerCandidate(
                primer_id=f"F{idx:04d}",
                sequence=cand.sequence,
                start=cand.start,
                end=cand.end,
                length=cand.length,
                strand=1,
                tm=cand.tm,
                gc=cand.gc,
                gc_clamp=cand.gc_clamp,
                score=cand.score,
                is_ideal=cand.is_ideal,
            )
        )

    for idx, cand in enumerate(top, start=1):
        rc = str(Seq(cand.sequence).reverse_complement())
        reverse_is_ideal = _is_ideal_primer_candidate(
            sequence=rc,
            tm=cand.tm,
            gc=cand.gc,
            tm_target=tm_target,
            gc_min=gc_min,
            gc_max=gc_max,
            len_min=len_min,
            len_max=len_max,
            gc_clamp_min=gc_clamp_min,
            gc_clamp_max=gc_clamp_max,
            repeat_run_limit=repeat_run_limit,
            ideal_tm_min=ideal_tm_min,
            ideal_tm_max=ideal_tm_max,
            ideal_tm_gap=ideal_tm_gap,
            ideal_repeat_unit_limit=ideal_repeat_unit_limit,
        )
        result.append(
            PrimerCandidate(
                primer_id=f"R{idx:04d}",
                sequence=rc,
                start=cand.start,
                end=cand.end,
                length=cand.length,
                strand=-1,
                tm=cand.tm,
                gc=cand.gc,
                gc_clamp=cand.gc_clamp,
                score=cand.score,
                is_ideal=reverse_is_ideal,
            )
    )

    return result, rejected_features


def validate_primer_pair(
    sequence: str,
    forward_seq: str,
    reverse_seq: str,
    *,
    product_min: int,
    product_max: int,
    tm_gap_fail: float = DEFAULT_MANUAL_TM_GAP_FAIL,
    hairpin_min_k: int = DEFAULT_MANUAL_HAIRPIN_MIN_K,
    hairpin_max_k: int = DEFAULT_MANUAL_HAIRPIN_MAX_K,
    self_dimer_min_overlap: int = DEFAULT_MANUAL_SELF_DIMER_MIN_OVERLAP,
    self_dimer_max_overlap: int = DEFAULT_MANUAL_SELF_DIMER_MAX_OVERLAP,
    pair_dimer_min_overlap: int = DEFAULT_MANUAL_PAIR_DIMER_MIN_OVERLAP,
    pair_dimer_max_overlap: int = DEFAULT_MANUAL_PAIR_DIMER_MAX_OVERLAP,
    pair_dimer_require_3p: bool = DEFAULT_MANUAL_REQUIRE_3P_DIMER,
    offtarget_seed_len: int = DEFAULT_MANUAL_OFFTARGET_SEED_LEN,
    offtarget_seed_warning_limit: int = DEFAULT_MANUAL_OFFTARGET_SEED_WARNING_LIMIT,
    self_dimer_exclude_identical_window: bool = DEFAULT_MANUAL_SELF_DIMER_EXCLUDE_IDENTICAL_WINDOW,
    tm_target: float,
    tm_tolerance: float,
) -> dict[str, Any]:
    sequence = sequence.upper()
    forward = _normalize_primer(forward_seq)
    reverse = _normalize_primer(reverse_seq)
    interference_details: list[dict[str, Any]] = []

    if not forward or not reverse:
        return {
            "valid": False,
            "errors": ["both forward and reverse primer sequences are required"],
            "interference_details": [],
            "pairs": [],
        }

    if not _VALID_DNA.fullmatch(forward) or not _VALID_DNA.fullmatch(reverse):
        return {
            "valid": False,
            "errors": ["primers must contain only A/C/G/T/N characters"],
            "interference_details": [],
            "pairs": [],
        }

    forward_tm = _calc_tm(forward)
    reverse_tm = _calc_tm(reverse)
    if forward_tm is None or reverse_tm is None:
        return {
            "valid": False,
            "errors": ["failed to calculate Tm"],
            "interference_details": [],
            "pairs": [],
        }

    if tm_gap_fail > 0 and abs(forward_tm - reverse_tm) > tm_gap_fail:
        return {
            "valid": False,
            "errors": [f"Tm gap too large ({abs(forward_tm - reverse_tm):.2f}): exceeds failure limit {tm_gap_fail:.2f}"],
            "interference_details": [
                {
                    "scope": "pair",
                    "type": "tm_gap",
                    "forward_tm": forward_tm,
                    "reverse_tm": reverse_tm,
                    "gap": abs(forward_tm - reverse_tm),
                    "limit": tm_gap_fail,
                }
            ],
            "pairs": [],
        }

    interference_errors: list[str] = []
    forward_hairpin = _hairpin_overlap_info(forward, min_k=hairpin_min_k, max_k=hairpin_max_k)
    if forward_hairpin:
        interference_errors.append("forward primer is likely to form hairpin")
        interference_details.append(
            {
                "scope": "single_primer",
                "primer": "forward",
                "type": "hairpin",
                "overlap_len": forward_hairpin["overlap_len"],
                "primer_tail": forward_hairpin["tail"],
                "internal_match": forward_hairpin["match_seq"],
                "match_index": forward_hairpin["match_start"],
            }
        )

    reverse_hairpin = _hairpin_overlap_info(reverse, min_k=hairpin_min_k, max_k=hairpin_max_k)
    if reverse_hairpin:
        interference_errors.append("reverse primer is likely to form hairpin")
        interference_details.append(
            {
                "scope": "single_primer",
                "primer": "reverse",
                "type": "hairpin",
                "overlap_len": reverse_hairpin["overlap_len"],
                "primer_tail": reverse_hairpin["tail"],
                "internal_match": reverse_hairpin["match_seq"],
                "match_index": reverse_hairpin["match_start"],
            }
        )

    forward_self_dimer = _pair_dimer_overlap_info(
        forward,
        forward,
        min_overlap=self_dimer_min_overlap,
        max_overlap=self_dimer_max_overlap,
        full_sequence_scan=True,
        allow_identical_window_match=not self_dimer_exclude_identical_window,
    )
    if forward_self_dimer:
        interference_errors.append("forward primer is likely to form self-dimer")
        interference_details.append(
            {
                "scope": "single_primer",
                "primer": "forward",
                "type": "self_dimer",
                "overlap_len": forward_self_dimer["overlap_len"],
                "primer_tail": forward_self_dimer["a_tail"],
                "complement_match": forward_self_dimer["match_seq"],
                "paired_index": forward_self_dimer["matched_index_in_b"],
                "pairing": forward_self_dimer["pairing"],
            }
        )

    reverse_self_dimer = _pair_dimer_overlap_info(
        reverse,
        reverse,
        min_overlap=self_dimer_min_overlap,
        max_overlap=self_dimer_max_overlap,
        full_sequence_scan=True,
        allow_identical_window_match=not self_dimer_exclude_identical_window,
    )
    if reverse_self_dimer:
        interference_errors.append("reverse primer is likely to form self-dimer")
        interference_details.append(
            {
                "scope": "single_primer",
                "primer": "reverse",
                "type": "self_dimer",
                "overlap_len": reverse_self_dimer["overlap_len"],
                "primer_tail": reverse_self_dimer["a_tail"],
                "complement_match": reverse_self_dimer["match_seq"],
                "paired_index": reverse_self_dimer["matched_index_in_b"],
                "pairing": reverse_self_dimer["pairing"],
            }
        )

    if interference_errors:
        return {
            "valid": False,
            "errors": interference_errors,
            "interference_details": list(interference_details),
            "pairs": [],
        }

    f_hits = _find_bindings(sequence, forward)
    r_hits = _find_bindings(sequence, reverse)
    if not f_hits:
        return {
            "valid": False,
            "errors": ["forward primer does not match template"],
            "interference_details": [
                {
                    "scope": "pair",
                    "type": "binding_miss",
                    "primer": "forward",
                    "message": "forward primer has no perfect match in template",
                }
            ],
            "pairs": [],
        }
    if not r_hits:
        return {
            "valid": False,
            "errors": ["reverse primer does not match template"],
            "interference_details": [
                {
                    "scope": "pair",
                    "type": "binding_miss",
                    "primer": "reverse",
                    "message": "reverse primer has no perfect match in template",
                }
            ],
            "pairs": [],
        }

    seed_len = max(1, int(offtarget_seed_len))
    warning_limit = int(offtarget_seed_warning_limit)
    product_size_samples: list[int] = []
    off_target_samples: list[dict[str, Any]] = []
    pair_dimer_samples: list[dict[str, Any]] = []

    pairs: list[PrimerPairBinding] = []
    off_target_pairs = 0
    total_combinations = len(f_hits) * len(r_hits)
    skipped_same_strand = 0
    skipped_size = 0
    skipped_dimer = 0
    for f in f_hits:
        for r in r_hits:
            if f[2] == r[2]:
                skipped_same_strand += 1
                continue
            forward_3 = _three_prime_pos(f[0], f[1], f[2])
            reverse_3 = _three_prime_pos(r[0], r[1], r[2])
            product_size = abs(reverse_3 - forward_3) + 1
            if product_size < product_min or product_size > product_max:
                skipped_size += 1
                if len(product_size_samples) < 10:
                    product_size_samples.append(product_size)
                continue

            dimer_risk = _pair_dimer_risk(
                forward,
                reverse,
                min_overlap=pair_dimer_min_overlap,
                max_overlap=pair_dimer_max_overlap,
                require_3p=pair_dimer_require_3p,
            )
            if dimer_risk:
                skipped_dimer += 1
                dimer_match = _pair_dimer_overlap_info(
                    forward,
                    reverse,
                    min_overlap=pair_dimer_min_overlap,
                    max_overlap=pair_dimer_max_overlap,
                    require_3p=pair_dimer_require_3p,
                )
                if dimer_match and len(pair_dimer_samples) < 10:
                    pair_dimer_samples.append(
                        {
                            "scope": "pair",
                            "type": "cross_dimer",
                            "forward_pos": f[0],
                            "reverse_pos": r[0],
                            "forward_strand": f[2],
                            "reverse_strand": r[2],
                            "overlap_len": dimer_match["overlap_len"],
                            "a_tail": dimer_match["a_tail"],
                            "b_tail": dimer_match["b_tail"],
                            "match_seq": dimer_match["match_seq"],
                            "pairing": dimer_match["pairing"],
                        }
                    )
                continue

            seed_warnings = max(_seed_offtarget_count(sequence, forward[-seed_len:]) - 1, 0) + max(
                _seed_offtarget_count(sequence, reverse[-seed_len:]) - 1,
                0,
            )
            if warning_limit > 0 and seed_warnings >= warning_limit:
                off_target_pairs += 1
                if len(off_target_samples) < 10:
                    off_target_samples.append(
                        {
                            "seed_len": seed_len,
                            "forward_seed": forward[-seed_len:],
                            "reverse_seed": reverse[-seed_len:],
                            "forward_seed_hits": _seed_offtarget_count(sequence, forward[-seed_len:]),
                            "reverse_seed_hits": _seed_offtarget_count(sequence, reverse[-seed_len:]),
                        }
                    )
            pairs.append(
                PrimerPairBinding(
                    forward_start=f[0],
                    forward_end=f[1],
                    forward_strand=f[2],
                    reverse_start=r[0],
                    reverse_end=r[1],
                    reverse_strand=r[2],
                    product_size=product_size,
                    tm_forward=forward_tm,
                    tm_reverse=reverse_tm,
                    tm_gap=abs(forward_tm - reverse_tm),
                    seed_warnings=seed_warnings,
                    dimer_risk=dimer_risk,
                )
            )

    pairs.sort(key=lambda item: (item.tm_gap, item.product_size, item.seed_warnings))
    opposite_pairs_total = total_combinations - skipped_same_strand
    summary_messages = [
        f"총 후보 결합쌍: {total_combinations}",
        f"반대가닥 쌍: {opposite_pairs_total}",
        f"same-strand 제외: {skipped_same_strand}",
        f"산물 크기 범위 제외: {skipped_size}",
        f"이량체 필터 제외: {skipped_dimer}",
        f"오프타겟 잠재성 표시: {off_target_pairs}",
    ]
    if warning_limit <= 0:
        summary_messages.append("off-target 경고 임계값 0: 경고 집계 비활성화")

    if product_size_samples:
        summary_messages.append(f"예시(크기 제외): {product_size_samples[:10]}")
    if off_target_samples:
        summary_messages.append(f"예시(오프타겟 시드): {off_target_samples[:10]}")

    if not pairs:
        interference_summary = list(interference_details)
        interference_summary.extend(pair_dimer_samples)
        interference_summary.extend(off_target_samples)
        failure_reasons: list[str] = []
        if opposite_pairs_total <= 0:
            failure_reasons.append("forward/reverse binding sites are on the same strand; reverse primer may need reverse-complement input")
        if skipped_size and skipped_size == opposite_pairs_total:
            failure_reasons.append(f"all opposite-strand pairs were outside product range {product_min}-{product_max}")
        elif skipped_size > 0:
            failure_reasons.append(f"some opposite-strand pairs were filtered by product size range {product_min}-{product_max}")
        if skipped_dimer and skipped_dimer == opposite_pairs_total:
            failure_reasons.append("all opposite-strand pairs were rejected by dimer risk filter")
        elif skipped_dimer > 0:
            failure_reasons.append("some opposite-strand pairs were rejected by dimer risk filter")
        if off_target_pairs > 0 and opposite_pairs_total > 0:
            if off_target_pairs >= opposite_pairs_total:
                failure_reasons.append("all candidate pairs show high off-target seed-seed potential")
            else:
                failure_reasons.append(f"{off_target_pairs} candidate pair(s) show high off-target seed-seed potential")
        if not failure_reasons:
            failure_reasons.append("no valid pair found after current validation filters")
        return {
            "valid": False,
            "summary_messages": summary_messages,
            "filter_summary": {
                "total_combinations": total_combinations,
                "opposite_pairs_total": opposite_pairs_total,
                "skipped_same_strand": skipped_same_strand,
                "skipped_size": skipped_size,
                "skipped_dimer": skipped_dimer,
                "off_target_pairs": off_target_pairs,
                "seed_len": seed_len,
                "warning_limit": warning_limit,
                "product_size_samples": product_size_samples,
                "off_target_samples": off_target_samples,
            },
            "interference_details": interference_summary,
            "errors": ["no valid product found under configured range", *failure_reasons],
            "pairs": [],
        }

    warnings: list[str] = []
    if abs(forward_tm - reverse_tm) > tm_tolerance:
        warnings.append(f"Tm gap {abs(forward_tm - reverse_tm):.2f} exceeds tolerance {tm_tolerance}")
    if abs(forward_tm - tm_target) > tm_tolerance:
        warnings.append(f"forward_tm_outside_target_{tm_target:.1f}")
    if abs(reverse_tm - tm_target) > tm_tolerance:
        warnings.append(f"reverse_tm_outside_target_{tm_target:.1f}")

    return {
        "valid": True,
        "interference_details": list(interference_details) + list(off_target_samples),
        "summary_messages": summary_messages,
        "filter_summary": {
            "total_combinations": total_combinations,
            "opposite_pairs_total": opposite_pairs_total,
            "skipped_same_strand": skipped_same_strand,
            "skipped_size": skipped_size,
            "skipped_dimer": skipped_dimer,
            "off_target_pairs": off_target_pairs,
            "seed_len": seed_len,
            "warning_limit": warning_limit,
            "product_size_samples": product_size_samples,
            "off_target_samples": off_target_samples,
        },
        "tm_forward": forward_tm,
        "tm_reverse": reverse_tm,
        "tm_gap": abs(forward_tm - reverse_tm),
        "warnings": warnings,
        "pairs": [
            {
                "forward_start": item.forward_start,
                "forward_end": item.forward_end,
                "forward_strand": item.forward_strand,
                "reverse_start": item.reverse_start,
                "reverse_end": item.reverse_end,
                "reverse_strand": item.reverse_strand,
                "product_size": item.product_size,
                "tm_gap": item.tm_gap,
                "seed_warnings": item.seed_warnings,
                "dimer_risk": item.dimer_risk,
            }
            for item in pairs[:50]
        ],
    }


def build_annotated_record(
    record: SeqRecord,
    interference_regions: list[NegativeFeature],
    primer_candidates: list[PrimerCandidate],
    palindrome_regions: list[NegativeFeature],
    *,
    plasmid_name: str,
) -> SeqRecord:
    output = SeqRecord(record.seq, id=record.id, name=record.name, description=record.description)
    output.annotations = dict(record.annotations)
    for item in record.features:
        output.features.append(item)

    for feature in sorted(interference_regions, key=lambda item: (item.start, item.end)):
        output.features.append(
            SeqFeature(
                FeatureLocation(
                    max(0, feature.start),
                    min(len(record.seq), feature.end),
                    strand=feature.strand or 0,
                ),
                type="misc_feature",
                qualifiers={
                    "label": [f"interference_{feature.feature_type}"],
                    "note": [feature.description],
                    "source": [feature.source],
                    "reason": [feature.feature_type],
                    "ApEinfo_fwdcolor": ["#ff9999"],
                    "ApEinfo_revcolor": ["#ff6666"],
                },
            )
        )

    for feature in sorted(palindrome_regions, key=lambda item: (item.start, item.end)):
        reason = str(feature.attributes.get("reason") or feature.feature_type).replace("palindrome_", "")
        output.features.append(
            SeqFeature(
                FeatureLocation(
                    max(0, feature.start),
                    min(len(record.seq), feature.end),
                    strand=feature.strand or 0,
                ),
                type="misc_feature",
                qualifiers={
                    "label": [f"excluded_{reason}"],
                    "note": [feature.description],
                    "source": [feature.source],
                    "reason": [reason],
                    "core_feature_type": [feature.feature_type],
                    "ApEinfo_fwdcolor": ["#ffb347"],
                    "ApEinfo_revcolor": ["#ffb347"],
                },
            )
        )

    for primer in primer_candidates:
        color = "#85dae9" if primer.strand == 1 else "#b3df8f"
        direction = "F" if primer.strand == 1 else "R"
        primer_label = _format_primer_label(
            plasmid_name=plasmid_name,
            primer=primer,
        )
        output.features.append(
            SeqFeature(
                FeatureLocation(primer.start, primer.end, strand=primer.strand),
                type="primer_bind",
                qualifiers={
                    "label": [primer_label],
                    "note": [
                        f"sequence={primer.sequence}",
                        f"Tm={primer.tm:.2f}",
                        f"GC={primer.gc:.2f}",
                        f"len={primer.length}",
                        f"primer_id={primer.primer_id}",
                        f"orientation={direction}",
                        f"strand={primer.strand}",
                    ],
                    "primer_strand": [str(primer.strand)],
                    "orientation": [direction],
                    "strand": [str(primer.strand)],
                    "gc_clamp": [str(primer.gc_clamp)],
                    "primer_name": [primer_label],
                    "ApEinfo_fwdcolor": [color],
                    "ApEinfo_revcolor": [color],
                },
            )
        )

    return output


def _iter_decoded_texts(gb_bytes: bytes | bytearray) -> list[str]:
    bytes_value = bytes(gb_bytes)
    encodings = [
        "utf-8-sig",
        "latin1",
        "utf-8",
        "utf-16",
        "utf-16-le",
        "utf-16-be",
    ]
    decoded: list[str] = []
    seen: set[str] = set()
    for encoding in encodings:
        try:
            text = bytes_value.decode(encoding)
        except UnicodeDecodeError:
            continue
        text = text.lstrip("\ufeff")
        if text not in seen:
            seen.add(text)
            decoded.append(text)
    if not decoded:
        decoded.append(bytes_value.decode("utf-8", errors="replace"))
    return decoded


def _try_parse_record(raw_seq_text: str) -> tuple[SeqRecord, int, int] | None:
    if not raw_seq_text.strip():
        return None
    try:
        genbank_records = list(SeqIO.parse(StringIO(raw_seq_text), "genbank"))
        if genbank_records:
            return genbank_records[0], 1, len(genbank_records)
    except Exception:
        pass

    try:
        fasta_records = list(SeqIO.parse(StringIO(raw_seq_text), "fasta"))
        if fasta_records:
            return fasta_records[0], 1, len(fasta_records)
    except Exception:
        pass
    return None


def _extract_raw_sequence(raw_seq_text: str) -> str:
    if not raw_seq_text:
        return ""

    text = raw_seq_text.replace("\r", "\n")
    if "ORIGIN" in text:
        payload = text.split("ORIGIN", 1)[1]
        payload = payload.split("//", 1)[0]
        return re.sub(r"[^ACGTNacgtn]", "", payload).upper()

    return _extract_sequence_text(raw_seq_text)


def _read_single_genbank_record(raw_seq_text: str) -> tuple[SeqRecord, int, int]:
    genbank_records = list(SeqIO.parse(StringIO(raw_seq_text), "genbank"))
    if genbank_records:
        return genbank_records[0], 1, len(genbank_records)

    fasta_records = list(SeqIO.parse(StringIO(raw_seq_text), "fasta"))
    if fasta_records:
        return fasta_records[0], 1, len(fasta_records)

    raw_sequence = _extract_sequence_text(raw_seq_text)
    if raw_sequence:
        if not _VALID_DNA.fullmatch(raw_sequence):
            raise PrimerPipelineError(
                "업로드 파일에서 A/C/G/T/N 이외의 문자만 발견되었습니다. "
                "GenBank/FASTA 형식이거나 염기 서열만 포함된 텍스트인지 확인하세요."
            )
        return (
            SeqRecord(
                Seq(raw_sequence),
                id="raw_sequence",
                description="converted_from_raw_sequence",
            ),
            1,
            1,
        )

    raise PrimerPipelineError(
        "유효한 GenBank(.gb/.gbk), FASTA(.fasta), 또는 염기서열 텍스트를 찾지 못했습니다. "
        "빈 파일이거나 파일 인코딩/헤더(LOCUS, >) 또는 포맷이 깨진 것 같습니다."
    )


def _extract_sequence_text(raw_seq_text: str) -> str:
    lines = [line.strip() for line in raw_seq_text.splitlines() if line.strip()]
    if not lines:
        return ""

    # Remove FASTA-style and annotation-like header lines
    body = []
    for line in lines:
        if line.startswith(">"):
            continue
        if line.upper().startswith("LOCUS"):
            continue
        if line.startswith("//"):
            continue
        body.append(line)

    compact = "".join(body)
    compact = re.sub(r"[^ACGTNacgtn]", "", compact)
    return compact.upper()

def _decode_gb_input(gb_bytes: bytes | bytearray) -> str:
    return _iter_decoded_texts(gb_bytes)[0]


def _to_genbank_text(record: SeqRecord) -> str:
    output = StringIO()
    SeqIO.write(record, output, "genbank")
    return output.getvalue()


def _format_primer_label(plasmid_name: str, primer: PrimerCandidate) -> str:
    safe_plasmid = re.sub(r"[^A-Za-z0-9._-]", "_", plasmid_name)
    safe_plasmid = (safe_plasmid or "plasmid").strip("._-")
    if primer.strand == 1:
        end_bp_1based = primer.end
    else:
        end_bp_1based = primer.start + 1
    tm_value = f"{primer.tm:.2f}".rstrip("0").rstrip(".")
    direction = "F" if primer.strand == 1 else "R"
    ideal_mark = "*" if primer.is_ideal else ""
    return f"{safe_plasmid}_{end_bp_1based}_{tm_value}_{direction}{ideal_mark}"


def _build_palindrome_interference_regions(
    rejected_primer_features: list[RejectedPrimerFeature],
    *,
    seq_len: int,
) -> list[NegativeFeature]:
    regions: list[NegativeFeature] = []
    for feature in rejected_primer_features:
        if feature.reason == "hairpin":
            regions.extend(_hairpin_core_regions(feature, seq_len=seq_len))
        elif feature.reason == "self_dimer":
            regions.extend(_self_dimer_core_regions(feature, seq_len=seq_len))
    return _merge_intervals([item for item in regions if item.end > item.start])


def _hairpin_core_regions(
    feature: RejectedPrimerFeature,
    *,
    seq_len: int,
) -> list[NegativeFeature]:
    details = feature.details or {}
    overlap_len = _safe_int(details.get("overlap_len"))
    match_start = _safe_int(details.get("match_start"), default=-1)
    match_end = _safe_int(details.get("match_end"), default=-1)
    regions: list[NegativeFeature] = []

    if overlap_len > 0:
        tail_start = max(0, feature.length - overlap_len)
        tail_end = feature.length
        tail_region = _core_region_from_local_span(
            feature,
            local_start=tail_start,
            local_end=tail_end,
            seq_len=seq_len,
            feature_type="palindrome_hairpin",
            description="hairpin-forming core region (3' tail)",
            attributes={"reason": "hairpin", "overlap_len": overlap_len, "core_role": "tail"},
        )
        if tail_region is not None:
            regions.append(tail_region)

    if match_start >= 0 and match_end > match_start:
        match_region = _core_region_from_local_span(
            feature,
            local_start=match_start,
            local_end=match_end,
            seq_len=seq_len,
            feature_type="palindrome_hairpin",
            description="hairpin-forming core region (internal match)",
            attributes={"reason": "hairpin", "overlap_len": overlap_len, "core_role": "internal"},
        )
        if match_region is not None:
            regions.append(match_region)

    if regions:
        return regions

    fallback = _core_region_from_local_span(
        feature,
        local_start=0,
        local_end=feature.length,
        seq_len=seq_len,
        feature_type="palindrome_hairpin",
        description="hairpin-forming primer region",
        attributes={"reason": "hairpin"},
    )
    return [fallback] if fallback is not None else []


def _self_dimer_core_regions(
    feature: RejectedPrimerFeature,
    *,
    seq_len: int,
) -> list[NegativeFeature]:
    details = feature.details or {}
    overlap_len = _safe_int(details.get("overlap_len"))
    paired_index_in_a = _safe_int(details.get("paired_index_in_a"), default=-1)
    paired_index_in_b = _safe_int(details.get("paired_index_in_b"), default=-1)
    regions: list[NegativeFeature] = []

    if overlap_len > 0:
        for local_start, core_role in (
            (paired_index_in_a, "segment_a"),
            (paired_index_in_b, "segment_b"),
        ):
            if local_start < 0:
                continue
            core_region = _core_region_from_local_span(
                feature,
                local_start=local_start,
                local_end=local_start + overlap_len,
                seq_len=seq_len,
                feature_type="palindrome_self_dimer",
                description="self-dimer-forming core region",
                attributes={
                    "reason": "self_dimer",
                    "overlap_len": overlap_len,
                    "core_role": core_role,
                },
            )
            if core_region is not None:
                regions.append(core_region)

    if regions:
        return regions

    fallback = _core_region_from_local_span(
        feature,
        local_start=0,
        local_end=feature.length,
        seq_len=seq_len,
        feature_type="palindrome_self_dimer",
        description="self-dimer-forming primer region",
        attributes={"reason": "self_dimer"},
    )
    return [fallback] if fallback is not None else []


def _core_region_from_local_span(
    feature: RejectedPrimerFeature,
    *,
    local_start: int,
    local_end: int,
    seq_len: int,
    feature_type: str,
    description: str,
    attributes: dict[str, Any],
) -> Optional[NegativeFeature]:
    start = max(0, feature.start + local_start)
    end = min(seq_len, feature.start + local_end)
    if end <= start:
        return None
    return NegativeFeature(
        feature_type=feature_type,
        start=start,
        end=end,
        description=description,
        source="primer_pipeline",
        strand=0,
        attributes=attributes,
    )


def _safe_int(value: object, *, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def _safe_name(value: str, max_len: int = 40) -> str:
    safe = re.sub(r"[^A-Za-z0-9._-]", "_", value)
    return (safe or "primer").strip("._-")[:max_len]


def _count_by_feature_type(features: list[NegativeFeature]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for item in features:
        counts[item.feature_type] = counts.get(item.feature_type, 0) + 1
    return counts


def _merge_intervals(features: list[NegativeFeature]) -> list[NegativeFeature]:
    if not features:
        return []
    merged: list[NegativeFeature] = []
    sorted_features = sorted(features, key=lambda item: (item.feature_type, item.start, item.end))

    for feature in sorted_features:
        if not merged:
            merged.append(feature)
            continue

        prev = merged[-1]
        if feature.feature_type == prev.feature_type and feature.start <= prev.end:
            merged[-1] = NegativeFeature(
                feature_type=prev.feature_type,
                start=min(prev.start, feature.start),
                end=max(prev.end, feature.end),
                description=f"{prev.description}; {feature.description}",
                source=prev.source,
                score=_max_score(prev.score, feature.score),
                strand=prev.strand if prev.strand is not None else feature.strand,
                attributes={**prev.attributes, **feature.attributes},
            )
        else:
            merged.append(feature)
    return merged


def _max_score(a: Optional[float], b: Optional[float]) -> Optional[float]:
    if a is None and b is None:
        return None
    if a is None:
        return b
    if b is None:
        return a
    return max(a, b)


def _intervals_from_features(features: list[NegativeFeature]) -> list[tuple[int, int]]:
    return [(item.start, item.end) for item in features]


def _location_intervals(location: Any, seq_len: int) -> list[tuple[int, int, int]]:
    """Split a feature location into one or more [start, end) genomic intervals.

    Handles simple locations, wrapped locations, and CompoundLocation from circular
    annotations.
    """
    if isinstance(location, CompoundLocation):
        intervals: list[tuple[int, int, int]] = []
        for part in location.parts:
            intervals.extend(_location_intervals(part, seq_len))
        return intervals

    strand = int(location.strand or 0)
    start = int(location.start)
    end = int(location.end)
    if start < 0 and end < 0:
        return []

    if start < 0:
        start = 0
    if end < 0:
        end = 0
    if start > end:
        return [
            (start, seq_len, strand),
            (0, end, strand),
        ]
    if start == end:
        return [(start, end, strand)]
    return [(start, end, strand)]


def _overlaps_any(start: int, end: int, intervals: list[tuple[int, int]]) -> bool:
    for a, b in intervals:
        if not (end <= a or start >= b):
            return True
    return False


def _has_repeat_run(sequence: str, max_run_limit: int) -> bool:
    if max_run_limit <= 1:
        return False
    for base in "ACGT":
        if base * max_run_limit in sequence:
            return True
    return False


def _has_tandem_repeat_unit(
    sequence: str,
    *,
    unit_size: int = 2,
    repeat_limit: int = 4,
) -> bool:
    if unit_size < 2 or repeat_limit <= 2:
        return False
    seq = sequence.upper()
    if len(seq) < unit_size * repeat_limit:
        return False

    max_start = len(seq) - unit_size + 1
    for start in range(max_start):
        unit = seq[start : start + unit_size]
        if len(unit) < unit_size:
            continue
        repeat_count = 1
        cursor = start + unit_size
        while cursor + unit_size <= len(seq) and seq[cursor : cursor + unit_size] == unit:
            repeat_count += 1
            if repeat_count >= repeat_limit:
                return True
            cursor += unit_size
    return False


def _is_ideal_primer_candidate(
    *,
    sequence: str,
    tm: float,
    gc: float,
    tm_target: float,
    gc_min: float,
    gc_max: float,
    len_min: int,
    len_max: int,
    gc_clamp_min: int,
    gc_clamp_max: int,
    repeat_run_limit: int,
    ideal_tm_min: float,
    ideal_tm_max: float,
    ideal_tm_gap: float,
    ideal_repeat_unit_limit: int,
) -> bool:
    _ = (
        tm,
        gc,
        tm_target,
        gc_min,
        gc_max,
        len_min,
        len_max,
        gc_clamp_min,
        gc_clamp_max,
        repeat_run_limit,
        ideal_tm_min,
        ideal_tm_max,
        ideal_tm_gap,
        ideal_repeat_unit_limit,
    )
    if not sequence:
        return False
    return sequence[-1].upper() in {"G", "C"}


def _gc_clamp_count(sequence: str) -> int:
    return sum(1 for c in sequence[-5:] if c in "GC")


def _gc_clamp_ok(sequence: str, min_count: int, max_count: int) -> bool:
    clamp = _gc_clamp_count(sequence)
    if not (min_count <= clamp <= max_count):
        return False
    trailing = 0
    for c in reversed(sequence):
        if c in "GC":
            trailing += 1
        else:
            break
    return trailing <= 2


def _has_hairpin_like(sequence: str, *, min_k: int = 4, max_k: int = 7) -> bool:
    return _hairpin_overlap_info(sequence, min_k=min_k, max_k=max_k) is not None


def _hairpin_overlap_info(
    sequence: str,
    *,
    min_k: int = 4,
    max_k: int = 7,
) -> Optional[dict[str, Any]]:
    seq = sequence.upper()
    if min_k < 1:
        min_k = 1
    if max_k < min_k:
        max_k = min_k
    for k in range(min_k, min(len(seq), max_k) + 1):
        tail = seq[-k:]
        rc = str(Seq(tail).reverse_complement())
        match_pos = seq[:-1].find(rc)
        if match_pos >= 0:
            return {
                "overlap_len": k,
                "tail": tail,
                "match_seq": rc,
                "match_start": match_pos,
                "match_end": match_pos + k,
            }
    return None


def _calc_tm(sequence: str) -> Optional[float]:
    if not sequence:
        return None
    try:
        return float(MeltingTemp.Tm_NN(Seq(sequence)))
    except Exception:
        return None


def _score_candidate(
    tm: float,
    tm_target: float,
    gc: float,
    gc_min: float,
    gc_max: float,
    length: int,
    len_min: int,
    len_max: int,
) -> float:
    target_gc = (gc_min + gc_max) / 2.0
    target_len = (len_min + len_max) / 2.0
    return abs(tm - tm_target) + abs(gc - target_gc) * 0.2 + abs(length - target_len) * 0.1


def _normalize_primer(value: str) -> str:
    return re.sub(r"\s+", "", (value or "").upper())


def _find_bindings(sequence: str, primer: str) -> list[tuple[int, int, int]]:
    primer_u = _normalize_primer(primer)
    primer_len = len(primer_u)
    if primer_len < 1 or primer_len > len(sequence):
        return []
    rc = str(Seq(primer_u).reverse_complement())

    matches: list[tuple[int, int, int]] = []
    seen: set[tuple[int, int, int]] = set()
    for start in range(0, len(sequence) - primer_len + 1):
        end = start + primer_len
        window = sequence[start:end]
        if window == primer_u:
            hit = (start, end, 1)
            if hit not in seen:
                seen.add(hit)
                matches.append(hit)
        if window == rc:
            hit = (start, end, -1)
            if hit not in seen:
                seen.add(hit)
                matches.append(hit)
    return matches


def _three_prime_pos(start: int, end: int, strand: int) -> int:
    return end - 1 if strand == 1 else start


def _seed_offtarget_count(sequence: str, seed: str) -> int:
    if len(seed) < 4:
        return 0
    return sum(1 for start in range(0, len(sequence) - len(seed) + 1) if sequence[start : start + len(seed)] == seed)


def _self_dimer_risk(sequence: str, *, min_overlap: int = 4, max_overlap: int = 7) -> bool:
    return _pair_dimer_overlap_info(
        sequence,
        sequence,
        min_overlap=min_overlap,
        max_overlap=max_overlap,
        full_sequence_scan=True,
        allow_identical_window_match=not DEFAULT_MANUAL_SELF_DIMER_EXCLUDE_IDENTICAL_WINDOW,
    ) is not None


def _pair_dimer_risk(
    seq_a: str,
    seq_b: str,
    *,
    min_overlap: int = 4,
    max_overlap: int = 7,
    require_3p: bool = False,
) -> bool:
    return _pair_dimer_overlap_info(
        seq_a,
        seq_b,
        min_overlap=min_overlap,
        max_overlap=max_overlap,
        require_3p=require_3p,
    ) is not None


def _pair_dimer_overlap_info(
    seq_a: str,
    seq_b: str,
    *,
    min_overlap: int = 4,
    max_overlap: int = 7,
    require_3p: bool = False,
    full_sequence_scan: bool = False,
    allow_identical_window_match: bool = True,
) -> Optional[dict[str, Any]]:
    a = seq_a.upper()
    b = seq_b.upper()
    min_overlap = max(1, min_overlap)
    if max_overlap < min_overlap:
        max_overlap = min_overlap
    max_len = min(max_overlap, len(a), len(b))
    if max_len <= 0:
        return None

    check_a_start = 0
    check_b_start = 0
    if require_3p:
        check_a_start = max(0, len(b) - 2 * max_len)
        check_b_start = max(0, len(a) - 2 * max_len)
    for k in range(min_overlap, max_len + 1):
        if full_sequence_scan:
            for start_a in range(len(a) - k + 1):
                a_seg = a[start_a : start_a + k]
                a_seg_rc = str(Seq(a_seg).reverse_complement())
                start_b = b.find(a_seg_rc)
                while start_b >= 0:
                    if (
                        not allow_identical_window_match
                        and start_b == start_a
                        and start_b + k <= len(b)
                        and a[start_a : start_a + k] == b[start_b : start_b + k]
                    ):
                        start_b = b.find(a_seg_rc, start_b + 1)
                        continue
                    if not require_3p or (start_a + k >= len(a) - max_len and start_b + k >= len(b) - max_len):
                        return {
                            "overlap_len": k,
                            "a_tail": a_seg,
                            "b_tail": b[start_b : start_b + k],
                            "match_seq": a_seg_rc,
                            "paired_index_in_a": start_a,
                            "paired_index_in_b": start_b,
                            "match_in": "b",
                            "pairing": "a_seg_to_b",
                        }
                    next_pos = b.find(a_seg_rc, start_b + 1)
                    start_b = next_pos

            for start_b in range(len(b) - k + 1):
                b_seg = b[start_b : start_b + k]
                b_seg_rc = str(Seq(b_seg).reverse_complement())
                start_a = a.find(b_seg_rc)
                while start_a >= 0:
                    if (
                        not allow_identical_window_match
                        and start_a == start_b
                        and start_a + k <= len(a)
                        and a[start_a : start_a + k] == b[start_b : start_b + k]
                    ):
                        start_a = a.find(b_seg_rc, start_a + 1)
                        continue
                    if not require_3p or (start_a + k >= len(a) - max_len and start_b + k >= len(b) - max_len):
                        return {
                            "overlap_len": k,
                            "a_tail": a[start_a : start_a + k],
                            "b_tail": b_seg,
                            "match_seq": b_seg_rc,
                            "paired_index_in_a": start_a,
                            "paired_index_in_b": start_b,
                            "match_in": "a",
                            "pairing": "b_seg_to_a",
                        }
                    next_pos = a.find(b_seg_rc, start_a + 1)
                    start_a = next_pos
            continue

        a_tail = a[-k:]
        b_tail = b[-k:]
        a_tail_rc = str(Seq(a_tail).reverse_complement())
        pos_in_b = b.find(a_tail_rc)
        if pos_in_b >= 0:
            if not require_3p or pos_in_b >= check_a_start:
                return {
                    "overlap_len": k,
                    "a_tail": a_tail,
                    "b_tail": b_tail,
                    "match_seq": a_tail_rc,
                    "paired_index_in_a": len(a) - k,
                    "paired_index_in_b": pos_in_b,
                    "match_in": "b",
                    "pairing": "a_tail_to_b",
                }

        b_tail_rc = str(Seq(b_tail).reverse_complement())
        pos_in_a = a.find(b_tail_rc)
        if pos_in_a >= 0:
            if not require_3p or pos_in_a >= check_b_start:
                return {
                    "overlap_len": k,
                    "a_tail": a_tail,
                    "b_tail": b_tail,
                    "match_seq": b_tail_rc,
                    "paired_index_in_a": pos_in_a,
                    "paired_index_in_b": len(b) - k,
                    "match_in": "a",
                    "pairing": "b_tail_to_a",
                }
    return None
