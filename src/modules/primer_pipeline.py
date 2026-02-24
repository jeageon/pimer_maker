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

    primer_candidates = find_candidates(
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
    )

    run_id = hashlib.sha1(gb_bytes).hexdigest()[:10]
    annotated_record = build_annotated_record(
        record=record,
        interference_regions=interference,
        primer_candidates=primer_candidates,
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
        "primer_candidates": [asdict(item) for item in primer_candidates],
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
) -> list[PrimerCandidate]:
    if len_min < 12:
        len_min = 12
    if len_max < len_min:
        len_max = len_min

    seq_len = len(sequence)
    forbidden = _intervals_from_features(interference_regions)
    forward_candidates: list[PrimerCandidate] = []
    seq = sequence.upper()

    if ideal_tm_min > ideal_tm_max:
        ideal_tm_min, ideal_tm_max = ideal_tm_max, ideal_tm_min

    for primer_len in range(len_min, len_max + 1):
        for start in range(0, seq_len - primer_len + 1):
            end = start + primer_len
            template_window = seq[start:end]
            if "N" in template_window:
                continue
            if _overlaps_any(start, end, forbidden):
                continue
            if len(template_window) != primer_len or _has_repeat_run(template_window, repeat_run_limit):
                continue
            if not _gc_clamp_ok(template_window, gc_clamp_min, gc_clamp_max):
                continue
            if _has_hairpin_like(template_window):
                continue

            tm = _calc_tm(template_window)
            if tm is None:
                continue
            if not (tm_target - tm_tolerance <= tm <= tm_target + tm_tolerance):
                continue

            gc = gc_percent(template_window)
            if not (gc_min <= gc <= gc_max):
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

    return result


def validate_primer_pair(
    sequence: str,
    forward_seq: str,
    reverse_seq: str,
    *,
    product_min: int,
    product_max: int,
    tm_target: float,
    tm_tolerance: float,
) -> dict[str, Any]:
    sequence = sequence.upper()
    forward = _normalize_primer(forward_seq)
    reverse = _normalize_primer(reverse_seq)

    if not forward or not reverse:
        return {
            "valid": False,
            "errors": ["both forward and reverse primer sequences are required"],
            "pairs": [],
        }

    if not _VALID_DNA.fullmatch(forward) or not _VALID_DNA.fullmatch(reverse):
        return {
            "valid": False,
            "errors": ["primers must contain only A/C/G/T/N characters"],
            "pairs": [],
        }

    forward_tm = _calc_tm(forward)
    reverse_tm = _calc_tm(reverse)
    if forward_tm is None or reverse_tm is None:
        return {
            "valid": False,
            "errors": ["failed to calculate Tm"],
            "pairs": [],
        }

    f_hits = _find_bindings(sequence, forward)
    r_hits = _find_bindings(sequence, reverse)
    if not f_hits:
        return {
            "valid": False,
            "errors": ["forward primer does not match template"],
            "pairs": [],
        }
    if not r_hits:
        return {
            "valid": False,
            "errors": ["reverse primer does not match template"],
            "pairs": [],
        }

    pairs: list[PrimerPairBinding] = []
    for f in f_hits:
        for r in r_hits:
            if f[2] == r[2]:
                continue
            forward_3 = _three_prime_pos(f[0], f[1], f[2])
            reverse_3 = _three_prime_pos(r[0], r[1], r[2])
            product_size = abs(reverse_3 - forward_3) + 1
            if product_size < product_min or product_size > product_max:
                continue

            seed_warnings = max(_seed_offtarget_count(sequence, forward[-5:]) - 1, 0) + max(
                _seed_offtarget_count(sequence, reverse[-5:]) - 1,
                0,
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
                    dimer_risk=_pair_dimer_risk(forward, reverse),
                )
            )

    pairs.sort(key=lambda item: (item.tm_gap, item.product_size, item.seed_warnings))
    if not pairs:
        return {
            "valid": False,
            "errors": ["no valid product found under configured range"],
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


def _has_hairpin_like(sequence: str) -> bool:
    seq = sequence.upper()
    for k in range(4, min(len(seq), 8) + 1):
        tail = seq[-k:]
        rc = str(Seq(tail).reverse_complement())
        if rc in seq[:-1]:
            return True
    return False


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


def _pair_dimer_risk(seq_a: str, seq_b: str) -> bool:
    a = seq_a.upper()
    b = seq_b.upper()
    max_len = min(7, len(a), len(b))
    for k in range(4, max_len + 1):
        if str(Seq(a[-k:]).reverse_complement()) in b:
            return True
        if str(Seq(b[-k:]).reverse_complement()) in a:
            return True
    return False
