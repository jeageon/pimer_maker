from __future__ import annotations

from dataclasses import dataclass


def ensembl_to_relative(
    feature_start_1based: int,
    feature_end_1based: int,
    ext_start_1based: int,
    seq_len: int,
) -> tuple[int, int]:
    rel_start = feature_start_1based - ext_start_1based
    rel_end = feature_end_1based - ext_start_1based + 1
    if rel_start < 0:
        rel_start = 0
    if rel_end > seq_len:
        rel_end = seq_len
    return rel_start, rel_end


def ensembl_to_rel0(
    feature_start_1based: int,
    feature_end_1based: int,
    ext_start_1based: int,
    seq_len: int,
) -> tuple[int, int]:
    return ensembl_to_relative(
        feature_start_1based=feature_start_1based,
        feature_end_1based=feature_end_1based,
        ext_start_1based=ext_start_1based,
        seq_len=seq_len,
    )


def apply_flank(
    gene_start_1based: int,
    gene_end_1based: int,
    flank_bp: int,
    flank_mode: str,
    strand: int,
    sequence_start_min: int = 1,
) -> tuple[int, int]:
    if flank_mode == "strand_relative":
        if strand not in (-1, 1):
            flank_mode = "genomic"
    if flank_mode in {"genomic", "strand_relative"}:
        start = max(sequence_start_min, gene_start_1based - flank_bp)
        end = gene_end_1based + flank_bp
    else:
        start = max(sequence_start_min, gene_start_1based - flank_bp)
        end = gene_end_1based + flank_bp
    return start, end


def build_region_string(chr_name: str, start_1based: int, end_1based: int, strand: int) -> str:
    return f"{chr_name}:{start_1based}..{end_1based}:{strand}"


def clamp_region_length(
    start_1based: int,
    end_1based: int,
    max_len: int,
) -> tuple[int, int, bool]:
    length = end_1based - start_1based + 1
    if length <= max_len:
        return start_1based, end_1based, False
    shrink = length - max_len
    new_start = start_1based + (shrink // 2)
    new_end = end_1based - (shrink - (shrink // 2))
    if new_start < 1:
        new_start = 1
        new_end = max(new_start, max_len)
    return new_start, new_end, True


@dataclass
class Chunk:
    start: int
    end: int


def build_chunks(region_start: int, region_end: int, chunk_size: int) -> list[Chunk]:
    if region_start > region_end:
        return []
    chunks: list[Chunk] = []
    current = region_start
    while current <= region_end:
        chunk_end = min(region_end, current + chunk_size - 1)
        chunks.append(Chunk(current, chunk_end))
        current = chunk_end + 1
    return chunks
