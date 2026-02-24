from __future__ import annotations

import re
from collections.abc import Iterator


_AMBIGUOUS_PATTERN = re.compile(r"[^ATGCatgc]+")


def count_invalid_bases(seq: str) -> int:
    return len(_AMBIGUOUS_PATTERN.findall(seq))


def scan_extreme_gc_windows(
    sequence: str,
    window_size: int,
    step: int,
    gc_min: float,
    gc_max: float,
) -> list[tuple[int, int, float]]:
    seq = sequence.upper()
    windows: list[tuple[int, int, float]] = []
    n = len(seq)
    if window_size <= 0 or step <= 0 or n < window_size:
        return windows

    for start in range(0, n - window_size + 1, step):
        win = seq[start : start + window_size]
        gc = (win.count("G") + win.count("C")) / max(len(win), 1) * 100.0
        if gc < gc_min or gc > gc_max:
            windows.append((start, start + window_size, gc))
    return windows


def merge_intervals_with_gap(intervals: list[tuple[int, int, float]], gap: int = 0) -> list[tuple[int, int, float]]:
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda item: item[0])
    merged: list[tuple[int, int, float]] = []
    cur_start, cur_end, cur_score = intervals[0]
    for start, end, score in intervals[1:]:
        if start <= cur_end + gap:
            cur_end = max(cur_end, end)
            cur_score = (cur_score + score) / 2.0
        else:
            merged.append((cur_start, cur_end, cur_score))
            cur_start, cur_end, cur_score = start, end, score
    merged.append((cur_start, cur_end, cur_score))
    return merged


def scan_homopolymers(sequence: str, at_run: int = 5, gc_run: int = 4) -> list[tuple[str, int, int]]:
    seq = sequence.upper()
    at_pattern = re.compile(rf"A{{{at_run},}}|T{{{at_run},}}", re.IGNORECASE)
    gc_pattern = re.compile(rf"G{{{gc_run},}}|C{{{gc_run},}}", re.IGNORECASE)

    hits: list[tuple[str, int, int]] = []
    for match in at_pattern.finditer(seq):
        start, end = match.span()
        hits.append((match.group(0)[0].upper(), start, end))
    for match in gc_pattern.finditer(seq):
        start, end = match.span()
        hits.append((match.group(0)[0].upper(), start, end))
    return sorted(hits, key=lambda x: x[1])


def scan_ambiguous(sequence: str) -> list[tuple[int, int]]:
    raw = [(m.start(), m.end()) for m in _AMBIGUOUS_PATTERN.finditer(sequence)]
    if not raw:
        return []
    merged: list[tuple[int, int]] = []
    start, end = raw[0]
    for s, e in raw[1:]:
        if s <= end:
            end = max(end, e)
        else:
            merged.append((start, end))
            start, end = s, e
    merged.append((start, end))
    return merged


def gc_percent(seq: str) -> float:
    if not seq:
        return 0.0
    gc = sum(1 for c in seq if c in "GgCc")
    return (gc / len(seq)) * 100.0


def sliding_windows(seq: str, window: int, step: int) -> Iterator[tuple[int, str]]:
    if window <= 0 or step <= 0:
        raise ValueError("window and step must be positive")
    for i in range(0, max(0, len(seq) - window + 1), step):
        yield i, seq[i : i + window]


def find_homopolymers(seq: str, at_threshold: int = 5, gc_threshold: int = 4) -> list[tuple[int, int, str, int]]:
    hits: list[tuple[int, int, str, int]] = []
    for base, start, end in scan_homopolymers(seq, at_run=at_threshold, gc_run=gc_threshold):
        hits.append((start, end, base, end - start))
    return sorted(hits, key=lambda item: (item[0], item[1]))


def find_ambiguous_runs(seq: str) -> list[tuple[int, int]]:
    return scan_ambiguous(seq)
