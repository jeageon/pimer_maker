from __future__ import annotations

from dataclasses import dataclass
from io import StringIO
import re
from typing import Any, Optional

from Bio import SeqIO

from ..config import DH5A_ACCESSION, DH5A_NAME, DH5A_TAXID, NCBI_EFETCH
from ..models.data_schemas import GenomicCoordinates
from ..utils.coord_utils import apply_flank
from ..utils.exceptions import NoMappingError, ToolError
from ..utils.api_client import ApiClient


@dataclass
class ResolverResult:
    coordinates: GenomicCoordinates
    warnings: list[str]


_GENE_FEATURE_TYPES = {
    "gene",
    "cds",
    "mrna",
    "exon",
    "rrna",
    "trna",
    "ncrna",
    "misc_feature",
    "mobile_element",
    "rep_origin",
    "repeat_region",
}


def _to_int(value: Any, default: Optional[int] = None) -> Optional[int]:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def _normalize_token(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, list):
        value = value[0] if value else ""
    token = str(value).strip().lower()
    token = token.replace("_", "")
    token = token.replace("-", "")
    token = re.sub(r"\s+", "", token)
    return token


def _tokenize_qualifier(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, list):
        tokens: list[str] = []
        for item in value:
            tokens.extend(_tokenize_qualifier(item))
        return [token for token in tokens if token]
    text = str(value).strip()
    if not text:
        return []
    raw = re.split(r"[,;/]|\band\b", text, flags=re.IGNORECASE)
    return [part.strip() for part in raw if part.strip()]


def _compact_qualifier(values: Any) -> str:
    if values is None:
        return ""
    if isinstance(values, list):
        merged = [str(v).strip() for v in values if str(v).strip()]
        return "; ".join(dict.fromkeys(merged))
    return str(values).strip()


class CoordinateResolver:
    def __init__(self, api_client: ApiClient) -> None:
        self.api = api_client
        self._record = None

    def resolve(
        self,
        uniprot_id: str,
        flank_bp: int,
        flank_mode: str = "genomic",
        assembly_preference: str = "auto",
        taxid_filter: Optional[int] = None,
    ) -> ResolverResult:
        del assembly_preference
        del taxid_filter
        query = str(uniprot_id or "").strip()
        if not query:
            raise NoMappingError("gene name is empty")

        record = self._load_reference_record()
        warnings: list[str] = []
        candidates = self._find_gene_candidates(record, query)
        if not candidates:
            suggestions = ", ".join(sorted(set(self._collect_common_gene_names(record)))[:20])
            warn = "No direct match found in CP076470"
            if suggestions:
                warn = f"{warn}. 가까운 후보: {suggestions}"
            raise NoMappingError(warn)

        best_score = candidates[0]["score"]
        best = candidates[0]["feature"]
        match_name = candidates[0]["matched"]
        if best_score < 50:
            warnings.append(
                f"No exact match for '{query}', selected '{match_name}' by partial match. "
                "기록자 이름/alias가 다를 경우 매칭 결과가 다를 수 있습니다."
            )

        gene_start = _to_int(best.location.start) + 1
        gene_end = _to_int(best.location.end)
        if gene_start is None or gene_end is None:
            raise NoMappingError(f"Invalid coordinate span for matched feature in CP076470: {match_name}")
        if gene_start > gene_end:
            gene_start, gene_end = gene_end, gene_start

        strand = _to_int(best.location.strand, 1)
        if strand not in (-1, 1):
            strand = 1
        ext_start, ext_end = apply_flank(gene_start, gene_end, flank_bp, flank_mode, strand)

        genome_len = len(record.seq)
        if ext_start < 1:
            warnings.append(f"flank start {ext_start} is below 1; clamped to 1")
            ext_start = 1
        if ext_end > genome_len:
            warnings.append(f"flank end {ext_end} exceeds genome length; clamped to {genome_len}")
            ext_end = genome_len
        if ext_start > ext_end:
            raise NoMappingError("Failed to build extraction region after clamping to genome boundaries")

        display_name = self._choose_display_name(best, match_name)
        annotations = self._collect_region_annotations(record=record, region_start=ext_start, region_end=ext_end)

        return ResolverResult(
            coordinates=GenomicCoordinates(
                uniprot_id=query,
                ensembl_gene_id=match_name,
                coordinate_source="ncbi",
                query_gene=query,
                query_type="gene_name",
                ncbi_accession=DH5A_ACCESSION,
                species=DH5A_NAME,
                assembly_name=DH5A_NAME,
                seq_region_name=DH5A_ACCESSION,
                gene_start_1based=gene_start,
                gene_end_1based=gene_end,
                strand=strand,
                display_name=display_name,
                taxid=DH5A_TAXID,
                ncbi_genome_length=genome_len,
                ncbi_annotations=annotations,
                ext_start_1based=ext_start,
                ext_end_1based=ext_end,
            ),
            warnings=warnings,
        )

    def _load_reference_record(self):
        if self._record is not None:
            return self._record

        response = self.api.get(
            NCBI_EFETCH,
            headers={"Accept": "text/plain"},
            params={
                "db": "nuccore",
                "id": DH5A_ACCESSION,
                "rettype": "gb",
                "retmode": "text",
            },
        )
        if not response.text:
            raise ToolError(f"No sequence returned for reference accession {DH5A_ACCESSION}")

        try:
            record = SeqIO.read(StringIO(response.text), "genbank")
        except Exception as exc:
            raise ToolError(f"Failed to parse GenBank record for {DH5A_ACCESSION}: {exc}")

        self._record = record
        return record

    def _find_gene_candidates(self, record, query: str) -> list[dict[str, Any]]:
        normalized_query = _normalize_token(query)
        bests: list[dict[str, Any]] = []

        for feature in record.features:
            if feature.type.lower() not in _GENE_FEATURE_TYPES:
                continue
            if not hasattr(feature, "location") or feature.location is None:
                continue

            candidates = self._collect_feature_tokens(feature)
            if not candidates:
                continue

            score, matched = self._score_feature(normalized_query, candidates)
            if score <= 0:
                continue
            bests.append({
                "feature": feature,
                "score": score,
                "matched": matched,
                "tokens": candidates,
            })

        bests.sort(
            key=lambda item: (
                item["score"],
                item["feature"].type.lower() == "gene",
                item["feature"].location.strand or 0,
            ),
            reverse=True,
        )
        return bests

    def _collect_feature_tokens(self, feature) -> list[str]:
        tokens: list[str] = []
        qualifiers = getattr(feature, "qualifiers", {}) or {}
        for key in (
            "gene",
            "locus_tag",
            "old_locus_tag",
            "protein_id",
            "product",
            "gene_synonym",
            "name",
            "note",
            "db_xref",
        ):
            values = qualifiers.get(key)
            for item in _tokenize_qualifier(values):
                token = item.strip()
                if token:
                    tokens.append(token)
        return list(dict.fromkeys(tokens))

    def _score_feature(self, query_norm: str, tokens: list[str]) -> tuple[int, str]:
        best_score = 0
        matched = ""

        for token in tokens:
            token_norm = _normalize_token(token)
            if not token_norm:
                continue

            if token_norm == query_norm:
                score = 200
                if score > best_score:
                    best_score = score
                    matched = token
                continue

            if query_norm in token_norm:
                score = 110
            elif token_norm in query_norm:
                score = 90
            elif token_norm.startswith(query_norm) or query_norm.startswith(token_norm):
                score = 80
            else:
                score = 0

            if score > best_score:
                best_score = score
                matched = token

        return best_score, matched


    def _choose_display_name(self, feature, fallback: str) -> str:
        qualifiers = getattr(feature, "qualifiers", {}) or {}
        for key in ("gene", "locus_tag", "old_locus_tag", "product"):
            values = qualifiers.get(key)
            if not values:
                continue
            for item in _tokenize_qualifier(values):
                if item:
                    return str(item).strip()
        return fallback

    def _collect_region_annotations(self, record, region_start: int, region_end: int) -> list[dict[str, Any]]:
        annotations: list[dict[str, Any]] = []
        for feature in record.features:
            if feature.type.lower() == "source":
                continue
            if feature.type.lower() not in _GENE_FEATURE_TYPES:
                continue
            if not hasattr(feature, "location") or feature.location is None:
                continue

            start = _to_int(feature.location.start, 1) + 1
            end = _to_int(feature.location.end)
            if start is None or end is None:
                continue
            if start > region_end or end < region_start:
                continue
            if start < 1:
                start = 1

            qualifiers = getattr(feature, "qualifiers", {}) or {}
            annotations.append(
                {
                    "feature_type": feature.type.lower(),
                    "start": max(start, region_start),
                    "end": min(end, region_end),
                    "strand": _to_int(feature.location.strand, 0),
                    "display_name": self._choose_display_name(feature, feature.type),
                    "qualifiers": {str(k): _compact_qualifier(v) for k, v in qualifiers.items() if _compact_qualifier(v)},
                }
            )
        return annotations

    def _collect_common_gene_names(self, record) -> list[str]:
        names: list[str] = []
        for feature in record.features:
            if feature.type.lower() not in {"gene", "cds"}:
                continue
            qualifiers = getattr(feature, "qualifiers", {}) or {}
            for key in ("gene", "locus_tag", "old_locus_tag", "product", "gene_synonym"):
                for item in _tokenize_qualifier(qualifiers.get(key)):
                    if item:
                        names.append(item.strip())
        return names

