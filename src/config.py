from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


DEFAULT_FEATURES = [
    "annotation",
    "extreme_gc",
    "homopolymer",
    "ambiguous",
]

DEFAULT_TIMEOUT = 20.0
DEFAULT_RETRIES = 5
DEFAULT_CACHE_TTL_HOURS = 24
DEFAULT_FLANK = 10_000

NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# Reference genome used only by optional legacy CLI resolver
REFERENCE_ACCESSION = "CP076470"
REFERENCE_TAXID = 511145
REFERENCE_NAME = "Escherichia coli reference"

ENSEMBL_SEQUENCE_REGION = "https://rest.ensembl.org/sequence/region/{species}/{region}"
ENSEMBL_SEQUENCE_ID = "https://rest.ensembl.org/sequence/id/{ensembl_id}"

OUTPUT_FILE_SUFFIX = ".negfeatures.gb"
PRIMER_OUTPUT_FILE_SUFFIX = ".primer_design.gb"
OUTPUT_DIR = Path("data/output")
CACHE_DIR = Path("data/cache")

ENSEMBL_SEQUENCE_MAX_BP = 10_000_000
ENSEMBL_SEQUENCE_SAFETY_BP = 9_500_000

USER_AGENT = "PrimerMaker/1.0.0 (+https://github.com/)"


# Primer maker defaults
DEFAULT_PRIMER_TM_TARGET = 55.0
DEFAULT_PRIMER_TM_TOLERANCE = 2.5
DEFAULT_IDEAL_PRIMER_TM_MIN = 50.0
DEFAULT_IDEAL_PRIMER_TM_MAX = 65.0
DEFAULT_IDEAL_REPEAT_UNIT_MAX = 4
DEFAULT_PRIMER_TM_GAP_MAX = 5.0
DEFAULT_PRIMER_LEN_MIN = 18
DEFAULT_PRIMER_LEN_MAX = 24
DEFAULT_GC_MIN = 40.0
DEFAULT_GC_MAX = 60.0
DEFAULT_GC_CLAMP_MIN = 1
DEFAULT_GC_CLAMP_MAX = 2
DEFAULT_MANUAL_TM_GAP_FAIL = 5.0
DEFAULT_MANUAL_HAIRPIN_MIN_K = 4
DEFAULT_MANUAL_HAIRPIN_MAX_K = 4
DEFAULT_MANUAL_SELF_DIMER_MIN_OVERLAP = 4
DEFAULT_MANUAL_SELF_DIMER_MAX_OVERLAP = 5
DEFAULT_MANUAL_SELF_DIMER_EXCLUDE_IDENTICAL_WINDOW = True
DEFAULT_MANUAL_PAIR_DIMER_MIN_OVERLAP = 4
DEFAULT_MANUAL_PAIR_DIMER_MAX_OVERLAP = 5
DEFAULT_MANUAL_OFFTARGET_SEED_LEN = 5
DEFAULT_MANUAL_OFFTARGET_SEED_WARNING_LIMIT = 2
DEFAULT_MANUAL_REQUIRE_3P_DIMER = True

DEFAULT_INTERFERENCE_GC_WINDOW = 50
DEFAULT_INTERFERENCE_GC_STEP = 10
DEFAULT_INTERFERENCE_GC_MIN = 30.0
DEFAULT_INTERFERENCE_GC_MAX = 70.0
DEFAULT_INTERFERENCE_HOMOPOLYMER_AT = 4
DEFAULT_INTERFERENCE_HOMOPOLYMER_GC = 4
DEFAULT_INTERFERENCE_HOMOPOLYMER_STEP = 1
DEFAULT_INTERFERENCE_REPEAT_RUN = 4
DEFAULT_MIN_PRODUCT_SIZE = 70
DEFAULT_MAX_PRODUCT_SIZE = 3000


@dataclass(frozen=True)
class FeatureScanOptions:
    maf_threshold: float = 0.01
    gc_window: int = 50
    gc_step: int = 10
    gc_min: float = 30.0
    gc_max: float = 70.0
    homopolymer_at: int = 5
    homopolymer_gc: int = 4


GENBANK_FEATURE_MAP = {
    "annotation": "misc_feature",
    "gene": "gene",
    "cds": "CDS",
    "tRNA": "tRNA",
    "rRNA": "rRNA",
    "nc_rna": "misc_RNA",
    "rep_origin": "rep_origin",
    "repeat_region": "repeat_region",
    "repeat": "repeat_region",
    "simple": "repeat_region",
    "variation": "variation",
    "structural_variation": "misc_feature",
    "extreme_gc": "misc_feature",
    "homopolymer": "misc_feature",
    "ambiguous": "misc_feature",
}
