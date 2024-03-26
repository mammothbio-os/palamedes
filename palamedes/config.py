GLOBAL_ALIGN_MODE: str = "global"

# default alignment params, meant to favor a single Substitution vs some Indel
DEFAULT_MATCH_SCORE: int = 1
DEFAULT_MISMATCH_SCORE: int = -1
DEFAULT_OPEN_GAP_SCORE: int = -1
DEFAULT_EXTEND_GAP_SCORE: float = -0.1

REF_SEQUENCE_ID: str = "ref"
ALT_SEQUENCE_ID: str = "alt"

ALIGNMENT_GAP_CHAR: str = "-"

VARIANT_BASE_MATCH: str = "M"
VARIANT_BASE_MISMATCH: str = "m"
VARIANT_BASE_DELETION: str = "d"
VARIANT_BASE_INSERTION: str = "i"

HGVS_VARIANT_TYPE_SUBSTITUTION: str = "substitution"
HGVS_VARIANT_TYPE_DELETION: str = "deletion"
HGVS_VARIANT_TYPE_EXTENSION: str = "extension"
HGVS_VARIANT_TYPE_DUPLICATION: str = "duplication"
HGVS_VARIANT_TYPE_REPEAT: str = "repeat"
HGVS_VARIANT_TYPE_INSERTION: str = "insertion"
HGVS_VARIANT_TYPE_DELETION_INSERTION: str = "deletion_insertion"

MOLECULE_TYPE_ANNOTATION_KEY: str = "molecule_type"
MOLECULE_TYPE_PROTEIN: str = "protein"
HGVS_TYPE_PROTEIN: str = "p"
