from Bio.Align import Alignment

from palamedes.models import VariantBlock
from palamedes.config import (
    HGVS_VARIANT_TYPE_SUBSTITUTION,
    HGVS_VARIANT_TYPE_DELETION,
    HGVS_VARIANT_TYPE_EXTENSION,
    HGVS_VARIANT_TYPE_DUPLICATION,
    HGVS_VARIANT_TYPE_REPEAT,
    HGVS_VARIANT_TYPE_INSERTION,
    HGVS_VARIANT_TYPE_DELETION_INSERTION,
    VARIANT_BASE_MISMATCH,
    VARIANT_BASE_DELETION,
    VARIANT_BASE_INSERTION,
)
from palamedes.utils import contains_repeated_substring, yield_repeating_substrings


def categorize_variant_block(alignment: Alignment, variant_block: VariantBlock) -> str:
    """
    First pass over variant blocks to categorize them with their base HGVS "type".

    The following rule-set is used:
    - A single position mismatch => Substitution
    - 1 or more deletion positions and nothing else => Deletion
    - 1 or more insertion positions and nothing else
        - If fully upstream of downstream of reference => Extension
        - If inserted bases perfectly match upstream bases => Duplication
        - If inserted bases are comprised of a repeating sub-string that perfectly matches upstream bases => Repeat
        - If nothing else => Insertion
    - Any other combo of calls is a Deletion-Insertion
    """
    if variant_block.alignment_block.bases == VARIANT_BASE_MISMATCH:
        return HGVS_VARIANT_TYPE_SUBSTITUTION

    if set(variant_block.alignment_block.bases) == set([VARIANT_BASE_DELETION]):
        return HGVS_VARIANT_TYPE_DELETION

    if set(variant_block.alignment_block.bases) == set([VARIANT_BASE_INSERTION]):
        if variant_block.alignment_block.start == 0 or variant_block.alignment_block.end == len(alignment[0]):
            return HGVS_VARIANT_TYPE_EXTENSION

        inserted_bases = variant_block.alternate_blocks[0].bases

        def get_upstream_reference_sequence(num_bases):
            start_pos = variant_block.alignment_block.start - num_bases
            end_pos = variant_block.alignment_block.start
            return alignment[0][start_pos:end_pos]

        if inserted_bases == get_upstream_reference_sequence(len(inserted_bases)):
            return HGVS_VARIANT_TYPE_DUPLICATION

        if contains_repeated_substring(inserted_bases):
            candidate_repeats = [
                substring
                for substring in yield_repeating_substrings(inserted_bases)
                if get_upstream_reference_sequence(len(substring)) == substring
            ]
            if len(candidate_repeats) > 0:
                return HGVS_VARIANT_TYPE_REPEAT

        return HGVS_VARIANT_TYPE_INSERTION

    return HGVS_VARIANT_TYPE_DELETION_INSERTION
