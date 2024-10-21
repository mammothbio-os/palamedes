import logging
from functools import reduce, partial

from Bio.Align import Alignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from palamedes.config import (
    ALIGNMENT_GAP_CHAR,
    VARIANT_BASE_DELETION,
    VARIANT_BASE_INSERTION,
    VARIANT_BASE_MATCH,
    VARIANT_BASE_MISMATCH,
    MOLECULE_TYPE_PROTEIN,
    MOLECULE_TYPE_ANNOTATION_KEY,
)
from palamedes.models import Block, VariantBlock

LOGGER = logging.getLogger(__name__)


def generate_seq_record(sequence: str, seq_id: str, molecule_type: str = MOLECULE_TYPE_PROTEIN) -> SeqRecord:
    """
    Helper function to generate a SeqRecord object from a raw input sequence. This also handles
    configuring the expected molecule_type annotation which is required for downstream steps.
    """
    return SeqRecord(
        Seq(sequence),
        id=seq_id,
        annotations={MOLECULE_TYPE_ANNOTATION_KEY: molecule_type},
    )


def reverse_seq_record(seq_record: SeqRecord) -> SeqRecord:
    """
    Helper function to copy a SeqRecord into a new one, with the sequence reversed. This is a best effort copy,
    which is only used internally to the alignment logic to more easily access the 3' end most representation of
    the best alignment.
    """
    return SeqRecord(
        Seq(seq_record.seq[::-1]),
        id=seq_record.id,
        name=seq_record.name,
        description=seq_record.description,
        dbxrefs=seq_record.dbxrefs[:],
        features=seq_record.features[:],
        annotations=seq_record.annotations.copy(),
        letter_annotations=seq_record.letter_annotations.copy(),
    )


def make_variant_base(ref_base: str, alt_base: str) -> str:
    """Helper function to generate the correct variant base given the ref and alt alignment bases"""
    if ref_base == alt_base:
        return VARIANT_BASE_MATCH

    if ref_base == ALIGNMENT_GAP_CHAR:
        return VARIANT_BASE_INSERTION

    if alt_base == ALIGNMENT_GAP_CHAR:
        return VARIANT_BASE_DELETION

    return VARIANT_BASE_MISMATCH


def can_merge_variant_blocks(left: VariantBlock, right: VariantBlock) -> bool:
    """
    Determine if 2 variant blocks can be merged, based on the following rules:
    - Must be directly adjacent (left.end == right.start)
    - No matches can appear in either
    """
    return all(
        [
            left.alignment_block.end == right.alignment_block.start,
            VARIANT_BASE_MATCH not in left.alignment_block.bases,
            VARIANT_BASE_MATCH not in right.alignment_block.bases,
        ]
    )


def merge_variant_blocks(left: VariantBlock, right: VariantBlock) -> VariantBlock:
    """
    Merge 2 variant blocks, returning the new one. This just requires that all
    sets of blocks (reference, alternate and alignment) be collapsed together
    and new reference and alt lists be created (since they can be empty).
    """
    reference_blocks = left.reference_blocks + right.reference_blocks
    if len(reference_blocks) > 0:
        new_reference_blocks = [Block.collapse(reference_blocks)]
    else:
        new_reference_blocks = []

    alternate_blocks = left.alternate_blocks + right.alternate_blocks
    if len(alternate_blocks) > 0:
        new_alternate_blocks = [Block.collapse(alternate_blocks)]
    else:
        new_alternate_blocks = []

    new_alignment_block = Block.collapse([left.alignment_block, right.alignment_block])

    return VariantBlock(new_alignment_block, new_reference_blocks, new_alternate_blocks)


def merge_reduce(
    blocks: list[VariantBlock], next_block: VariantBlock, split_consecutive_mismatches: bool = False
) -> list[VariantBlock]:
    """
    Helper function to implement merge or append, passed to functools.reduce to merge blocks. Optional
    split_consecutive_mismatches will keep consecutive mismatches from being merged, allowing for a
    chain of mismatches vs a delins. Note that this is non-standard behavior according to HGVS.
    """
    peek_block = blocks[-1]

    # base check for merge-ability
    base_can_merge = can_merge_variant_blocks(peek_block, next_block)

    # extra check for flag, merge if false OR merge if true and one of  left and right
    # are not single bp mismatches. If true and both are single base mismatches
    # then we do not merge, this 'splitting'
    should_split_mismatches = all(
        [
            split_consecutive_mismatches,
            peek_block.alignment_block.bases == VARIANT_BASE_MISMATCH,
            next_block.alignment_block.bases == VARIANT_BASE_MISMATCH,
        ]
    )

    if base_can_merge and not should_split_mismatches:
        blocks[-1] = merge_variant_blocks(peek_block, next_block)
    else:
        blocks.append(next_block)

    return blocks


def generate_variant_blocks(alignment: Alignment, split_consecutive_mismatches: bool = False) -> list[VariantBlock]:
    """
    Given a BioPython.Alignment object, parse the alignment to generate a list of VariantBlock objects.
    A VariantBlock is an internal object which represents a contiguous run of positions within the alignment
    which are not matches (mismatch, del or ins). These blocks will be categorized and converted into HGVS
    objects in a later step, but this intermediate representation is useful for debugging and testing.

    This function uses functools.reduce in order to structure the problem of merging blocks
    into a more understandable approach:
    - We define the logic for when 2 VariantBlocks can be merged together, see can_merge_variant_blocks
    - We define the logic for how to merge 2 VariantBlocks together, see merge_variant_blocks
    - We define a reduce function that treats the accumulator like a stack, peeking the head, checking
      for merge-ability and merging if possible at each iteration.
    """
    reference_indices = alignment.indices[0].tolist()
    alternate_indices = alignment.indices[1].tolist()
    bases_zipper = zip(alignment[0], alignment[1])

    # First generate a VariantBlock for all positions in the Alignment along with the relevant sequence Blocks
    single_position_variant_blocks = [
        VariantBlock(
            Block(idx, idx + 1, make_variant_base(reference_base, alternate_base)),
            []
            if reference_base == ALIGNMENT_GAP_CHAR
            else [
                Block(
                    reference_indices[idx],
                    reference_indices[idx] + 1,
                    reference_base,
                )
            ],
            []
            if alternate_base == ALIGNMENT_GAP_CHAR
            else [
                Block(
                    alternate_indices[idx],
                    alternate_indices[idx] + 1,
                    alternate_base,
                )
            ],
        )
        for idx, (reference_base, alternate_base) in enumerate(bases_zipper)
    ]

    # use reduce with a helper function to merge valid/adjacent blocks together
    wrapped_merge_func = partial(merge_reduce, split_consecutive_mismatches=split_consecutive_mismatches)
    merged_blocks = reduce(wrapped_merge_func, single_position_variant_blocks[1:], [single_position_variant_blocks[0]])

    # filter out the matches
    return [
        merged_block for merged_block in merged_blocks if VARIANT_BASE_MATCH not in merged_block.alignment_block.bases
    ]


def get_upstream_reference_sequence(alignment: Alignment, anchor_position: int, num_bases: int) -> str:
    """
    Get num_bases of the reference sequence directly upstream of a variant block. This is done by treating
    the starting position of the variant block as an anchor into the aligned/gapped reference sequence.
    A slice is then taken from 0 -> the anchor position, and gaps are removed from this slice. Finally,
    return the last num_bases of the substring. This is done to avoid cases where an upstream gap
    may cause a more naive approach to miss a duplication or repeat (since a gap might appear in the checked sequence).
    """
    ungapped_ref = alignment[0][:anchor_position].replace(ALIGNMENT_GAP_CHAR, "")
    return ungapped_ref[-num_bases:]
