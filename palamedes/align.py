import logging
from functools import reduce
from typing import List, NamedTuple, Optional, Union

from Bio.Align import Alignment, PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from palamedes.config import (
    ALIGNMENT_GAP_CHAR,
    ALT_SEQUENCE_ID,
    DEFAULT_EXTEND_GAP_SCORE,
    DEFAULT_MATCH_SCORE,
    DEFAULT_MISMATCH_SCORE,
    DEFAULT_OPEN_GAP_SCORE,
    GLOBAL_ALIGN_MODE,
    REF_SEQUENCE_ID,
    VARIANT_BASE_DELETION,
    VARIANT_BASE_INSERTION,
    VARIANT_BASE_MATCH,
    VARIANT_BASE_MISMATCH,
)

LOGGER = logging.getLogger(__name__)


class Block(NamedTuple):
    """
    General purpose block for slices of the alignment. Note (start, end) should be ZBHO
    """

    start: int
    end: int
    bases: str

    @classmethod
    def collapse(cls, blocks: List["Block"]) -> "Block":
        """
        Given a variable length list of blocks, collapse all blocks into one or raise a ValueError if the list is empty.
        Blocks will be sorted but must be adjacent to each-other.
        """
        if len(blocks) == 0:
            raise ValueError("Cannot collapse empty list of Blocks")

        sorted_blocks = sorted(blocks)
        for idx, block in enumerate(sorted_blocks[:-1]):
            if sorted_blocks[idx + 1].start != block.end:
                raise ValueError(f"Cannot collapse Blocks, they must be adjacent. Blocks {idx} and {idx + 1} are not!")

        new_start = sorted_blocks[0].start
        new_end = sorted_blocks[-1].end
        new_bases = "".join([block.bases for block in sorted_blocks])
        return Block(new_start, new_end, new_bases)


class VariantBlock(NamedTuple):
    """
    Internal representation of a slice of an alignment, based on a contiguous
    run of positions that are not matches. The object tracks the start and
    end positions from each sequence which are covered by the block, as well
    as the bases. It also tracks the start and end of the block within the
    alignment and the "variant bases", which is similar to a CIGAR string
    but without the short-hand notation.

    Note that:
    - reference_blocks and alternate_blocks store Block objects corresponding to the block from either
      sequence in the alignment. The list should have 0 or 1 elements. The 0 is if no part of the sequence
      is included in the alignment (like an insertion upstream of the first base).
    - Alignment block will always be set, and should have a sequence with the VARIANT_BASE characters
    """

    alignment_block: Block
    reference_blocks: List[Block]
    alternate_blocks: List[Block]


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


def merge_reduce(blocks: List[VariantBlock], next_block: VariantBlock) -> List[VariantBlock]:
    """Helper function to implement merge or append, passed to functools.reduce to merge blocks"""
    peek_block = blocks[-1]

    if can_merge_variant_blocks(peek_block, next_block):
        blocks[-1] = merge_variant_blocks(peek_block, next_block)
    else:
        blocks.append(next_block)

    return blocks


def generate_alignment(
    reference_sequence: Union[str, SeqRecord],
    alternate_sequence: Union[str, SeqRecord],
    aligner: Optional[PairwiseAligner] = None,
) -> Alignment:
    """
    Using biopython's PairwiseAligner, generate an alignment object representing the best alignment
    between the 2 inputs (which can either be raw string sequences or biopython SeqRecords).

    By default the function creates an aligner object using the defaults, but the caller may provide their
    own pre-configured aligner. This aligner must be set to 'global' mode.

    A warning is generated if a tie is detected between the 1st and 2nd best alignments.
    """
    if aligner is not None:
        if aligner.mode != GLOBAL_ALIGN_MODE:
            raise ValueError(f"Custom PairwiseAligner must be set to global mode, got: {aligner.mode}")
    else:
        aligner = PairwiseAligner(
            mode=GLOBAL_ALIGN_MODE,
            match_score=DEFAULT_MATCH_SCORE,
            mismatch_score=DEFAULT_MISMATCH_SCORE,
            open_gap_score=DEFAULT_OPEN_GAP_SCORE,
            extend_gap_score=DEFAULT_EXTEND_GAP_SCORE,
        )

    ref_seq_obj = (
        reference_sequence
        if isinstance(reference_sequence, SeqRecord)
        else SeqRecord(Seq(reference_sequence), id=REF_SEQUENCE_ID)
    )
    alt_seq_obj = (
        alternate_sequence
        if isinstance(alternate_sequence, SeqRecord)
        else SeqRecord(Seq(alternate_sequence), id=ALT_SEQUENCE_ID)
    )

    alignments = aligner.align(ref_seq_obj, alt_seq_obj)
    best_alignment = next(alignments)

    try:
        next_alignnment = next(alignments)
        if best_alignment.score == next_alignnment.score:
            LOGGER.warning("Tie detected between best and second best alignment!")
    except StopIteration:
        pass

    return best_alignment


def generate_variant_blocks(alignment: Alignment) -> List[VariantBlock]:
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
    reference_indices, alternate_indices = alignment.indices
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
    merged_blocks = reduce(merge_reduce, single_position_variant_blocks[1:], [single_position_variant_blocks[0]])

    # filter out the matches
    return [
        merged_block for merged_block in merged_blocks if VARIANT_BASE_MATCH not in merged_block.alignment_block.bases
    ]
