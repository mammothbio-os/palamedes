import logging
from functools import reduce

from Bio.Align import Alignment, PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from palamedes.config import (
    ALIGNMENT_GAP_CHAR,
    DEFAULT_EXTEND_GAP_SCORE,
    DEFAULT_MATCH_SCORE,
    DEFAULT_MISMATCH_SCORE,
    DEFAULT_OPEN_GAP_SCORE,
    GLOBAL_ALIGN_MODE,
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


def merge_reduce(blocks: list[VariantBlock], next_block: VariantBlock) -> list[VariantBlock]:
    """Helper function to implement merge or append, passed to functools.reduce to merge blocks"""
    peek_block = blocks[-1]

    if can_merge_variant_blocks(peek_block, next_block):
        blocks[-1] = merge_variant_blocks(peek_block, next_block)
    else:
        blocks.append(next_block)

    return blocks


def generate_alignment(
    reference_seq_record: SeqRecord,
    alternate_seq_record: SeqRecord,
    molecule_type: str = MOLECULE_TYPE_PROTEIN,
    aligner: PairwiseAligner | None = None,
) -> Alignment:
    """
    Using biopython's PairwiseAligner, generate an alignment object representing the best alignment
    between the 2 biopython SeqRecords. Note that the molecule_type argument will be used to specify a molecule_type
    annotation on the SeqRecord.

    By default the function creates an aligner object using the defaults, but the caller may provide their
    own pre-configured aligner. This aligner must be set to 'global' mode.

    Note that it is possible for multiple alignments to be returned with the same max score. Synthetic testing
    has shown that generally the highest scoring alignments are returned in "left to right" order. Take the
    following example, which was for testing HGVS duplications. This was the "intended" alignment:

        ATC---GGGGGGGG
        ATCATCGGGGGGGG

    However, there are multiple ways to get a three base insertion, all with the same score:
        ipdb> for aln in aligner.align('ATCGGGGGGGG', 'ATCATCGGGGGGGG'):
        print(aln)

        target            0 ---ATCGGGGGGGG 11
                          0 ---||||||||||| 14
        query             0 ATCATCGGGGGGGG 14

        target            0 A---TCGGGGGGGG 11
                          0 |---|||||||||| 14
        query             0 ATCATCGGGGGGGG 14

        target            0 AT---CGGGGGGGG 11
                          0 ||---||||||||| 14
        query             0 ATCATCGGGGGGGG 14

        target            0 ATC---GGGGGGGG 11
                          0 |||---|||||||| 14
        query             0 ATCATCGGGGGGGG 14

    Due to the HGVS rule of always indicating that the 3' most base is considered the modified one, a best effort
    attempt is made to return the most "right aligned" alignment, by returning the last alignment with the highest
    score. This way not yield the ideal results in more complicated cases.
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

    if (ref_molecule_type := reference_seq_record.annotations.get("molecule_type")) != molecule_type:
        raise ValueError(
            "Cannot generate alignment, reference_seq_record is a SeqRecord an invalid molecule_type annotation "
            f"got: {ref_molecule_type}, expected: {molecule_type}!"
        )

    if (alt_molecule_type := alternate_seq_record.annotations.get("molecule_type")) != molecule_type:
        raise ValueError(
            "Cannot generate alignment, alternate_seq_record is a SeqRecord an invalid molecule_type annotation "
            f"got: {alt_molecule_type}, expected: {molecule_type}!"
        )

    alignments = aligner.align(reference_seq_record, alternate_seq_record)
    best_alignment = next(alignments)
    best_alignment_score = best_alignment.score

    other_alignments_with_best = [alignment for alignment in alignments if alignment.score == best_alignment_score]
    if len(other_alignments_with_best) > 0:
        LOGGER.debug(
            "Found %s alignments with max score, returning last in the list (3' end rule)",
            len(other_alignments_with_best) + 1,
        )
        best_alignment = other_alignments_with_best[-1]

    return best_alignment


def generate_variant_blocks(alignment: Alignment) -> list[VariantBlock]:
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
    merged_blocks = reduce(merge_reduce, single_position_variant_blocks[1:], [single_position_variant_blocks[0]])

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
