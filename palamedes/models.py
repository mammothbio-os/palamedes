from typing import NamedTuple


class Block(NamedTuple):
    """
    General purpose block for slices of the alignment. Note (start, end) should be ZBHO
    """

    start: int
    end: int
    bases: str

    @classmethod
    def collapse(cls, blocks: list["Block"]) -> "Block":
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
    reference_blocks: list[Block]
    alternate_blocks: list[Block]
