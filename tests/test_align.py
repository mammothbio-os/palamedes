from palamedes.align import (
    make_variant_base,
    can_merge_variant_blocks,
    merge_variant_blocks,
    generate_seq_record,
    generate_variant_blocks,
)
from palamedes.config import (
    REF_SEQUENCE_ID,
    VARIANT_BASE_MATCH,
    ALIGNMENT_GAP_CHAR,
    VARIANT_BASE_INSERTION,
    VARIANT_BASE_DELETION,
    VARIANT_BASE_MISMATCH,
    MOLECULE_TYPE_PROTEIN,
    MOLECULE_TYPE_ANNOTATION_KEY,
)
from palamedes.models import VariantBlock, Block
from tests.base import PalamedesBaseCase


class GenerateSeqRecordsTestCase(PalamedesBaseCase):
    def test_generate_seq_record(self):
        seq = "AAA"
        seq_rec = generate_seq_record(seq, REF_SEQUENCE_ID)

        self.assertEqual(seq_rec.seq, seq)
        self.assertEqual(seq_rec.id, REF_SEQUENCE_ID)
        self.assertEqual(seq_rec.annotations, {MOLECULE_TYPE_ANNOTATION_KEY: MOLECULE_TYPE_PROTEIN})

    def test_generate_seq_record_custom_molecule_type(self):
        seq = "AAA"
        custom_molecule_type = "dna"

        seq_rec = generate_seq_record(seq, REF_SEQUENCE_ID, molecule_type=custom_molecule_type)
        self.assertEqual(seq_rec.annotations, {MOLECULE_TYPE_ANNOTATION_KEY: custom_molecule_type})


class MakeVariantBaseTestCase(PalamedesBaseCase):
    def test_make_variant_base_match(self):
        self.assertEqual(make_variant_base("A", "A"), VARIANT_BASE_MATCH)

    def test_make_variant_base_mismatch(self):
        self.assertEqual(make_variant_base("A", "T"), VARIANT_BASE_MISMATCH)

    def test_make_variant_base_deletion(self):
        self.assertEqual(make_variant_base("A", ALIGNMENT_GAP_CHAR), VARIANT_BASE_DELETION)

    def test_make_variant_base_insertion(self):
        self.assertEqual(make_variant_base(ALIGNMENT_GAP_CHAR, "T"), VARIANT_BASE_INSERTION)


class CanMergeVariantBlocksTestCase(PalamedesBaseCase):
    def make_variant_block(self, start: int, end: int, bases: str) -> VariantBlock:
        return VariantBlock(Block(start, end, bases), [], [])

    def test_can_merge_variant_blocks_not_adjacent(self):
        left = self.make_variant_block(0, 4, VARIANT_BASE_DELETION * 4)
        right = self.make_variant_block(6, 8, VARIANT_BASE_MISMATCH * 2)
        self.assertFalse(can_merge_variant_blocks(left, right))

    def test_can_merge_variant_blocks_left_has_match(self):
        left = self.make_variant_block(0, 4, VARIANT_BASE_DELETION * 3 + VARIANT_BASE_MATCH)
        right = self.make_variant_block(4, 6, VARIANT_BASE_MISMATCH * 2)
        self.assertFalse(can_merge_variant_blocks(left, right))

    def test_can_merge_variant_blocks_right_has_match(self):
        left = self.make_variant_block(0, 4, VARIANT_BASE_DELETION * 3 + VARIANT_BASE_MISMATCH)
        right = self.make_variant_block(4, 6, VARIANT_BASE_MISMATCH + VARIANT_BASE_MATCH)
        self.assertFalse(can_merge_variant_blocks(left, right))

    def test_can_merge_variant_blocks_both_have_match(self):
        left = self.make_variant_block(0, 4, VARIANT_BASE_DELETION * 3 + VARIANT_BASE_MATCH)
        right = self.make_variant_block(4, 6, VARIANT_BASE_MISMATCH + VARIANT_BASE_MATCH)
        self.assertFalse(can_merge_variant_blocks(left, right))

    def test_can_merge_variant_blocks_can_merge(self):
        left = self.make_variant_block(0, 4, VARIANT_BASE_DELETION * 3 + VARIANT_BASE_MISMATCH)
        right = self.make_variant_block(4, 6, VARIANT_BASE_MISMATCH + VARIANT_BASE_INSERTION)
        self.assertTrue(can_merge_variant_blocks(left, right))


class MergeVariantBlocksTestCase(PalamedesBaseCase):
    def test_merge_variant_blocks_simple(self):
        """Test merging 2 single base mismatches"""
        left = VariantBlock(
            Block(0, 1, VARIANT_BASE_MISMATCH),
            [Block(0, 1, "A")],
            [Block(0, 1, "T")],
        )
        right = VariantBlock(
            Block(1, 2, VARIANT_BASE_MISMATCH),
            [Block(1, 2, "G")],
            [Block(1, 2, "C")],
        )

        merged = merge_variant_blocks(left, right)

        self.assertEqual(merged.alignment_block.start, left.alignment_block.start)
        self.assertEqual(merged.alignment_block.end, right.alignment_block.end)
        self.assertEqual(merged.alignment_block.bases, VARIANT_BASE_MISMATCH + VARIANT_BASE_MISMATCH)

        self.assertEqual(len(merged.reference_blocks), 1)
        self.assertEqual(merged.reference_blocks[0].start, left.reference_blocks[0].start)
        self.assertEqual(merged.reference_blocks[0].end, right.reference_blocks[0].end)
        self.assertEqual(
            merged.reference_blocks[0].bases, left.reference_blocks[0].bases + right.reference_blocks[0].bases
        )

        self.assertEqual(len(merged.alternate_blocks), 1)
        self.assertEqual(merged.alternate_blocks[0].start, left.alternate_blocks[0].start)
        self.assertEqual(merged.alternate_blocks[0].end, right.alternate_blocks[0].end)
        self.assertEqual(
            merged.alternate_blocks[0].bases, left.alternate_blocks[0].bases + right.alternate_blocks[0].bases
        )

    def test_merge_variant_blocks_complex(self):
        """Test merging 2 complex VariantBlocks based on the following hypothetical alignment:
         T T T -
         - - C G
        0 1 2 3 4 GLOBAL
        0 1 2 3 . REF
        . . 0 1 2 ALT
        The left side will already have ddm, we are merging in the insertion for the test.
        """
        left = VariantBlock(
            Block(0, 3, "".join([VARIANT_BASE_DELETION, VARIANT_BASE_DELETION, VARIANT_BASE_MISMATCH])),
            [Block(0, 3, "TTT")],
            [Block(0, 1, "C")],
        )
        right = VariantBlock(
            Block(3, 4, VARIANT_BASE_INSERTION),
            [],
            [Block(1, 2, "G")],
        )

        merged = merge_variant_blocks(left, right)

        self.assertEqual(merged.alignment_block.start, left.alignment_block.start)
        self.assertEqual(merged.alignment_block.end, right.alignment_block.end)
        self.assertEqual(
            merged.alignment_block.bases,
            "".join([VARIANT_BASE_DELETION, VARIANT_BASE_DELETION, VARIANT_BASE_MISMATCH, VARIANT_BASE_INSERTION]),
        )

        # there is no right ref block, so everything is on the left one
        self.assertEqual(len(merged.reference_blocks), 1)
        self.assertEqual(merged.reference_blocks[0], left.reference_blocks[0])

        self.assertEqual(len(merged.alternate_blocks), 1)
        self.assertEqual(merged.alternate_blocks[0].start, left.alternate_blocks[0].start)
        self.assertEqual(merged.alternate_blocks[0].end, right.alternate_blocks[0].end)
        self.assertEqual(
            merged.alternate_blocks[0].bases, left.alternate_blocks[0].bases + right.alternate_blocks[0].bases
        )


class GenerateVariantBlocksTestCase(PalamedesBaseCase):
    def test_generate_variant_blocks_all_matches(self):
        alignment = self.make_alignment("A" * 5, "A" * 5)
        self.assertEqual(generate_variant_blocks(alignment), [])

    def test_generate_variant_blocks_all_matches_single_mismatch(self):
        alignment = self.make_alignment("ACT", "AGT")
        variant_blocks = generate_variant_blocks(alignment)

        self.assertEqual(len(variant_blocks), 1)
        self.assertEqual(variant_blocks[0].alignment_block, Block(1, 2, VARIANT_BASE_MISMATCH))
        self.assertEqual(variant_blocks[0].reference_blocks, [Block(1, 2, "C")])
        self.assertEqual(variant_blocks[0].alternate_blocks, [Block(1, 2, "G")])

    def test_generate_variant_blocks_merging(self):
        """
        Test an alignment with a 3 position variant in order, deletion -> mismatch -> insertion
            A T C - T
            A - G A T
           0 1 2 3 4 5 GLOBAL
           0 1 2 3 . 4 REF
           0 . 1 2 3 4 ALT
        """
        alignment = self.make_alignment(
            "ATC-T",
            "A-GAT",
        )
        variant_blocks = generate_variant_blocks(alignment)

        self.assertEqual(len(variant_blocks), 1)
        self.assertEqual(
            variant_blocks[0].alignment_block,
            Block(1, 4, "".join([VARIANT_BASE_DELETION, VARIANT_BASE_MISMATCH, VARIANT_BASE_INSERTION])),
        )
        self.assertEqual(variant_blocks[0].reference_blocks, [Block(1, 3, "TC")])
        self.assertEqual(variant_blocks[0].alternate_blocks, [Block(1, 3, "GA")])

    def test_generate_variant_blocks_complex(self):
        """
        Test an alignment with a number of variants, some merged and some not
            A T C T - - T
            A - C G A A T
           0 1 2 3 4 5 6 7 GLOBAL
           0 1 2 3 4 . . 6 REF
           0 . 1 2 3 4 5 6 ALT

        The full "variant bases" would be: MdMmiiM
        """
        alignment = self.make_alignment(
            "ATCT--T",
            "A-CGAAT",
        )
        variant_blocks = generate_variant_blocks(alignment)

        # mii gets merged so 2 total
        self.assertEqual(len(variant_blocks), 2)
        self.assertEqual(
            variant_blocks[0].alignment_block,
            Block(1, 2, VARIANT_BASE_DELETION),
        )
        self.assertEqual(variant_blocks[0].reference_blocks, [Block(1, 2, "T")])
        self.assertEqual(variant_blocks[0].alternate_blocks, [])

        self.assertEqual(
            variant_blocks[1].alignment_block,
            Block(3, 6, "".join([VARIANT_BASE_MISMATCH, VARIANT_BASE_INSERTION, VARIANT_BASE_INSERTION])),
        )
        self.assertEqual(variant_blocks[1].reference_blocks, [Block(3, 4, "T")])
        self.assertEqual(variant_blocks[1].alternate_blocks, [Block(2, 5, "GAA")])

    def test_generate_variant_blocks_split_subs(self):
        alignment = self.make_alignment(
            "AAAA",
            "TTTT",
        )
        variant_blocks = generate_variant_blocks(alignment, split_consecutive_mismatches=True)

        self.assertEqual(len(variant_blocks), 4)
        for idx, block in enumerate(variant_blocks):
            expected_start = idx
            expecte_end = idx + 1
            self.assertEqual(block.alignment_block, Block(expected_start, expecte_end, VARIANT_BASE_MISMATCH))
            self.assertEqual(block.reference_blocks, [Block(expected_start, expecte_end, "A")])
            self.assertEqual(block.alternate_blocks, [Block(expected_start, expecte_end, "T")])

    def test_generate_variant_blocks_no_split_indels(self):
        """Test split flag does not change how actual indels are treated"""
        alignment = self.make_alignment(
            "AAAA---",
            "---TTTT",
        )
        variant_blocks = generate_variant_blocks(alignment, split_consecutive_mismatches=True)
        self.assertEqual(len(variant_blocks), 1)

        self.assertEqual(
            variant_blocks[0].alignment_block,
            Block(
                0,
                7,
                "".join(
                    [
                        VARIANT_BASE_DELETION,
                        VARIANT_BASE_DELETION,
                        VARIANT_BASE_DELETION,
                        VARIANT_BASE_MISMATCH,
                        VARIANT_BASE_INSERTION,
                        VARIANT_BASE_INSERTION,
                        VARIANT_BASE_INSERTION,
                    ]
                ),
            ),
        )
