from palamedes import generate_variant_blocks, generate_alignment
from palamedes.config import (
    VARIANT_BASE_INSERTION,
    VARIANT_BASE_DELETION,
    VARIANT_BASE_MISMATCH,
    MOLECULE_TYPE_PROTEIN,
    MOLECULE_TYPE_ANNOTATION_KEY,
    GLOBAL_ALIGN_MODE,
)
from Bio.Align import Alignment, PairwiseAligner
from palamedes.models import Block
from tests.base import PalamedesBaseCase


class GenerateAlignmentTestCase(PalamedesBaseCase):
    def test_generate_alignment_missing_molecule_type_error(self):
        ref, alt = self.make_seq_records("A", "A", molecule_type="foobar")
        del ref.annotations[MOLECULE_TYPE_ANNOTATION_KEY]
        with self.assertRaisesRegex(ValueError, "got: None"):
            generate_alignment(ref, alt)

    def test_generate_alignment_wrong_molecule_type_error(self):
        ref, alt = self.make_seq_records("A", "A", molecule_type="foobar")

        with self.assertRaisesRegex(ValueError, f"expected: {MOLECULE_TYPE_PROTEIN}"):
            generate_alignment(ref, alt)

    def test_generate_alignment_custom_aligner_mode_error(self):
        local_mode = "local"
        custom_aligner = PairwiseAligner(mode=local_mode)
        ref, alt = self.make_seq_records("A", "T")
        with self.assertRaisesRegex(ValueError, f"got: {local_mode}"):
            generate_alignment(ref, alt, aligner=custom_aligner)

    def test_generate_alignment(self):
        ref, alt = self.make_seq_records("A", "T")

        alignment = generate_alignment(ref, alt)

        self.assertTrue(isinstance(alignment, Alignment))
        self.assertIs(alignment.target, ref)
        self.assertIs(alignment.query, alt)

    def test_generate_alignment_custom_aligner(self):
        ref, alt = self.make_seq_records("A", "A")

        custom_match_score = 10_000
        custom_aligner = PairwiseAligner(mode=GLOBAL_ALIGN_MODE, match_score=custom_match_score)

        alignment = generate_alignment(ref, alt, aligner=custom_aligner)
        self.assertTrue(alignment.score, custom_match_score)

    def test_generate_alignment_three_prime_end_most_del(self):
        ref, alt = self.make_seq_records(
            "T" + "A" * 10 + "G",
            "T" + "A" * 9 + "G",
        )

        alignment = generate_alignment(ref, alt)

        # ensure we get the ref unchanged
        # and the last A deleted, not any others
        self.assertEqual(alignment[0], ref.seq)
        self.assertEqual(alignment[1], "T" + "A" * 9 + "-G")

    def test_generate_alignment_three_prime_end_most_double_del(self):
        ref, alt = self.make_seq_records(
            "ATTGCCA",
            "ATGCA",
        )

        alignment = generate_alignment(ref, alt)

        # ensure we get the ref unchanged
        # and the second of each paired based deleted
        self.assertEqual(alignment[0], ref.seq)
        self.assertEqual(alignment[1], "AT-GC-A")

    def test_generate_alignment_three_prime_end_most_ins(self):
        ref, alt = self.make_seq_records(
            "T" + "A" * 9 + "G",
            "T" + "A" * 10 + "G",
        )

        alignment = generate_alignment(ref, alt)

        # ensure we get the alt unchanged
        # and the insertion happens after the last A in the ref
        self.assertEqual(alignment[0], "T" + "A" * 9 + "-G")
        self.assertEqual(alignment[1], alt.seq)

    def test_generate_alignment_three_prime_end_most_double_ins(self):
        ref, alt = self.make_seq_records(
            "ATGCA",
            "ATTGCCA",
        )
        alignment = generate_alignment(ref, alt)

        # ensure we get the alt unchanged
        # and the insertion happens after the last A in the ref
        self.assertEqual(alignment[0], "AT-GC-A")
        self.assertEqual(alignment[1], alt.seq)


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
