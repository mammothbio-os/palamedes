from unittest import TestCase

from palamedes.hgvs_utils import categorize_variant_block
from palamedes.config import (
    HGVS_VARIANT_TYPE_SUBSTITUTION,
    HGVS_VARIANT_TYPE_DELETION,
    HGVS_VARIANT_TYPE_EXTENSION,
    HGVS_VARIANT_TYPE_DUPLICATION,
    HGVS_VARIANT_TYPE_REPEAT,
    HGVS_VARIANT_TYPE_INSERTION,
    HGVS_VARIANT_TYPE_DELETION_INSERTION,
    VARIANT_BASE_MATCH,
    VARIANT_BASE_MISMATCH,
    VARIANT_BASE_DELETION,
    VARIANT_BASE_INSERTION,
)
from palamedes.models import Block, VariantBlock


class CategorizeVariantBlockTestCase(TestCase):
    def test_categorize_variant_block_match_error(self):
        variant_block = VariantBlock(
            Block(0, 1, VARIANT_BASE_MATCH),
            [Block(0, 1, "A")],
            [Block(0, 1, "A")],
        )

        with self.assertRaisesRegex(ValueError, "with matching bases"):
            categorize_variant_block(variant_block, "A")

    def test_categorize_variant_block_single_substitution(self):
        variant_block = VariantBlock(
            Block(0, 1, VARIANT_BASE_MISMATCH),
            [Block(0, 1, "A")],
            [Block(0, 1, "T")],
        )

        self.assertEqual(categorize_variant_block(variant_block, "A"), HGVS_VARIANT_TYPE_SUBSTITUTION)

    def test_categorize_variant_block_single_base_deletion(self):
        variant_block = VariantBlock(
            Block(0, 1, VARIANT_BASE_DELETION),
            [Block(0, 1, "A")],
            [],
        )

        self.assertEqual(categorize_variant_block(variant_block, "AAAA"), HGVS_VARIANT_TYPE_DELETION)

    def test_categorize_variant_block_multi_base_deletion(self):
        variant_block = VariantBlock(
            Block(0, 1, VARIANT_BASE_DELETION * 3),
            [Block(0, 3, "AAA")],
            [],
        )

        self.assertEqual(categorize_variant_block(variant_block, "AAAA"), HGVS_VARIANT_TYPE_DELETION)

    def test_categorize_variant_block_upstream_insertion(self):
        variant_block = VariantBlock(
            Block(0, 3, VARIANT_BASE_INSERTION * 3),
            [],
            [Block(0, 3, "AAA")],
        )

        self.assertEqual(categorize_variant_block(variant_block, "---TTTT"), HGVS_VARIANT_TYPE_EXTENSION)

    def test_categorize_variant_block_downstream_insertion(self):
        """
        Assume the following fake alignment:

            T T T T - - -
            - - - - A A A
           0 1 2 3 4 5 6 7 ALIGNMENT
           0 1 2 3 4 . . . REF
           . . . . 0 1 2 3 ALT
        """
        aligned_reference_sequence = "TTTT---"
        variant_block = VariantBlock(
            Block(4, 7, VARIANT_BASE_INSERTION * 3),
            [],
            [Block(0, 3, "AAA")],
        )

        self.assertEqual(
            categorize_variant_block(variant_block, aligned_reference_sequence), HGVS_VARIANT_TYPE_EXTENSION
        )

    def test_categorize_variant_block_duplication(self):
        """
        Using the following alignment:
            T T T T - - - - A T G
           0 1 2 3 4 5 6 7 8 9 A B
        Set the insert to TTTT, which is directly upstream so duplication
        """
        aligned_reference_sequence = "TTTT----ATG"
        variant_block = VariantBlock(
            Block(4, 8, VARIANT_BASE_INSERTION * 4),
            [],
            [Block(0, 4, "TTTT")],
        )

        self.assertEqual(
            categorize_variant_block(variant_block, aligned_reference_sequence), HGVS_VARIANT_TYPE_DUPLICATION
        )

    def test_categorize_variant_block_repeat(self):
        """
        Using the following alignment:
            T T T A - - - - A T G
           0 1 2 3 4 5 6 7 8 9 A B
        Set the insert to TATA, TA  is directly upstream so repeat
        """
        aligned_reference_sequence = "TTTA----ATG"
        variant_block = VariantBlock(
            Block(4, 8, VARIANT_BASE_INSERTION * 4),
            [],
            [Block(0, 4, "TATA")],
        )

        self.assertEqual(categorize_variant_block(variant_block, aligned_reference_sequence), HGVS_VARIANT_TYPE_REPEAT)

    def test_categorize_variant_block_insert(self):
        """
        Using the following alignment:
            T T T A - - - - A T G
           0 1 2 3 4 5 6 7 8 9 A B
        Set the insert to GCTA, just an insertion
        """
        aligned_reference_sequence = "TTTA----ATG"
        variant_block = VariantBlock(
            Block(4, 8, VARIANT_BASE_INSERTION * 4),
            [],
            [Block(0, 4, "GCTA")],
        )

        self.assertEqual(
            categorize_variant_block(variant_block, aligned_reference_sequence), HGVS_VARIANT_TYPE_INSERTION
        )

    def test_categorize_variant_block_double_mismatch(self):
        """2 mismatches is a delins"""
        aligned_reference_sequence = "AAT"
        variant_block = VariantBlock(
            Block(0, 2, VARIANT_BASE_MISMATCH * 2),
            [Block(0, 2, "AA")],
            [Block(0, 2, "CG")],
        )

        self.assertEqual(
            categorize_variant_block(variant_block, aligned_reference_sequence), HGVS_VARIANT_TYPE_DELETION_INSERTION
        )

    def test_categorize_variant_block_mixture(self):
        """
        Using the following alignment:
            A - A T
            - G A T
           0 1 2 3 4
        """
        aligned_reference_sequence = "A-AT"
        variant_block = VariantBlock(
            Block(0, 2, VARIANT_BASE_DELETION + VARIANT_BASE_INSERTION),
            [Block(0, 1, "A")],
            [Block(0, 1, "G")],
        )

        self.assertEqual(
            categorize_variant_block(variant_block, aligned_reference_sequence), HGVS_VARIANT_TYPE_DELETION_INSERTION
        )
