from hgvs.sequencevariant import SequenceVariant

from palamedes.hgvs.builders import HgvsProteinBuilder
from palamedes.models import VariantBlock, Block
from palamedes.config import (
    VARIANT_BASE_MISMATCH,
    HGVS_VARIANT_TYPE_SUBSTITUTION,
    REF_SEQUENCE_ID,
    HGVS_TYPE_PROTEIN,
    HGVS_VARIANT_TYPE_DELETION,
    VARIANT_BASE_DELETION,
    VARIANT_BASE_INSERTION,
    HGVS_VARIANT_TYPE_INSERTION,
    HGVS_VARIANT_TYPE_EXTENSION,
    HGVS_VARIANT_TYPE_REPEAT,
    HGVS_VARIANT_TYPE_DUPLICATION,
    HGVS_VARIANT_TYPE_DELETION_INSERTION,
)

from tests.base import PalamedesBaseCase


class HgvsProteinBuilderTestCase(PalamedesBaseCase):
    def assert_variant_string_matches(self, hgvs_variant: SequenceVariant, variant_string: str) -> None:
        self.assertEqual(
            hgvs_variant.format(conf={"p_3_letter": False}),
            f"{REF_SEQUENCE_ID}:{HGVS_TYPE_PROTEIN}.{variant_string}",
        )

    def test_hgvs_protein_builder_build_substitution(self):
        alignment = self.make_alignment("FFF", "FSF")
        variant_block = VariantBlock(
            Block(1, 2, VARIANT_BASE_MISMATCH),
            [Block(1, 2, "F")],
            [Block(1, 2, "S")],
        )

        as_hgvs = HgvsProteinBuilder(alignment).build(variant_block, HGVS_VARIANT_TYPE_SUBSTITUTION)
        self.assert_variant_string_matches(as_hgvs, "F2S")

    def test_hgvs_protein_builder_build_substitution_first_position(self):
        alignment = self.make_alignment("FFF", "SFF")
        variant_block = VariantBlock(
            Block(0, 1, VARIANT_BASE_MISMATCH),
            [Block(0, 1, "F")],
            [Block(0, 1, "S")],
        )

        as_hgvs = HgvsProteinBuilder(alignment).build(variant_block, HGVS_VARIANT_TYPE_SUBSTITUTION)
        self.assert_variant_string_matches(as_hgvs, "F1S")

    def test_hgvs_protein_builder_build_substitution_last_position(self):
        alignment = self.make_alignment("FFF", "FFS")
        variant_block = VariantBlock(
            Block(2, 3, VARIANT_BASE_MISMATCH),
            [Block(2, 3, "F")],
            [Block(2, 3, "S")],
        )

        as_hgvs = HgvsProteinBuilder(alignment).build(variant_block, HGVS_VARIANT_TYPE_SUBSTITUTION)
        self.assert_variant_string_matches(as_hgvs, "F3S")

    def test_hgvs_protein_builder_build_deletion_single_base(self):
        alignment = self.make_alignment("FFF", "F-F")
        variant_block = VariantBlock(
            Block(1, 2, VARIANT_BASE_DELETION),
            [Block(1, 2, "F")],
            [],
        )

        as_hgvs = HgvsProteinBuilder(alignment).build(variant_block, HGVS_VARIANT_TYPE_DELETION)
        self.assert_variant_string_matches(as_hgvs, "F2del")

    def test_hgvs_protein_builder_build_deletion_multi_base(self):
        """
         F S S S S F
         F - - - - F
        0 1 2 3 4 5 6
        """
        alignment = self.make_alignment("FSSSSF", "F----F")
        variant_block = VariantBlock(
            Block(1, 5, VARIANT_BASE_DELETION * 4),
            [Block(1, 5, "S" * 4)],
            [],
        )

        as_hgvs = HgvsProteinBuilder(alignment).build(variant_block, HGVS_VARIANT_TYPE_DELETION)
        self.assert_variant_string_matches(as_hgvs, "S2_S5del")

    def test_hgvs_protein_builder_build_deletion_multi_base_covering_start(self):
        """
         S S S S F
         - - - - F
        0 1 2 3 4 5
        """
        alignment = self.make_alignment("SSSSF", "----F")
        variant_block = VariantBlock(
            Block(0, 4, VARIANT_BASE_DELETION * 4),
            [Block(0, 4, "S" * 4)],
            [],
        )

        as_hgvs = HgvsProteinBuilder(alignment).build(variant_block, HGVS_VARIANT_TYPE_DELETION)
        self.assert_variant_string_matches(as_hgvs, "S1_S4del")

    def test_hgvs_protein_builder_build_deletion_multi_base_covering_end(self):
        """
         F S S S S
         F - - - -
        0 1 2 3 4 5
        """
        alignment = self.make_alignment("FSSSS", "F----")
        variant_block = VariantBlock(
            Block(1, 5, VARIANT_BASE_DELETION * 4),
            [Block(1, 5, "S" * 4)],
            [],
        )

        as_hgvs = HgvsProteinBuilder(alignment).build(variant_block, HGVS_VARIANT_TYPE_DELETION)
        self.assert_variant_string_matches(as_hgvs, "S2_S5del")

    def test_hgvs_protein_builder_build_insertion_single_base(self):
        alignment = self.make_alignment("F-F", "FSF")
        variant_block = VariantBlock(
            Block(1, 2, VARIANT_BASE_INSERTION),
            [],
            [Block(1, 2, "S")],
        )

        as_hgvs = HgvsProteinBuilder(alignment).build(variant_block, HGVS_VARIANT_TYPE_INSERTION)
        self.assert_variant_string_matches(as_hgvs, "F1_F2insS")

    def test_hgvs_protein_builder_build_insertion_multiple_base(self):
        """
         F - - - - - F
         F S S L Q N F
        0 1 2 3 4 5 6 7
        """
        alignment = self.make_alignment("F-----F", "FSSLQNF")
        variant_block = VariantBlock(
            Block(1, 6, VARIANT_BASE_INSERTION * 5),
            [],
            [Block(1, 6, "SSLQNF")],
        )

        as_hgvs = HgvsProteinBuilder(alignment).build(variant_block, HGVS_VARIANT_TYPE_INSERTION)
        self.assertEqual(
            as_hgvs.format(conf={"p_3_letter": False}), f"{REF_SEQUENCE_ID}:{HGVS_TYPE_PROTEIN}.F1_F2insSSLQNF"
        )

    def test_hgvs_protein_builder_build_extension_upstream_single_base(self):
        """
         - F F
         S F F
        0 1 2 3
        """
        alignment = self.make_alignment("-FF", "SFF")
        variant_block = VariantBlock(
            Block(0, 1, VARIANT_BASE_INSERTION),
            [],
            [Block(0, 1, "S")],
        )

        as_hgvs = HgvsProteinBuilder(alignment).build(variant_block, HGVS_VARIANT_TYPE_EXTENSION)
        self.assert_variant_string_matches(as_hgvs, "F1extS-1")

    def test_hgvs_protein_builder_build_extension_upstream_multiple_base(self):
        """
         - - - F F
         S K L F F
        0 1 2 3 5 6
        """
        alignment = self.make_alignment("---FF", "SKLFF")
        variant_block = VariantBlock(
            Block(0, 3, VARIANT_BASE_INSERTION * 3),
            [],
            [Block(0, 3, "SKL")],
        )

        as_hgvs = HgvsProteinBuilder(alignment).build(variant_block, HGVS_VARIANT_TYPE_EXTENSION)
        self.assertEqual(
            as_hgvs.format(conf={"p_3_letter": False}), f"{REF_SEQUENCE_ID}:{HGVS_TYPE_PROTEIN}.F1extSKL-3"
        )

    def test_hgvs_protein_builder_build_extension_downstream_single_base(self):
        """
         F F -
         F F S
        0 1 2 3
        """
        alignment = self.make_alignment("FF-", "FFS")
        variant_block = VariantBlock(
            Block(2, 3, VARIANT_BASE_INSERTION),
            [],
            [Block(2, 3, "S")],
        )

        as_hgvs = HgvsProteinBuilder(alignment).build(variant_block, HGVS_VARIANT_TYPE_EXTENSION)
        self.assert_variant_string_matches(as_hgvs, "F2extS1")

    def test_hgvs_protein_builder_build_extension_downstream_multiple_bases(self):
        """
         F F - - - -
         F F S L L Q
        0 1 2 3 4 5 6
        """
        alignment = self.make_alignment("FF----", "FFSLLQ")
        variant_block = VariantBlock(
            Block(2, 6, VARIANT_BASE_INSERTION * 4),
            [],
            [Block(2, 6, "SLLQ")],
        )

        as_hgvs = HgvsProteinBuilder(alignment).build(variant_block, HGVS_VARIANT_TYPE_EXTENSION)
        self.assertEqual(
            as_hgvs.format(conf={"p_3_letter": False}), f"{REF_SEQUENCE_ID}:{HGVS_TYPE_PROTEIN}.F2extSLLQ4"
        )

    def test_hgvs_protein_builder_build_duplication_single_base(self):
        """
         F F - A
         F F F A
        0 1 2 3 4
        """
        alignment = self.make_alignment("FF-A", "FFFA")
        variant_block = VariantBlock(
            Block(2, 3, VARIANT_BASE_INSERTION),
            [],
            [Block(2, 3, "F")],
        )

        as_hgvs = HgvsProteinBuilder(alignment).build(variant_block, HGVS_VARIANT_TYPE_DUPLICATION)
        self.assert_variant_string_matches(as_hgvs, "F2dup")

    def test_hgvs_protein_builder_build_duplication_multiple_bases(self):
        """
         F F S A - - - - K
         F F S A F F S A K
        0 1 2 3 4 5 6 7 8 9
        """
        alignment = self.make_alignment("FFSA----K", "FFSAFFSAK")
        variant_block = VariantBlock(
            Block(4, 8, VARIANT_BASE_INSERTION * 4),
            [],
            [Block(4, 8, "FFSA")],
        )

        as_hgvs = HgvsProteinBuilder(alignment).build(variant_block, HGVS_VARIANT_TYPE_DUPLICATION)
        self.assert_variant_string_matches(as_hgvs, "F1_A4dup")

    def test_hgvs_protein_builder_build_repeat_single_base(self):
        """
         F - - A
         F F F A
        0 1 2 3 4
        """
        alignment = self.make_alignment("F--A", "FFFA")
        variant_block = VariantBlock(
            Block(1, 3, VARIANT_BASE_INSERTION * 2),
            [],
            [Block(1, 3, "FF")],
        )

        as_hgvs = HgvsProteinBuilder(alignment).build(variant_block, HGVS_VARIANT_TYPE_REPEAT)
        self.assert_variant_string_matches(as_hgvs, "F1[2]")

    def test_hgvs_protein_builder_build_repeat_multiple_bases(self):
        """
         F L S - - - - - - - - - A
         F L S F L S F L S F L S A
        0 1 2 3 4 5 6 7 8 9 A B C D
        """
        alignment = self.make_alignment("FLS---------A", "FLSFLSFLSFLSA")
        variant_block = VariantBlock(
            Block(3, 12, VARIANT_BASE_INSERTION * 9),
            [],
            [Block(3, 12, "FLSFLSFLS")],
        )

        as_hgvs = HgvsProteinBuilder(alignment).build(variant_block, HGVS_VARIANT_TYPE_REPEAT)
        self.assert_variant_string_matches(as_hgvs, "F1_S3[3]")

    def test_hgvs_protein_builder_build_deletion_insertion(self):
        """
         F L S -
         F K - A
        0 1 2 3 4 ALIGNMENT
        0 1 2 3 . REF
        0 1 . 2 3 ALT
        """
        alignment = self.make_alignment("FLS-", "FK-A")
        variant_block = VariantBlock(
            Block(1, 4, VARIANT_BASE_MISMATCH + VARIANT_BASE_DELETION + VARIANT_BASE_INSERTION),
            [Block(1, 3, "LS")],
            [Block(1, 3, "KA")],
        )

        as_hgvs = HgvsProteinBuilder(alignment).build(variant_block, HGVS_VARIANT_TYPE_DELETION_INSERTION)
        self.assert_variant_string_matches(as_hgvs, "L2_S3delinsKA")
