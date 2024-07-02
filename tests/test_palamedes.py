from palamedes import generate_alignment, generate_hgvs_variants_from_alignment
from palamedes.config import (
    MOLECULE_TYPE_PROTEIN,
    MOLECULE_TYPE_ANNOTATION_KEY,
    GLOBAL_ALIGN_MODE,
)
from Bio.Align import Alignment, PairwiseAligner
from tests.base import PalamedesBaseCase
from tests.hgvs.test_builders import HgvsProteinBuilderTestCase


class GenerateHGVSVariantsFromAlignmentTestCase(HgvsProteinBuilderTestCase):
    def test_generate_hgvs_variants_from_alignment(self):
        ref, alt = self.make_seq_records(
            "ATGCA",
            "ATTGCCA",
        )
        alignment = generate_alignment(ref, alt)

        # ensure we get the alt unchanged
        # and the insertion happens after the last A in the ref
        self.assertEqual(alignment[0], "AT-GC-A")
        self.assertEqual(alignment[1], alt.seq)

        variant_blocks = generate_hgvs_variants_from_alignment(alignment)

        self.assert_variant_string_matches(
            variant_blocks[0],
            "T2dup",
        )

        self.assert_variant_string_matches(
            variant_blocks[1],
            "C4dup",
        )

    def test_generate_hgvs_variants_from_alignment_molecule_type_error(self):
        bad_molecule_type = "FAKE"
        with self.assertRaisesRegex(
            NotImplementedError, f"No HGVS builder is defined for molecule_type: {bad_molecule_type}!"
        ):
            ref, alt = self.make_seq_records(
                "ATGCA",
                "ATTGCCA",
            )
            alignment = generate_alignment(ref, alt)
            generate_hgvs_variants_from_alignment(alignment, molecule_type=bad_molecule_type)


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
