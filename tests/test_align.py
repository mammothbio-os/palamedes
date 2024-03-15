from unittest import TestCase

from Bio.Align import Alignment, PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from palamedes.align import generate_alignment
from palamedes.config import ALT_SEQUENCE_ID, GLOBAL_ALIGN_MODE, REF_SEQUENCE_ID


class GenerateAlignmentTestCase(TestCase):
    def test_generate_alignment_custom_aligner_mode_error(self):
        local_mode = 'local'
        custom_aligner = PairwiseAligner(mode=local_mode)
        with self.assertRaisesRegex(ValueError, f'got: {local_mode}'):
            generate_alignment('A', 'T', aligner=custom_aligner)

    def test_generate_alignment(self):
        ref_seq = 'A'
        alt_seq = 'T'

        alignment = generate_alignment(ref_seq, alt_seq)
        self.assertTrue(isinstance(alignment, Alignment))

        self.assertEqual(alignment.target.id, REF_SEQUENCE_ID)
        self.assertEqual(str(alignment.target.seq), ref_seq)

        self.assertEqual(alignment.query.id, ALT_SEQUENCE_ID)
        self.assertEqual(str(alignment.query.seq), alt_seq)

    def test_generate_alignment_custom_aligner(self):
        ref_seq = 'A'
        alt_seq = 'A'

        custom_match_score = 10_000
        custom_aligner = PairwiseAligner(mode=GLOBAL_ALIGN_MODE, match_score=custom_match_score)

        alignment = generate_alignment(ref_seq, alt_seq, aligner=custom_aligner)
        self.assertTrue(alignment.score, custom_match_score)

    def test_generate_alignment_input_seq_records(self):
        ref_seq_rec = SeqRecord(Seq('A'), id='custom-ref')
        alt_seq_rec = SeqRecord(Seq('T'), id='custom-alt')

        alignment = generate_alignment(ref_seq_rec, alt_seq_rec)
        self.assertIs(alignment.target, ref_seq_rec)
        self.assertIs(alignment.query, alt_seq_rec)
