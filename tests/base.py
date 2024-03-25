from unittest import TestCase

from Bio.Align import Alignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from palamedes.config import (
    MOLECULE_TYPE_PROTEIN,
    MOLECULE_TYPE_ANNOTATION_KEY,
    REF_SEQUENCE_ID,
    ALT_SEQUENCE_ID,
    ALIGNMENT_GAP_CHAR,
)


class PalamedesBaseCase(TestCase):
    def make_seq_records(
        self, ref_seq: str, alt_seq: str, molecule_type: str = MOLECULE_TYPE_PROTEIN
    ) -> tuple[SeqRecord, SeqRecord]:
        ref = SeqRecord(Seq(ref_seq), id=REF_SEQUENCE_ID, annotations={MOLECULE_TYPE_ANNOTATION_KEY: molecule_type})
        alt = SeqRecord(Seq(alt_seq), id=ALT_SEQUENCE_ID, annotations={MOLECULE_TYPE_ANNOTATION_KEY: molecule_type})
        return ref, alt

    def make_alignment(self, ref_aligned_bases: str, alt_aligned_bases: str) -> Alignment:
        coords = Alignment.infer_coordinates([ref_aligned_bases, alt_aligned_bases])
        ref, alt = self.make_seq_records(
            ref_aligned_bases.replace(ALIGNMENT_GAP_CHAR, ""), alt_aligned_bases.replace(ALIGNMENT_GAP_CHAR, "")
        )
        return Alignment([ref, alt], coords)
