import logging
from typing import Union

from Bio.Align import Alignment, PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from palamedes.config import (
    ALT_SEQUENCE_ID,
    DEFAULT_EXTEND_GAP_SCORE,
    DEFAULT_MATCH_SCORE,
    DEFAULT_MISMATCH_SCORE,
    DEFAULT_OPEN_GAP_SCORE,
    GLOBAL_ALIGN_MODE,
    REF_SEQUENCE_ID,
)

LOGGER = logging.getLogger(__name__)


def generate_alignment(
    reference_sequence: Union[str, SeqRecord],
    alternate_sequence: Union[str, SeqRecord],
    aligner: PairwiseAligner = None,
) -> Alignment:
    '''
    Using biopython's PairwiseAligner, generate an alignment object representing the best alignment
    between the 2 inputs (which can either be raw string sequences or biopython SeqRecords).

    By default the function creates an aligner object using the defaults, but the caller may provide their
    own pre-configured aligner. This aligner must be set to 'global' mode.

    A warning is generated if a tie is detected between the 1st and 2nd best alignments.
    '''
    if aligner is not None:
        if aligner.mode != GLOBAL_ALIGN_MODE:
            raise ValueError(f'Custom PairwiseAligner must be set to global mode, got: {aligner.mode}')
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
            LOGGER.warning('Tie detected between best and second best alignment!')
    except StopIteration:
        pass

    return best_alignment
