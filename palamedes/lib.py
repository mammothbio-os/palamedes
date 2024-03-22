from typing import List, Union
from Bio.Align import Alignment, PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from hgvs.sequencevariant import SequenceVariant

from palamedes.config import MOLECULE_TYPE_PROTEIN

def generate_variants(
    reference_sequence: Union[str, SeqRecord],
    alternate_sequence: Union[str, SeqRecord],
    molecule_type=MOLECULE_TYPE_PROTEIN,
    aligner=None,
) -> List[]

    
