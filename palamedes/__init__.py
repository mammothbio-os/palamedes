from Bio.Align import PairwiseAligner
from Bio.SeqRecord import SeqRecord
from hgvs.sequencevariant import SequenceVariant

from palamedes.config import MOLECULE_TYPE_PROTEIN, ALT_SEQUENCE_ID, REF_SEQUENCE_ID
from palamedes.align import generate_seq_record, generate_alignment, generate_variant_blocks
from palamedes.hgvs.utils import categorize_variant_block
from palamedes.hgvs.builders import BUILDER_CONFIG


__version__ = "0.0.3"


def generate_hgvs_variants(
    reference_sequence: str | SeqRecord,
    alternate_sequence: str | SeqRecord,
    molecule_type: str = MOLECULE_TYPE_PROTEIN,
    aligner: PairwiseAligner | None = None,
) -> list[SequenceVariant]:
    """
    Given the reference and alternate sequences, as either strings or `Bio.SeqRecord.SeqRecord` objects, compute the
    alignment, find the variants and build `hgvs.sequencevariant.SequenceVariant` objects from them. Returning the list
    of objects at the end. Note the following things:

    - Variants and their HGVS objects are based on consecutive alignment positions that are non matches, any match will split or disrupt the variant.
    - Variant and their HGVS objects are handled in isolation from others, they do not "account" for things that happen upstream of their position in the alignment.

    A custom `PairwiseAligner` may be provided to tweak the alignment score, but the default one uses these scores:

    - mode: `global` (this is required, even for a custom one)
    - match score: `1`
    - mismatch score: `-1`
    - open gap score: `-1`
    - extend gap score: `-0.1`

    If using pre-built `SeqRecord` objects, be sure to set the `molecule_type` annotation key to a supported molecule type
    and pass in the corresponding molecule_type as an input. At the time of writing, only `protein` is supported as a
    molecule type.

    Example using a raw string:

    .. code-block:: python

        >>> from palamedes import generate_hgvs_variants
        >>> generate_hgvs_variants("PFKISIHL", "TPFKISIH")
        [
            SequenceVariant(ac=ref, type=p, posedit=Pro1extThr-1, gene=None),
            SequenceVariant(ac=ref, type=p, posedit=Leu8del, gene=None),
        ]

    Example using pre-configured `SeqRecord` objects:

    .. code-block:: python

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> from palamedes import generate_hgvs_variants
        >>> ref = SeqRecord(Seq("PFKISIHL"), id="Jelleine-I", annotations={"molecule_type": "protein"})
        >>> alt = SeqRecord(Seq("TPFKISIH"), id="Jelleine-IV", annotations={"molecule_type": "protein"})
        >>> generate_hgvs_variants(ref, alt)
        [
            SequenceVariant(ac=Jelleine-I, type=p, posedit=Pro1extThr-1, gene=None),
            SequenceVariant(ac=Jelleine-I, type=p, posedit=Leu8del, gene=None)
        ]
    """

    if molecule_type not in BUILDER_CONFIG:
        raise NotImplementedError(f"No HGVS builder is defined for molecule_type: {molecule_type}!")

    ref_seq_record = (
        generate_seq_record(reference_sequence, REF_SEQUENCE_ID, molecule_type=molecule_type)
        if isinstance(reference_sequence, str)
        else reference_sequence
    )

    alt_seq_record = (
        generate_seq_record(alternate_sequence, ALT_SEQUENCE_ID, molecule_type=molecule_type)
        if isinstance(alternate_sequence, str)
        else alternate_sequence
    )

    alignment = generate_alignment(ref_seq_record, alt_seq_record, molecule_type=molecule_type, aligner=aligner)
    variant_blocks = generate_variant_blocks(alignment)
    builder = BUILDER_CONFIG[molecule_type](alignment)
    return [
        builder.build(
            variant_block,
            categorize_variant_block(variant_block, alignment),
        )
        for variant_block in variant_blocks
    ]
