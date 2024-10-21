from Bio.Align import PairwiseAligner, Alignment
from Bio.SeqRecord import SeqRecord
from hgvs.sequencevariant import SequenceVariant

from palamedes.align import generate_seq_record, generate_variant_blocks, reverse_seq_record
from palamedes.hgvs.utils import categorize_variant_block
from palamedes.hgvs.builders import BUILDER_CONFIG
from palamedes.config import (
    DEFAULT_EXTEND_GAP_SCORE,
    DEFAULT_MATCH_SCORE,
    DEFAULT_MISMATCH_SCORE,
    DEFAULT_OPEN_GAP_SCORE,
    GLOBAL_ALIGN_MODE,
    MOLECULE_TYPE_PROTEIN,
    ALT_SEQUENCE_ID,
    REF_SEQUENCE_ID,
)

__version__ = "0.0.9"


def generate_hgvs_variants_from_alignment(
    alignment: Alignment, use_non_standard_substitution_rules: bool = False, molecule_type: str = MOLECULE_TYPE_PROTEIN
) -> list[SequenceVariant]:
    """
    Given a pairwise alignment object and a molecule type, generate a list of HGVS SequenceVariants.

    - Alignment: Generated via BioPython.PairwiseAligner - (See `generate_alignment` for more information)

    - molecule_type: Currently only molecule type 'protein' is supported.

    - An optional flag: `use_non_standard_substitution_rules` is a boolean flag which will enable logic that treats multiple consecutive mismatches as separate subsitutions, vs merging together into a delins. This is against HGVS
    spec but has utility for some use cases.

    .. code-block:: python

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> from palamedes import generate_hgvs_variants_from_alignment, generate_alignment
        >>> ref = SeqRecord(Seq("PFKISIHL"), id="Jelleine-I", annotations={"molecule_type": "protein"})
        >>> alt = SeqRecord(Seq("TPFKISIH"), id="Jelleine-IV", annotations={"molecule_type": "protein"})
        >>> alignment = generate_alignment(ref, alt)
        >>> generate_hgvs_variants_from_alignment(alignment)
        [
            SequenceVariant(ac=Jelleine-I, type=p, posedit=Pro1extThr-1, gene=None),
            SequenceVariant(ac=Jelleine-I, type=p, posedit=Leu8del, gene=None)
        ]
    """
    if molecule_type not in BUILDER_CONFIG:
        raise NotImplementedError(f"No HGVS builder is defined for molecule_type: {molecule_type}!")

    variant_blocks = generate_variant_blocks(
        alignment,
        split_consecutive_mismatches=use_non_standard_substitution_rules,
    )
    builder = BUILDER_CONFIG[molecule_type](alignment)
    return [
        builder.build(
            variant_block,
            categorize_variant_block(variant_block, alignment),
        )
        for variant_block in variant_blocks
    ]


def generate_alignment(
    reference_seq_record: SeqRecord,
    alternate_seq_record: SeqRecord,
    molecule_type: str = MOLECULE_TYPE_PROTEIN,
    aligner: PairwiseAligner | None = None,
) -> Alignment:
    """
    Using biopython's PairwiseAligner, generate an alignment object representing the best alignment
    between the 2 biopython SeqRecords. Note that the molecule_type argument will be used to specify a molecule_type
    annotation on the SeqRecord.

    By default the function creates an aligner object using the defaults, but the caller may provide their
    own pre-configured aligner. This aligner must be set to 'global' mode.

    Note that it is possible for multiple alignments to be returned with the same max score. Synthetic testing
    has shown that generally the highest scoring alignments are returned in "left to right" order. Take the
    following example, which was for testing HGVS duplications. This was the "intended" alignment:

        ATC---GGGGGGGG
        ATCATCGGGGGGGG

    However, there are multiple ways to get a three base insertion, all with the same score:
        ipdb> for aln in aligner.align('ATCGGGGGGGG', 'ATCATCGGGGGGGG'):
        print(aln)

        target            0 ---ATCGGGGGGGG 11
                          0 ---||||||||||| 14
        query             0 ATCATCGGGGGGGG 14

        target            0 A---TCGGGGGGGG 11
                          0 |---|||||||||| 14
        query             0 ATCATCGGGGGGGG 14

        target            0 AT---CGGGGGGGG 11
                          0 ||---||||||||| 14
        query             0 ATCATCGGGGGGGG 14

        target            0 ATC---GGGGGGGG 11
                          0 |||---|||||||| 14
        query             0 ATCATCGGGGGGGG 14

    Due to the HGVS rule of always indicating that the 3' most base is considered the modified one, a best effort
    attempt is made to return the most "right aligned" alignment, by returning the last alignment with the highest
    score. This way not yield the ideal results in more complicated cases.

    .. code-block:: python

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> from palamedes import generate_alignment
        >>> ref = SeqRecord(Seq("PFKISIHL"), id="Jelleine-I", annotations={"molecule_type": "protein"})
        >>> alt = SeqRecord(Seq("TPFKISIH"), id="Jelleine-IV", annotations={"molecule_type": "protein"})
        >>> generate_alignment(ref, alt)
        <Alignment object (2 rows x 9 columns) at ...>
    """
    if aligner is not None:
        if aligner.mode != GLOBAL_ALIGN_MODE:
            raise ValueError(f"Custom PairwiseAligner must be set to global mode, got: {aligner.mode}")
    else:
        aligner = PairwiseAligner(
            mode=GLOBAL_ALIGN_MODE,
            match_score=DEFAULT_MATCH_SCORE,
            mismatch_score=DEFAULT_MISMATCH_SCORE,
            open_gap_score=DEFAULT_OPEN_GAP_SCORE,
            extend_gap_score=DEFAULT_EXTEND_GAP_SCORE,
        )

    if (ref_molecule_type := reference_seq_record.annotations.get("molecule_type")) != molecule_type:
        raise ValueError(
            "Cannot generate alignment, reference_seq_record is a SeqRecord an invalid molecule_type annotation "
            f"got: {ref_molecule_type}, expected: {molecule_type}!"
        )

    if (alt_molecule_type := alternate_seq_record.annotations.get("molecule_type")) != molecule_type:
        raise ValueError(
            "Cannot generate alignment, alternate_seq_record is a SeqRecord an invalid molecule_type annotation "
            f"got: {alt_molecule_type}, expected: {molecule_type}!"
        )

    # reverse the sequences and then align the reversed, keeping the first best alignment
    reversed_ref_seq_record = reverse_seq_record(reference_seq_record)
    reversed_alt_seq_record = reverse_seq_record(alternate_seq_record)
    reversed_alignments = aligner.align(reversed_ref_seq_record, reversed_alt_seq_record)
    reversed_alignment = reversed_alignments[0]

    # undo the reversal, to recover the "last" highest scoring alignment for the forward
    # which should correspond to the 3' end most alignment and follow HGSV spec
    # note 'infer_coordinates' will soon be deprecated, update to 'parse_printed_alignment' once its live
    forward_coordinates = Alignment.infer_coordinates(
        [
            reversed_alignment[0][::-1],
            reversed_alignment[1][::-1],
        ]
    )
    forward_alignment = Alignment(
        [reference_seq_record, alternate_seq_record],
        forward_coordinates,
    )

    # score is not technically an attribute on the class
    setattr(forward_alignment, "score", reversed_alignment.score)

    return forward_alignment


def generate_hgvs_variants(
    reference_sequence: str | SeqRecord,
    alternate_sequence: str | SeqRecord,
    molecule_type: str = MOLECULE_TYPE_PROTEIN,
    aligner: PairwiseAligner | None = None,
    use_non_standard_substitution_rules: bool = False,
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

    An optional flag: `use_non_standard_substitution_rules` can be passed as `True`, which will enable logic that treats
    multiple consecutive mismatches as separate subsitutions, vs merging together into a delins. This is against HGVS
    spec but has utility for some use cases.

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
        raise NotImplementedError(
            f"Type {molecule_type} unsupported or unrecognized! Current supported types: {','.join(BUILDER_CONFIG.keys())}"
        )

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
    return generate_hgvs_variants_from_alignment(alignment, use_non_standard_substitution_rules, molecule_type)
