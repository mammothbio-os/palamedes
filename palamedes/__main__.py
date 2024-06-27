import logging
from argparse import ArgumentParser

from Bio.Align import PairwiseAligner

from palamedes import generate_alignment, generate_variant_blocks
from palamedes.align import generate_seq_record
from palamedes.config import MOLECULE_TYPE_PROTEIN, ALT_SEQUENCE_ID, REF_SEQUENCE_ID
from palamedes import __version__
from palamedes.hgvs.utils import categorize_variant_block
from palamedes.hgvs.builders import BUILDER_CONFIG
from palamedes.utils import configure_logging
from palamedes.config import (
    GLOBAL_ALIGN_MODE,
    DEFAULT_MATCH_SCORE,
    DEFAULT_MISMATCH_SCORE,
    DEFAULT_OPEN_GAP_SCORE,
    DEFAULT_EXTEND_GAP_SCORE,
)

LOGGER = logging.getLogger(__name__)


def main() -> None:
    parser = ArgumentParser(
        description="Generate HGVS objects for all variants found in the alignment between 2 sequences",
    )
    parser.add_argument(
        "ref",
        help="Reference sequence",
        type=str,
    )

    parser.add_argument(
        "alt",
        help="Alternate sequence",
        type=str,
    )

    parser.add_argument(
        "--ref-id",
        help="Identifier for reference sequence id",
        type=str,
        default=REF_SEQUENCE_ID,
    )
    parser.add_argument(
        "--alt-id",
        help="Identifier for alternate sequence id",
        type=str,
        default=ALT_SEQUENCE_ID,
    )
    parser.add_argument(
        "--molecule-type",
        help="Molecule type to use",
        choices=(MOLECULE_TYPE_PROTEIN,),
        default=MOLECULE_TYPE_PROTEIN,
    )
    parser.add_argument(
        "--use-non-standard-substitution-rules",
        help=(
            "Flag to enable non-standard substitution rules, "
            "replacing balanced delins calls with multiple substitutions"
        ),
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--match-score",
        help="Match score to use for base alignment",
        type=int,
        default=DEFAULT_MATCH_SCORE,
    )
    parser.add_argument(
        "--mismatch-score",
        help="Mismatch score to use for base alignment",
        type=int,
        default=DEFAULT_MISMATCH_SCORE,
    )
    parser.add_argument(
        "--gap-open-score",
        help="Gap open score to use for base alignment",
        type=int,
        default=DEFAULT_OPEN_GAP_SCORE,
    )
    parser.add_argument(
        "--gap-extend-score",
        help="Gap extend score to use for base alignment",
        type=int,
        default=DEFAULT_EXTEND_GAP_SCORE,
    )
    parser.add_argument(
        "--debug",
        help="Enable debug logging to stderr",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    args = parser.parse_args()
    configure_logging(args.debug)

    LOGGER.debug("Running with args: %s", args)

    aligner = PairwiseAligner(
        mode=GLOBAL_ALIGN_MODE,
        match_score=args.match_score,
        mismatch_score=args.mismatch_score,
        open_gap_score=args.gap_open_score,
        extend_gap_score=args.gap_extend_score,
    )

    ref_seq_record = generate_seq_record(args.ref, args.ref_id, molecule_type=args.molecule_type)
    alt_seq_record = generate_seq_record(args.alt, args.alt_id, molecule_type=args.molecule_type)
    alignment = generate_alignment(
        ref_seq_record,
        alt_seq_record,
        molecule_type=args.molecule_type,
        aligner=aligner,
    )

    LOGGER.debug("Found best alignment with score = %s", getattr(alignment, "score"))
    LOGGER.debug("Alignment:\n%s", str(alignment))

    variant_blocks = generate_variant_blocks(
        alignment,
        split_consecutive_mismatches=args.use_non_standard_substitution_rules,
    )
    LOGGER.debug("%s Variant blocks generated", len(variant_blocks))

    builder = BUILDER_CONFIG[args.molecule_type](alignment)
    for variant_block in variant_blocks:
        category = categorize_variant_block(
            variant_block,
            alignment,
        )
        LOGGER.debug("%s, categorized as: %s", variant_block, category)

        hgvs = builder.build(variant_block, category)
        print(hgvs.format())


if __name__ == "__main__":
    main()
