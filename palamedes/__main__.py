import logging
from argparse import ArgumentParser

from palamedes.align import generate_seq_record, generate_alignment, generate_variant_blocks
from palamedes.config import MOLECULE_TYPE_PROTEIN, ALT_SEQUENCE_ID, REF_SEQUENCE_ID

from palamedes import __version__
from palamedes.hgvs.utils import categorize_variant_block
from palamedes.hgvs.builders import BUILDER_CONFIG
from palamedes.utils import configure_logging

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

    ref_seq_record = generate_seq_record(args.ref, args.ref_id, molecule_type=args.molecule_type)
    alt_seq_record = generate_seq_record(args.alt, args.alt_id, molecule_type=args.molecule_type)
    alignment = generate_alignment(ref_seq_record, alt_seq_record, molecule_type=args.molecule_type)

    LOGGER.debug("Found best alignment with score = %s", getattr(alignment, "score"))
    LOGGER.debug("Alignment:\n%s", str(alignment))

    variant_blocks = generate_variant_blocks(alignment)
    LOGGER.debug("%s Variant blocks generated", len(variant_blocks))

    builder = BUILDER_CONFIG[args.molecule_type](alignment)
    for variant_block in variant_blocks:
        category = categorize_variant_block(variant_block, alignment)
        LOGGER.debug("%s, categorized as: %s", variant_block, category)

        hgvs = builder.build(variant_block, category)
        print(hgvs.format())


if __name__ == "__main__":
    main()
