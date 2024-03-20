import logging
from argparse import ArgumentParser

from palamedes.align import generate_seq_records, generate_alignment, generate_variant_blocks
from palamedes.config import MOLECULE_TYPE_PROTEIN
from palamedes.hgvs_utils import categorize_variant_block
from palamedes.utils import configure_logging

LOGGER = logging.getLogger(__name__)


def main() -> None:
    configure_logging()
    parser = ArgumentParser(description="Generate a variant summary string for an alignment between 2 sequences")
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
        "--molecule-type",
        help="Molecule type to use",
        choices=(MOLECULE_TYPE_PROTEIN,),
        default=MOLECULE_TYPE_PROTEIN,
    )

    args = parser.parse_args()
    LOGGER.info("Running with args: %s", args)

    ref_seq_record, alt_seq_record = generate_seq_records(args.ref, args.alt, molecule_type=args.molecule_type)
    alignment = generate_alignment(ref_seq_record, alt_seq_record, molecule_type=args.molecule_type)
    LOGGER.info("Found best alignment with score = %s", getattr(alignment, "score"))
    LOGGER.info("Alignment:\n%s", str(alignment))

    variant_blocks = generate_variant_blocks(alignment)
    LOGGER.info("%s Variant blocks generated", len(variant_blocks))
    for variant_block in variant_blocks:
        category = categorize_variant_block(variant_block, alignment)
        LOGGER.info("%s, categorized as: %s", variant_block, category)


if __name__ == "__main__":
    main()
