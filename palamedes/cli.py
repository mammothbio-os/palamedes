import logging
from argparse import ArgumentParser

from palamedes.align import generate_alignment, generate_variant_blocks
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
    args = parser.parse_args()
    LOGGER.info("Running with args: %s", args)

    alignment = generate_alignment(args.ref, args.alt)
    LOGGER.info("Found best alignment with score = %s", getattr(alignment, "score"))
    LOGGER.info("Alignment:\n%s", str(alignment))

    variant_blocks = generate_variant_blocks(alignment)
    LOGGER.info("%s Variant blocks generated", len(variant_blocks))
    for variant_block in variant_blocks:
        category = categorize_variant_block(alignment, variant_block)
        LOGGER.info("%s, categorized as: %s", variant_block, category)


if __name__ == "__main__":
    main()
