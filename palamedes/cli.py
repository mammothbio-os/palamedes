import logging
from argparse import ArgumentParser

from palamedes.utils import configure_logging

LOGGER = logging.getLogger(__name__)


def main():
    configure_logging()
    parser = ArgumentParser(description='Generate a variant summary string for an alignment between 2 sequences')
    parser.add_argument(
        'ref',
        help='Reference sequence',
        type=str,
    )

    parser.add_argument(
        'alt',
        help='Alternate sequence',
        type=str,
    )
    args = parser.parse_args()
    LOGGER.info('Running with args: %s', args)


if __name__ == '__main__':
    main()
