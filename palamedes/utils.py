import logging
import sys
from typing import Iterable


def configure_logging(debug: bool = False) -> None:
    stderr_handler = logging.StreamHandler(sys.stderr)
    logging.basicConfig(
        level=logging.DEBUG if debug else logging.WARNING,
        format="[%(asctime)s] {%(filename)s:%(lineno)d} %(levelname)s - %(message)s",
        handlers=[stderr_handler],
    )


def contains_repeated_substring(input_string: str) -> bool:
    """
    Common leetcode problem to determine if an input_string has a repeating sub-string which
    can be repeated 1 or more times to recreate the input string. This leverages the fact
    that such a string is a non-trivial rotation of itself and a string A is a rotation of
    another string B if A is a sub-string of B+B.

    We use find on B+B to see the location of the sub-string. It cannot be -1 (not found) or
    the length of the input sequence (since that is just the double).

    This is adapted from:
    https://stackoverflow.com/questions/55823298/how-do-i-check-if-a-string-is-entirely-made-of-the-same-substring
    """
    doubled = input_string + input_string
    return doubled.find(input_string, 1) not in {-1, len(input_string)}


def yield_repeating_substrings(input_string: str) -> Iterable[str]:
    """
    Generator to yield sub-strings of the input sequence which can be repeated some number of times
    to regenerate the input sequence. This is similar in nature to contains_repeated_substring, but
    the actual sub-strings are returned. The logic is like so:
    - Start with a length 1 sub-string (the first character)
    - Check if 1 (the length) is a divisor of the input string length (remainder of divmod == 0)
    - If so, check if the sub-string repeated q times (the quotient) equals the input_string (and yield if so)
    - Repeat the above with larger and larger sub-strings (starting from zero, but adding a new character each time)
    - Stop when you sub-string doubled in length would be larger then the input sequence
    """
    input_len = len(input_string)
    idx = 0

    while ((idx + 1) * 2) <= input_len:
        substr_for_idx = input_string[0 : idx + 1]
        quotient, remainder = divmod(input_len, len(substr_for_idx))
        if remainder == 0 and substr_for_idx * quotient == input_string:
            yield substr_for_idx

        idx += 1


def zb_to_ob(zb: int) -> int:
    """
    Convert a single zero based position (or sequence index) into one based. This means simply adding 1 to it.
    This function simply exists to make things more explicit and avoid having magic "+ 1"s everywhere.
    """
    return zb + 1


def zbho_to_obfc(zbho_start: int, zbho_end: int) -> tuple[int, int]:
    """
    Convert a zero based half open (start, end) into a one based fully closed, based on the following:
    - add 1 to both positions to account for zero ==> one
    - subtract 1 to the end position to account for half open to fully closed
    """
    return (
        zb_to_ob(zbho_start),
        zb_to_ob(zbho_end) - 1,
    )


def zb_position_to_end_coordinate(zb_position: int) -> int:
    """
    Convert a zero based position (or index) into the end coordinate of a ZBHO interval including
    this position as the ending position. This is done by simply adding 1 to the position to make it
    half open. Again this is just to make this clearer and more explicit.
    """
    return zb_position + 1
