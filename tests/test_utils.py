from sys import stderr
from unittest import TestCase
from unittest.mock import ANY, patch

from palamedes.utils import configure_logging, contains_repeated_substring, yield_repeating_substrings


class ConfigureLoggingTestCase(TestCase):
    def setUp(self):
        self.logging_patch = patch("palamedes.utils.logging")
        self.logging_mock = self.logging_patch.start()
        self.addCleanup(self.logging_patch.stop)

    def test_configure_logging(self):
        configure_logging()
        self.logging_mock.StreamHandler.assert_called_once_with(stderr)
        self.logging_mock.basicConfig.assert_called_once_with(
            level=self.logging_mock.WARNING,
            format=ANY,
            handlers=[self.logging_mock.StreamHandler()],
        )

    def test_configure_logging_debug(self):
        configure_logging(debug=True)
        self.logging_mock.StreamHandler.assert_called_once_with(stderr)
        self.logging_mock.basicConfig.assert_called_once_with(
            level=self.logging_mock.DEBUG,
            format=ANY,
            handlers=[self.logging_mock.StreamHandler()],
        )


class ContainsRepeatedSubstringTestCase(TestCase):
    def test_contains_repeated_substring_empty_string(self):
        self.assertFalse(contains_repeated_substring(""))

    def test_contains_repeated_substring_single_char(self):
        self.assertFalse(contains_repeated_substring("A"))

    def test_contains_repeated_substring_single_char_repeated(self):
        self.assertTrue(contains_repeated_substring("AA"))

    def test_contains_repeated_substring_multi_char_not_repeated(self):
        self.assertFalse(contains_repeated_substring("ATGATGCATG"))

    def test_contains_repeated_substring_multi_char_repeated(self):
        self.assertTrue(contains_repeated_substring("ATGATGATG"))


class YieldRepeatingSubstringsTestCase(TestCase):
    def test_yield_repeating_substrings_empty_string(self):
        substrings = [_ for _ in yield_repeating_substrings("")]
        self.assertEqual(substrings, [])

    def test_yield_repeating_substrings_single_char(self):
        substrings = [_ for _ in yield_repeating_substrings("A")]
        self.assertEqual(substrings, [])

    def test_yield_repeating_substrings_single_char_repeated_even(self):
        """Single char repeated 4 times, should yield the single and the double"""
        substrings = [_ for _ in yield_repeating_substrings("AAAA")]
        self.assertEqual(substrings, ["A", "AA"])

    def test_yield_repeating_substrings_single_char_repeated_odd(self):
        """Single char repeated 5 times, should yield only the single"""
        substrings = [_ for _ in yield_repeating_substrings("AAAAA")]
        self.assertEqual(substrings, ["A"])

    def test_yield_repeating_substrings_multi_char_not_repeated(self):
        """Multi character string with no repeating units"""
        substrings = [_ for _ in yield_repeating_substrings("ATGCAATGCATAGGACATGACACAC")]
        self.assertEqual(substrings, [])

    def test_yield_repeating_substrings_multi_single_repeat(self):
        """Multi character string with 1 repeating unit (ATG * 3)"""
        substrings = [_ for _ in yield_repeating_substrings("ATG" * 3)]
        self.assertEqual(substrings, ["ATG"])

    def test_yield_repeating_substrings_multi_repeat_even(self):
        """
        Multi character string with several repeating units (even): ATGATGATGATGATGATGATGATG
        - ATG * 8
        - ATGATG * 4
        - ATGATGATGATG * 2
        """
        substrings = [_ for _ in yield_repeating_substrings("ATGATG" * 4)]
        self.assertEqual(substrings, ["ATG", "ATGATG", "ATGATGATGATG"])
