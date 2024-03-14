from sys import stdout
from unittest import TestCase
from unittest.mock import ANY, patch

from palamedes.utils import configure_logging


class ConfigureLoggingTests(TestCase):
    def setUp(self):
        self.logging_patch = patch('palamedes.utils.logging')
        self.logging_mock = self.logging_patch.start()
        self.addCleanup(self.logging_patch.stop)

    def test_configure_logging(self):
        configure_logging()
        self.logging_mock.StreamHandler.assert_called_once_with(stdout)
        self.logging_mock.basicConfig.assert_called_once_with(
            level=self.logging_mock.INFO,
            format=ANY,
            handlers=[self.logging_mock.StreamHandler()],
            force=False,
        )

    def test_configure_logging_force(self):
        configure_logging(force=True)
        self.logging_mock.StreamHandler.assert_called_once_with(stdout)
        self.logging_mock.basicConfig.assert_called_once_with(
            level=self.logging_mock.INFO,
            format=ANY,
            handlers=[self.logging_mock.StreamHandler()],
            force=True,
        )
