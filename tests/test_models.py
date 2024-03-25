from unittest import TestCase

from palamedes.models import Block


class BlockTestCase(TestCase):
    def test_block_create(self):
        bases = "AAAA"
        start = 0
        end = 4
        block = Block(start, end, bases)

        self.assertEqual(block.bases, bases)
        self.assertEqual(block.start, start)
        self.assertEqual(block.end, end)

    def test_block_collapse_empty_error(self):
        with self.assertRaisesRegex(ValueError, "empty list"):
            Block.collapse([])

    def test_block_collapse_not_adjacent_error(self):
        first_block = Block(0, 4, "A" * 4)
        second_block = Block(5, 10, "G" * 6)

        with self.assertRaisesRegex(ValueError, "they must be adjacent."):
            Block.collapse([first_block, second_block])

    def test_block_collapse_single(self):
        block = Block(0, 4, "A" * 4)
        collapsed_block = Block.collapse([block])
        self.assertEqual(collapsed_block, block)

    def test_block_collapse(self):
        first_block = Block(0, 4, "A" * 4)
        second_block = Block(4, 10, "G" * 6)
        third_block = Block(10, 15, "C" * 5)

        collapsed_block = Block.collapse([first_block, second_block, third_block])
        self.assertEqual(collapsed_block.start, first_block.start)
        self.assertEqual(collapsed_block.end, third_block.end)
        self.assertEqual(collapsed_block.bases, "A" * 4 + "G" * 6 + "C" * 5)

    def test_block_collapse_unsorted(self):
        first_block = Block(0, 4, "A" * 4)
        second_block = Block(4, 10, "G" * 6)
        third_block = Block(10, 15, "C" * 5)

        # shuffle the order
        collapsed_block = Block.collapse([third_block, first_block, second_block])
        self.assertEqual(collapsed_block.start, first_block.start)
        self.assertEqual(collapsed_block.end, third_block.end)
        self.assertEqual(collapsed_block.bases, "A" * 4 + "G" * 6 + "C" * 5)
