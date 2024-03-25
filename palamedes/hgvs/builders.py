from Bio.Align import Alignment
from hgvs.edit import (
    Repeat,
    Dup,
    AAExt,
    AARefAlt,
    AASub,
)
from hgvs.posedit import PosEdit
from hgvs.location import Interval, AAPosition
from hgvs.sequencevariant import SequenceVariant

from palamedes.align import get_upstream_reference_sequence
from palamedes.models import VariantBlock
from palamedes.config import (
    HGVS_VARIANT_TYPE_SUBSTITUTION,
    HGVS_VARIANT_TYPE_DELETION,
    HGVS_VARIANT_TYPE_EXTENSION,
    HGVS_VARIANT_TYPE_DUPLICATION,
    HGVS_VARIANT_TYPE_REPEAT,
    HGVS_VARIANT_TYPE_INSERTION,
    HGVS_VARIANT_TYPE_DELETION_INSERTION,
    HGVS_TYPE_PROTEIN,
    MOLECULE_TYPE_PROTEIN,
)
from palamedes.utils import yield_repeating_substrings, zbho_to_obfc, zb_to_ob, zb_position_to_end_coordinate


class HgvsProteinBuilder:
    def __init__(self, alignment: Alignment) -> None:
        self._alignment = alignment

    def build(self, variant_block: VariantBlock, hgvs_type: str) -> SequenceVariant:
        pos_edit_builder_funcs = {
            HGVS_VARIANT_TYPE_SUBSTITUTION: self._build_substitution,
            HGVS_VARIANT_TYPE_DELETION: self._build_deletion,
            HGVS_VARIANT_TYPE_INSERTION: self._build_insertion,
            HGVS_VARIANT_TYPE_EXTENSION: self._build_extension,
            HGVS_VARIANT_TYPE_DUPLICATION: self._build_duplication,
            HGVS_VARIANT_TYPE_REPEAT: self._build_repeat,
            HGVS_VARIANT_TYPE_DELETION_INSERTION: self._build_deletion_insertion,
        }
        pos_edit = pos_edit_builder_funcs[hgvs_type](variant_block)

        return SequenceVariant(ac=self._alignment.target.id, type=HGVS_TYPE_PROTEIN, posedit=pos_edit)

    def _build_substitution(self, variant_block: VariantBlock) -> PosEdit:
        """
        Protein substitution build logic, this is a fairly simple case since the variant block has ref + alt data and
        the variant is always a single position. Convert the coordinates to OBFC and return the PosEdit.
        """
        reference_block = variant_block.reference_blocks[0]
        alternate_block = variant_block.alternate_blocks[0]
        start_obfc, end_obfc = zbho_to_obfc(reference_block.start, reference_block.end)

        return PosEdit(
            pos=Interval(
                start=AAPosition(
                    base=start_obfc,
                    aa=reference_block.bases,
                ),
                end=AAPosition(
                    base=end_obfc,
                    aa=reference_block.bases,
                ),
            ),
            edit=AASub(
                ref=reference_block.bases,
                alt=alternate_block.bases,
            ),
        )

    def _build_deletion(self, variant_block: VariantBlock) -> PosEdit:
        """
        Protein deletion build logic, this is a fairly simple case since the variant block has the ref data,
        including the positions + sequence that was deleted. Convert the coordinates to OBFC and return the PosEdit.
        """
        reference_block = variant_block.reference_blocks[0]
        start_obfc, end_obfc = zbho_to_obfc(reference_block.start, reference_block.end)

        return PosEdit(
            pos=Interval(
                start=AAPosition(
                    base=start_obfc,
                    aa=reference_block.bases[0],
                ),
                end=AAPosition(
                    base=end_obfc,
                    aa=reference_block.bases[-1],
                ),
            ),
            edit=AARefAlt(
                ref=reference_block.bases,
                alt=None,
            ),
        )

    def _build_insertion(self, variant_block: VariantBlock) -> PosEdit:
        """
        Protein insertion build logic, this is a more complicated example since the ref data does not exist on the
        variant block and has to be looked up from the Alignment. We leverage the .indices field on the Alignment
        which maps alignment coordinates back to sequence coordinates (gaps have -1). This is a numpy array so
        all values are cast to python integers before being saved into an HGVS object.

        The Interval we want is: (1 base upstream of the insert, 1 base downstream of the insert). To get there we
        do the following:
        - Get the start position of the insert in the alignment. This start is the index for that position within the
          .indices array. Subtract 1 from that to get the index of the next base upstream (the flanking ref base). This
          could out of bounds if the variant has not been categorized properly (should be an extension in this case).
        - Lookup the ref start index from .indices using the computed index above
        - Lookup the ref end index from .indices using the end position of the insert in the alignment. No offset is
          needed here since the end of a ZBHO interval is the same as the index for the next base.
        - Convert the 2 zero based indices into one based indices.
        - Build and return the object, using the zero based indices to get the anchor bases from the ref sequence
        """
        upstream_ref_base_index = variant_block.alignment_block.start - 1
        ref_flanking_start_position = int(self._alignment.indices[0][upstream_ref_base_index])
        ref_flanking_end_position = int(self._alignment.indices[0][variant_block.alignment_block.end])

        ref_flanking_start_position_ob = zb_to_ob(ref_flanking_start_position)
        ref_flanking_end_position_ob = zb_to_ob(ref_flanking_end_position)

        return PosEdit(
            pos=Interval(
                start=AAPosition(
                    base=ref_flanking_start_position_ob,
                    aa=self._alignment.target.seq[ref_flanking_start_position],
                ),
                end=AAPosition(
                    base=ref_flanking_end_position_ob,
                    aa=self._alignment.target.seq[ref_flanking_end_position],
                ),
            ),
            edit=AARefAlt(
                ref=None,
                alt=variant_block.alternate_blocks[0].bases,
            ),
        )

    def _build_extension(self, variant_block: VariantBlock) -> PosEdit:
        """
        Protein extension build logic. This case is simple since the Interval is always the first or the last base.
        Check if it is a start or end extension and get the first or last reference base from the ref sequence. The
        position is either 1 or the last base (len of the ref sequence). Since the values are static, they are defined
        in OB already.
        """
        is_start = variant_block.alignment_block.start == 0
        ref_base = self._alignment.target.seq[0] if is_start else self._alignment.target.seq[-1]
        ref_position = 1 if is_start else len(self._alignment.target.seq)

        return PosEdit(
            pos=Interval(
                start=AAPosition(
                    base=ref_position,
                    aa=ref_base,
                ),
                end=AAPosition(
                    base=ref_position,
                    aa=ref_base,
                ),
            ),
            edit=AAExt(
                ref=ref_base,
                aaterm=variant_block.alternate_blocks[0].bases,
                length=len(variant_block.alternate_blocks[0].bases) * (-1 if is_start else 1),
            ),
        )

    def _build_duplication(self, variant_block: VariantBlock) -> PosEdit:
        """
        Protein duplication build logic. This is another complicated case that requires some index munging with the
        reference sequence. It is similar to the insertion except that we want the duplicated bases/positions, not
        the bases/positions that flank the insertion. So we get the index for the first upstream base and convert it
        to an end coordinate. We then subtract the length of the insertion from the end to get the start. Finally,
        convert to OBFC. The edit object is just an empty Dup().
        """
        upstream_ref_base_index = variant_block.alignment_block.start - 1
        ref_duplication_end_position = zb_position_to_end_coordinate(
            self._alignment.indices[0][upstream_ref_base_index]
        )
        ref_duplication_start_position = ref_duplication_end_position - len(variant_block.alternate_blocks[0].bases)
        ref_duplication_start_position_obfc, ref_duplication_end_position_obfc = zbho_to_obfc(
            ref_duplication_start_position, ref_duplication_end_position
        )

        return PosEdit(
            pos=Interval(
                start=AAPosition(
                    base=ref_duplication_start_position_obfc,
                    aa=variant_block.alternate_blocks[0].bases[0],
                ),
                end=AAPosition(
                    base=ref_duplication_end_position_obfc,
                    aa=variant_block.alternate_blocks[0].bases[-1],
                ),
            ),
            edit=Dup(),
        )

    def _build_repeat(self, variant_block: VariantBlock) -> PosEdit:
        """
        Protein repeat build logic. This is the most complicated case, and it builds on the duplication case. The
        key difference is that instead of taking the insert sequence as the matching upstream element, we first
        have to find the largest repeating sub-string that matches the directly upstream bases. Once that is found,
        the interval covering the reference bases being repeated are found using the same logic as with duplication, but
        using the sub-string for the length. The last thing is to find the repeat number, which can be computed based on
        integer division between the length of the insert and the length of the largest repeat.
        """
        largest_upstream_repeat = [
            substring
            for substring in yield_repeating_substrings(variant_block.alternate_blocks[0].bases)
            if get_upstream_reference_sequence(self._alignment, variant_block.alignment_block.start, len(substring))
            == substring
        ][-1]

        upstream_ref_base_index = variant_block.alignment_block.start - 1
        ref_duplication_end_position = zb_position_to_end_coordinate(
            self._alignment.indices[0][upstream_ref_base_index]
        )
        ref_duplication_start_position = ref_duplication_end_position - len(largest_upstream_repeat)
        start_obfc, end_obfc = zbho_to_obfc(ref_duplication_start_position, ref_duplication_end_position)

        repeat_value = len(variant_block.alternate_blocks[0].bases) // len(largest_upstream_repeat)

        return PosEdit(
            pos=Interval(
                start=AAPosition(
                    base=start_obfc,
                    aa=largest_upstream_repeat[0],
                ),
                end=AAPosition(
                    base=end_obfc,
                    aa=largest_upstream_repeat[-1],
                ),
            ),
            edit=Repeat(min=repeat_value, max=repeat_value),
        )

    def _build_deletion_insertion(self, variant_block: VariantBlock) -> PosEdit:
        """
        Protein delins logic. Fairly simple logic by taking the ref + alt blocks and combining them together.
        """
        reference_block = variant_block.reference_blocks[0]
        alternate_block = variant_block.alternate_blocks[0]

        start_obfc, end_obfc = zbho_to_obfc(reference_block.start, reference_block.end)
        return PosEdit(
            pos=Interval(
                start=AAPosition(
                    base=start_obfc,
                    aa=reference_block.bases[0],
                ),
                end=AAPosition(
                    base=end_obfc,
                    aa=reference_block.bases[-1],
                ),
            ),
            edit=AARefAlt(
                ref=reference_block.bases,
                alt=alternate_block.bases,
            ),
        )


BUILDER_CONFIG = {
    MOLECULE_TYPE_PROTEIN: HgvsProteinBuilder,
}
