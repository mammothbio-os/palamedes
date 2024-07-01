# palamedes

![Merge - Passing](https://github.com/mammothbio-os/palamedes/actions/workflows/merge.yaml/badge.svg) ![Release - Passing](https://github.com/mammothbio-os/palamedes/actions/workflows/release.yaml/badge.svg) ![PyPI - Version](https://img.shields.io/pypi/v/palamedes) ![Read The Docs - Version](https://readthedocs.org/projects/mammothbio-os-palamedes/badge/?version=stable)

This repo contains a python package and CLI entrypoint which can be used to generate a list of [HGVS](https://github.com/biocommons/hgvs) variants representing the difference between 2 sequences, using a global alignment. The idea is to leverage the HGVS spec for more applications, since it provides a solid framework for maintaining a consistent set of rules and logic
around variants.

## Documentation

Documentation for the project can be found [here](https://mammothbio-os-palamedes.readthedocs.io/en/stable/)

## Installing

Palamedes uses the [hgvs](https://github.com/biocommons/hgvs) package as a dependency. At this time, `hgvs` requires `postgresql` system dependencies to be installed before use. Please see the [README](https://github.com/biocommons/hgvs/blob/main/README.md#installing-hgvs-locally) for install instructions if needed.

Palamedes itself is packaged in PyPI, to install simply run: `pip install palamedes`

## Usage - CLI

Palamedes includes a CLI entrypoint, which is mostly useful for debugging and exploration. Once installed, simply run `palamedes` and provide a reference and alternate sequence:
```shell
palamedes PFKISIHL TPFKISIH
[2024-03-25 11:40:59,904] {cli.py:39} INFO - Running with args: Namespace(ref='PFKISIHL', alt='TPFKISIH', molecule_type='protein')
[2024-03-25 11:40:59,906] {align.py:179} INFO - Found 2 alignments with max score, returning last in the list (3\' end rule)
[2024-03-25 11:40:59,907] {cli.py:45} INFO - Found best alignment with score = 5.0
[2024-03-25 11:40:59,909] {cli.py:46} INFO - Alignment:
ref               0 -PFKISIHL 8
                  0 -|||||||- 9
alt               0 TPFKISIH- 8

[2024-03-25 11:40:59,909] {cli.py:49} INFO - 2 Variant blocks generated
[2024-03-25 11:40:59,909] {cli.py:54} INFO - VariantBlock(alignment_block=Block(start=0, end=1, bases='i'), reference_blocks=[], alternate_blocks=[Block(start=0, end=1, bases='T')]), categorized as: extension
[2024-03-25 11:40:59,910] {cli.py:57} INFO - As HGVS: ref:p.Pro1extThr-1
[2024-03-25 11:40:59,910] {cli.py:54} INFO - VariantBlock(alignment_block=Block(start=8, end=9, bases='d'), reference_blocks=[Block(start=7, end=8, bases='L')], alternate_blocks=[]), categorized as: deletion
[2024-03-25 11:40:59,910] {cli.py:57} INFO - As HGVS: ref:p.Leu8del
```

## Usage - Python

Palamedes currently includes public functions in `__init__.py` that provide the alignment to HGVS functionality. All other functions should be treated as internal and private. No assurances are offered for their consistency and functionality from version to version. 

```python
>>> from Bio.Seq import Seq
>>> from Bio.SeqRecord import SeqRecord
>>> from palamedes import generate_hgvs_variants

# using raw strings
>>> generate_hgvs_variants("PFKISIHL", "TPFKISIH")
[
    SequenceVariant(ac=ref, type=p, posedit=Pro1extThr-1, gene=None),
    SequenceVariant(ac=ref, type=p, posedit=Leu8del, gene=None),
]

# using SeqRecord objects, note the molecule_type annotation must be provided and set to a supported type
# protein is the default and only supported option at this time
>>> ref = SeqRecord(Seq("PFKISIHL"), id="Jelleine-I", annotations={"molecule_type": "protein"})
>>> alt = SeqRecord(Seq("TPFKISIH"), id="Jelleine-IV", annotations={"molecule_type": "protein"})
>>> generate_hgvs_variants(ref, alt)
[
    SequenceVariant(ac=Jelleine-I, type=p, posedit=Pro1extThr-1, gene=None),
    SequenceVariant(ac=Jelleine-I, type=p, posedit=Leu8del, gene=None)
]
```

The `generate_hgvs_variants` also accepts a pre-built biopython `PairwiseAligner` instance to give the caller more control of the alignment parameters via the `aligner` keyword argument. The default settings are:
- mode: "global" (note this must be set to global on a custom aligner to ensure and end to end alignment)
- match_score: 1
- mismatch_score: -1
- open_gap_score: -1
- extend_gap_score: -0.1

## Name

The package is named after [Palamedes](https://en.wikipedia.org/wiki/Palamedes_(mythology)), a figure from Greek mythology. Palamedes was associated with the invention of the Greek letters and alphabet as well as with the invention of dice. Palamedes dedicated the first set of dice to the Greek goddess Tyche, who was the goddess of chance and randomness.
