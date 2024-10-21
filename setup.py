from setuptools import find_packages, setup

version = "0.0.9"

setup(
    name="palamedes",
    version=version,
    url="https://github.com/mammothbio-os/palamedes",
    python_requires=">=3.11",
    description="Palamedes: HGVS variants from a sequence alignment",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author_email="open-source@mammothbiosci.com",
    entry_points={
        "console_scripts": ["palamedes=palamedes.__main__:main"],
    },
    packages=find_packages(),
    install_requires=["biopython~=1.83", "hgvs~=1.5.4"],
    keywords=[
        "bioinformatics",
        "alignment",
        "hgvs",
        "variant",
        "amino acid",
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.11",
  ],
)
