from setuptools import find_packages, setup

version = '0.0.0'

setup(
    name='palamedes',
    version=version,
    description='Palamedes: HGVS variants from a sequence alignment',
    author_email='software-prod@mammothbiosci.com',
    entry_points={
        'console_scripts': ['palamedes=palamedes.cli:main'],
    },
    packages=find_packages(),
    install_requires=['biopython~=1.83', 'hgvs~=1.5.4'],
)
