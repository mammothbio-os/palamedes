from setuptools import find_packages, setup

version = '0.0.0'

setup(
    name='palamedes',
    version=version,
    python_requires='>=3.11',
    description='Palamedes: HGVS variants from a sequence alignment',
    author_email='open-source@mammothbiosci.com',
    entry_points={
        'console_scripts': ['palamedes=palamedes.__main__:main'],
    },
    packages=find_packages(),
    install_requires=['biopython~=1.83', 'hgvs~=1.5.4'],
)
