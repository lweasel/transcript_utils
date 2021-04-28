from setuptools import setup

import transcript_utils

setup(
    name='transcript_utils',
    version=transcript_utils.__version__,
    url='https://github.com/lweasel/transcript_utils',
    license='MIT License',
    author='Owen Dando',
    author_email='owen.dando@ed.ac.uk',
    packages=['transcript_utils'],
    install_requires=[
        'docopt',
        'numpy',
        'pandas',
        'pysam',
        'python-dateutil',
        'pytz',
        'schema',
        'six',
    ],
    scripts=[
        'bin/create_unspliced_transcripts_gtf',
        'bin/get_gene_lengths',
        'bin/transcripts_to_reads',
    ]
)
