#!/usr/bin/env python

"""Usage:
    create_unspliced_transcripts_gtf [--log-level=<log-level>] <transcript-gtf-file>

-h --help                    Show this message.
-v --version                 Show version.
--log-level=<log-level>      Set logging level (one of {log_level_vals})
                             [default: info].
<transcript-gtf-file>        File containing transcript definitions in GTF
                             format.

Write a GTF file containing transcripts corresponding to those in the input GTF
file, but with all introns removed.
"""

import docopt
import schema
import sys

from . import gtf
from . import log
from . import options as opt
from .__init__ import __version__

LOG_LEVEL = "--log-level"
LOG_LEVEL_VALS = str(log.LOG_LEVELS.keys())

TRANSCRIPT_GTF_FILE = "<transcript-gtf-file>"


def _validate_command_line_options(options):
    try:
        opt.validate_dict_option(
            options[LOG_LEVEL], log.LOG_LEVELS, "Invalid log level")
        opt.validate_file_option(
            options[TRANSCRIPT_GTF_FILE], "Transcript GTF file must exist")
    except schema.SchemaError as exc:
        exit(exc.code)


def _write_unspliced_gtf(transcript_info):
    gtf_rows = []

    for gene in transcript_info.values():
        for transcript in gene.transcripts.values():
            transcript_start = sys.maxsize
            transcript_end = -sys.maxsize - 1

            for exon in transcript.exons:
                if exon.start < transcript_start:
                    transcript_start = exon.start
                if exon.end > transcript_end:
                    transcript_end = exon.end

            gtf_rows.append(gtf.GtfRow.from_values(
                gene.seqname, "exon", transcript_start, transcript_end,
                transcript.strand, gene.name, transcript.name))

    for row in gtf_rows:
        print(row)


def create_unspliced_transcripts_gtf(args):
    docstring = __doc__.format(log_level_vals=LOG_LEVEL_VALS)
    options = docopt.docopt(
        docstring, version="create_unspliced_transcripts_gtf v" + __version__)
    _validate_command_line_options(options)

    logger = log.get_logger(sys.stderr, options[LOG_LEVEL])

    gtf_info = gtf.GtfInfo(options[TRANSCRIPT_GTF_FILE], logger)
    transcript_info = gtf_info.get_transcript_info()
    _write_unspliced_gtf(transcript_info)
