#!/usr/bin/env python

"""Usage:
    get_gene_lengths [--log-level=<log-level>] <transcript-gtf-file>

-h --help                    Show this message.
-v --version                 Show version.
--log-level=<log-level>      Set logging level (one of {log_level_vals})
                             [default: info].
<transcript-gtf-file>        File containing transcript definitions in GTF
                             format.

Write per-gene maximum transcript lengths and gene lengths (the latter defined
as the number of bases contained in the union of all exons of the gene,
including 3' and 5' UTRs).
"""

import docopt
import schema
import sys

from collections import defaultdict

from . import gtf
from . import log
from . import options as opt
from .__init__ import __version__

LOG_LEVEL = "--log-level"
LOG_LEVEL_VALS = str(log.LOG_LEVELS.keys())

TRANSCRIPT_GTF_FILE = "<transcript-gtf-file>"


def validate_command_line_options(options):
    try:
        opt.validate_dict_option(
            options[LOG_LEVEL], log.LOG_LEVELS, "Invalid log level")
        opt.validate_file_option(
            options[TRANSCRIPT_GTF_FILE], "Transcript GTF file must exist")
    except schema.SchemaError as exc:
        exit(exc.code)


def _calculate_gene_lengths(transcript_info, logger):
    logger.info("Calculating gene lengths...")
    gene_lengths = {}

    for gene_name, gene in transcript_info.iteritems():
        gene_start = sys.maxint
        gene_end = -sys.maxint - 1
        exon_starts = defaultdict(list)
        exon_ends = defaultdict(list)

        for transcript in gene.transcripts.values():
            for exon in transcript.exons:
                if exon.start < gene_start:
                    gene_start = exon.start
                if exon.end > gene_end:
                    gene_end = exon.end

                exon_starts[exon.start].append(exon)
                exon_ends[exon.end].append(exon)

        locus = gene_start - 1
        depth = 0
        gene_length = 0

        while locus < gene_end:
            locus += 1

            num_starts = 0
            if locus in exon_starts:
                num_starts = len(exon_starts[locus])

            depth += num_starts

            if depth > 0:
                gene_length += 1

            num_ends = 0
            if locus in exon_ends:
                num_ends = len(exon_ends[locus])

            depth -= num_ends

        gene_lengths[gene_name] = gene_length

        if len(gene_lengths) % 1000 == 0:
            logger.info("...processed {g} genes.".format(g=len(gene_lengths)))

    return gene_lengths


def _calculate_max_transcript_lengths(transcript_info, logger):
    logger.info("Calculating maximum transcript lengths...")
    max_transcript_lengths = {}

    for gene_name, gene in transcript_info.iteritems():
        max_transcript_length = 0

        for transcript in gene.transcripts.values():
            transcript_length = 0

            for exon in transcript.exons:
                transcript_length += exon.length()

            logger.debug("Transcript {t} length {l}".format(
                t=transcript.name, l=transcript_length))

            if transcript_length > max_transcript_length:
                max_transcript_length = transcript_length

        max_transcript_lengths[gene_name] = max_transcript_length

        if len(max_transcript_lengths) % 1000 == 0:
            logger.info("...processed {g} genes.".format(
                g=len(max_transcript_lengths)))

    return max_transcript_lengths


def _print_lengths(gene_lengths, max_transcript_lengths, logger):
    logger.info("Printing gene and maximum transcript lengths...")

    print("gene,gene_length,max_transcript_length")

    for gene in sorted(gene_lengths.keys()):
        gene_length = gene_lengths[gene]
        max_transcript_length = max_transcript_lengths[gene]

        if max_transcript_length > gene_length:
            exit(("{g}: Max transcript length ({mtl}) cannot exceed " +
                  "gene length ({gl})").format(
                g=gene, mtl=max_transcript_length, gl=gene_length))

        print("{g},{gl},{mtl}".format(
            g=gene, mtl=max_transcript_length, gl=gene_length))


def get_gene_lengths(args):
    docstring = __doc__.format(log_level_vals=LOG_LEVEL_VALS)
    options = docopt.docopt(docstring,
                            version="get_gene_lengths v" + __version__)
    validate_command_line_options(options)

    logger = log.get_logger(sys.stderr, options[LOG_LEVEL])

    gtf_info = gtf.GtfInfo(options[TRANSCRIPT_GTF_FILE], logger)
    transcript_info = gtf_info.get_transcript_info()
    gene_lengths = _calculate_gene_lengths(transcript_info, logger)
    max_transcript_lengths = _calculate_max_transcript_lengths(
        transcript_info, logger)

    _print_lengths(gene_lengths, max_transcript_lengths, logger)
