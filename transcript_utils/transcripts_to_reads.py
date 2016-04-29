#!/usr/bin/env python

"""Usage:
    transcripts_to_reads create
        [--log-level=<log-level>]
        [--read-length=<read-length> --insert-size=<insert-size> --paired-end]
        <transcript-gtf-file> <transcript-fasta-file>
    transcripts_to_reads count
        [--log-level=<log-level>]
        [--read-length=<read-length> --insert-size=<insert-size> --paired-end]
        <transcript-gtf-file> <transcript-fasta-file>
    transcripts_to_reads count_mapped
        [--log-level=<log-level>]
        <gene-sequence-counts-file> <mapped-reads-bam-file>

-h --help                    Show this message.
-v --version                 Show version.
--log-level=<log-level>      Set logging level (one of {log_level_vals})
                             [default: info].
--read-length=<read-length>  Length of reads to create [default: 50].
--insert-size=<insert-size>  Length of fragments from which paired end reads
                             are created [default: 150].
--paired-end                 If specified, the reads produced will be paired-
                             rather than single-end.
<transcript-gtf-file>        File containing transcript definitions in GTF
                             format. This file is used to extract transcript
                             to gene mappings.
<transcript-fasta-file>      File containing transcript sequences in FASTA
                             format.
<gene-sequence-counts-file>  File containing per-gene counts of theoretical
                             reads, produced by the "count" command.
<mapped-reads-bam-file>      BAM file containing mapped reads.

Creates or counts theoretical RNA-seq reads created from a specified set of
transcripts. Can also count reads mapped in a BAM file.

Commands:
* create: For a specified set of input transcripts, create all per-gene
  theoretical reads of a certain length (and, in the case of paired-end reads,
  with a certain insert size), and output in FASTA format.
* count: For a specified set of input transcripts, count the number of
  theoretical reads of a certain length and insert size that could be created
  per-gene.
* count_mapped: For a specified set of input transcripts, count per-gene
  numbers of reads mapped in a BAM file.
"""

import collections
import docopt
import pysam
import schema
import sys

from . import gtf
from . import log
from . import options as opt
from . import subsequences as subs
from .__init__ import __version__

CREATE_READS_FASTA = "create"
COUNT_READS_FOR_GENES = "count"
COUNT_MAPPED_READS_FOR_GENES = "count_mapped"
LOG_LEVEL = "--log-level"
LOG_LEVEL_VALS = str(log.LOG_LEVELS.keys())
READ_LENGTH = "--read-length"
INSERT_SIZE = "--insert-size"
PAIRED_END = "--paired-end"
TRANSCRIPT_FASTA_FILE = "<transcript-fasta-file>"
TRANSCRIPT_GTF_FILE = "<transcript-gtf-file>"
GENE_SEQUENCE_COUNTS_FILE = "<gene-sequence-counts-file>"
MAPPED_READS = "<mapped-reads-bam-file>"


def validate_command_line_options(options):
    try:
        opt.validate_dict_option(
            options[LOG_LEVEL], log.LOG_LEVELS, "Invalid log level")
        opt.validate_file_option(
            options[TRANSCRIPT_GTF_FILE], "Transcript GTF file must exist")
        opt.validate_file_option(
            options[TRANSCRIPT_FASTA_FILE], "Transcript FASTA file must exist")
        opt.validate_file_option(
            options[GENE_SEQUENCE_COUNTS_FILE],
            "Gene sequence counts file must exist")
        opt.validate_file_option(
            options[MAPPED_READS], "Mapped reads BAM file must exist")
        options[READ_LENGTH] = opt.validate_int_option(
            options[READ_LENGTH], "Read length must be a positive integer", 1)
        options[INSERT_SIZE] = opt.validate_int_option(
            options[INSERT_SIZE],
            "Insert size must be greater than 2 * read length",
            min_val=2 * options[READ_LENGTH] + 1)
    except schema.SchemaError as exc:
        exit(exc.code)


def get_subsequence_generator(read_length, insert_size, paired_end):
    if paired_end:
        return subs.PairedEndSubsequences(read_length, insert_size)
    else:
        return subs.SingleEndSubsequences(read_length)


def get_sequences_for_genes(
        fasta_file, read_length, insert_size, paired_end,
        transcript_to_gene_mappings, logger):

    logger.info("Getting sequences of length {l} from {f}.".
                format(l=read_length, f=fasta_file))

    identifier = True
    current_gene = None
    gene_sequences = collections.defaultdict(set)

    subseq_generator = get_subsequence_generator(
        read_length, insert_size, paired_end)

    for line in open(fasta_file, 'r'):
        line = line.strip()
        if identifier:
            transcript = line[1:]
            if transcript not in transcript_to_gene_mappings:
                exit("{t} not found in transcript to gene mappings.".
                     format(t=transcript))
            current_gene = transcript_to_gene_mappings[transcript]
        else:
            subsequences = subseq_generator.get_subsequences(line)
            gene_sequences[current_gene] |= set(subsequences)

        identifier = not identifier

    return gene_sequences


def print_gene_sequence_fasta(gene_sequences, read_length,
                              insert_size, paired_end, logger):

    logger.info("Printing sequences for {n} genes.".
                format(n=len(gene_sequences)))

    subseq_generator = get_subsequence_generator(
        read_length, insert_size, paired_end)

    for gene, sequences in gene_sequences.items():
        if len(sequences) == 0:
            continue
        for i, sequence in enumerate(sequences):
            subseq_generator.print_subsequence(gene, i, sequence)


def print_gene_sequence_counts(gene_sequences, logger):
    logger.info("Printing sequence counts for {n} genes.".
                format(n=len(gene_sequences)))

    print("gene,seq_count")
    for gene, sequences in gene_sequences.items():
        if len(sequences) == 0:
            continue
        print("{g},{c}".format(g=gene, c=len(sequences)))


def read_sequence_counts_for_genes(counts_file, logger):
    logger.info("Reading sequence counts for genes from " + counts_file)
    sequence_counts = {}

    header = True
    for line in open(counts_file):
        if header:
            header = False
            continue
        gene, seq_count = line.strip().split(",")
        sequence_counts[gene] = int(seq_count)

    return sequence_counts


def get_mapped_reads(mapped_reads_file, logger):
    logger.info("Reading mapped sequences from " + mapped_reads_file)
    mapped_sequences = collections.defaultdict(set)

    samfile = pysam.Samfile(mapped_reads_file, "rb")
    for read in samfile.fetch(until_eof=True):
        read_elems = read.qname.split(":")
        gene = read_elems[0]
        read_id = gene + ":" + read_elems[1]
        mapped_sequences[gene].add(read_id)

    return mapped_sequences


def print_mapped_sequence_counts(sequence_counts, mapped_sequences, logger):
    logger.info("Printing mapped sequence counts for {n} genes".
                format(n=len(sequence_counts)))
    print("gene,seq_count,mapped_seq_count,mapped_seq_frac")
    for gene, count in sequence_counts.items():
        mapped_count = len(mapped_sequences[gene]) \
            if gene in mapped_sequences else 0
        mapped_frac = float(mapped_count) / count
        print("{g},{c},{mc},{mf}".format(
            g=gene, c=count, mc=mapped_count, mf=mapped_frac))


def create_reads_from_transcripts(
        gtf_file, fasta_file, read_length, insert_size,
        paired_end, logger):

    logger.info("Creating reads from transcripts...")

    transcript_to_gene_mappings = \
        gtf.GtfInfo(gtf_file, logger).get_transcript_to_gene_mappings()
    gene_sequences = get_sequences_for_genes(
        fasta_file, read_length, insert_size, paired_end,
        transcript_to_gene_mappings, logger)
    print_gene_sequence_fasta(
        gene_sequences, read_length,
        insert_size, paired_end, logger)


def count_reads_for_genes(
        gtf_file, fasta_file, read_length, insert_size,
        paired_end, logger):

    logger.info("Counting reads for genes...")

    transcript_to_gene_mappings = \
        gtf.GtfInfo(gtf_file, logger).get_transcript_to_gene_mappings()
    gene_sequences = get_sequences_for_genes(
        fasta_file, read_length, insert_size, paired_end,
        transcript_to_gene_mappings, logger)
    print_gene_sequence_counts(gene_sequences, logger)


def count_mapped_reads_for_genes(counts_file, mapped_reads, logger):
    logger.info("Counting mapped reads for genes...")
    sequence_counts = read_sequence_counts_for_genes(counts_file, logger)
    mapped_reads = get_mapped_reads(mapped_reads, logger)
    print_mapped_sequence_counts(sequence_counts, mapped_reads, logger)


def transcripts_to_reads(args):
    docstring = __doc__.format(log_level_vals=LOG_LEVEL_VALS)
    options = docopt.docopt(docstring, argv=args,
                            version="transcripts_to_reads v" + __version__)
    validate_command_line_options(options)

    logger = log.get_logger(sys.stderr, options[LOG_LEVEL])

    gtf_file = options[TRANSCRIPT_GTF_FILE]
    fasta_file = options[TRANSCRIPT_FASTA_FILE]
    read_length = options[READ_LENGTH]
    insert_size = options[INSERT_SIZE]
    paired_end = options[PAIRED_END]

    if options[CREATE_READS_FASTA]:
        create_reads_from_transcripts(
            gtf_file, fasta_file, read_length, insert_size,
            paired_end, logger)
    elif options[COUNT_READS_FOR_GENES]:
        count_reads_for_genes(
            gtf_file, fasta_file, read_length, insert_size,
            paired_end, logger)
    elif options[COUNT_MAPPED_READS_FOR_GENES]:
        count_mapped_reads_for_genes(
            options[GENE_SEQUENCE_COUNTS_FILE],
            options[MAPPED_READS], logger)
