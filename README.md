transcript-utils
================

A collection of simple scripts to perform various tasks for gene transcripts and their sequences:

* **create_unspliced_transcripts_gtf**: Writes a GTF file containing transcripts corresponding to those in the input GTF file, but with all introns removed (i.e. each transcript contains a single "exon"; only GTF lines for exon features are written).
* **get_gene_lengths**: Calculate a couple of measures of gene "size". Given a GTF file as input, writes in CSV format the gene ID, maximum transcript length among isoforms of the gene, and gene "length" (the latter being defined as the number of bases contained in the union of all exons of all isoforms of the gene, including 3' and 5' UTRs).
* **transcripts_to_reads**: Create all possible unique theoretical RNA-seq reads (i.e. with no errors) from a set of transcript sequences. Reads can be single- or paired-end (with a fixed insert size in the latter case). The script can also count, per-gene, the number of such reads mapped in a BAM file.

Installation
============

As these scripts have a number of dependencies on other Python packages, it is recommended to install this package in an isolated environment using the [virtualenv](http://virtualenv.readthedocs.org/en/latest/index.html>) tool. The [virtualenvwrapper](http://virtualenvwrapper.readthedocs.org/en/latest/install.html>) tool makes managing multiple virtual environments easier.

After setting up ``virtualenv`` and ``virtualenvwrapper``, create and work in a virtual environment using the ``virtualenvwrapper`` tool:

```
mkproject transcript-utils
```

Install the *transcript-utils* package and its dependencies, into the virtual environment by running:

```
pip install git+https://github.com/sidbdri/transcript-utils.git
```

in the top-level project directory.

create_unspliced_transcripts_gtf
================================

Usage:

```
    create_unspliced_transcripts_gtf [--log-level=<log-level>] <transcript-gtf-file>
```

Single-exon transcripts are written in GTF format to standard out.


get_gene_lengths
================

Usage:

```
    get_gene_lengths [--log-level=<log-level>] <transcript-gtf-file>
```

Gene "size" information is written in CSV format to standard out.

transcripts_to_reads
====================

This script has three commands:

Create reads
------------

Usage:

```
    transcripts_to_reads create
        [--log-level=<log-level>]
        [--read-length=<read-length> --insert-size=<insert-size> --paired-end]
        <transcript-gtf-file> <transcript-fasta-file>
```

Writes, per-gene, all possible unique theoretical RNA-seq reads of a given length (default: 50) that could be created from transcripts of that gene. Reads are single-end by default; specifying ``--paired-end`` causes paired-end reads to be written, in which case ``--insert-size`` (default: 150) determines the (fixed) size of fragments from the ends of which reads will be created.

Reads are written to standard out in FASTA format. In the case of paired-end reads, left and right reads for each pair are interleaved.

The following must be supplied:

* A FASTA file containing transcript sequences
* A GTF file from which transcript -> gene mappings will be extracted.

Count reads
-----------

Usage:

```
    transcripts_to_reads count
        [--log-level=<log-level>]
        [--read-length=<read-length> --insert-size=<insert-size> --paired-end]
        <transcript-gtf-file> <transcript-fasta-file>
```

As above, but instead simply counts, per-gene, the number of unique reads that would be created. Counts are written to standard out in CSV format::

    <gene>,<unique_read_count>

Again, a FASTA transcript sequence file and a GTF annotation file must be supplied.

Count mapped reads
------------------

Usage:

```
    transcripts_to_reads count_mapped
        [--log-level=<log-level>]
        <gene-sequence-counts-file> <mapped-reads-bam-file>
```

Counts, per-gene, the number of unique theoretical RNA-seq reads that have been mapped in a BAM file, and also calculates the fraction of all possible such reads that appear in the BAM file. Counts are written to standard out in CSV format::

    <gene>,<unique_read_count>,<mapped_read_count>,<mapped_read_fraction>

The following must be supplied:

* A CSV file containing per-gene counts of unique theoretical RNA-seq reads, as output by the ``count`` command above.
* A BAM file containing mapped theoretical reads.
