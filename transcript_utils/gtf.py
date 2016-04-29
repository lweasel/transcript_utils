import pandas

GTF_FEATURE_COL = 2
GTF_ATTRIBUTES_COL = 8
GTF_GENE_ID_ATTRIBUTE = "gene_id"
GTF_TRANSCRIPT_ID_ATTRIBUTE = "transcript_id"


def get_gtf_attributes_dict(attributes_str):
    strip_quotes = lambda x: x.replace('"', '')
    return {attr: strip_quotes(val) for attr, val in
            [av.split(" ", 1) for av in attributes_str.split("; ")]}


def get_transcript_to_gene_mappings(gtf_file, logger):
    logger.info("Getting transcript->gene mappings from " + gtf_file)
    gtf_info = pandas.read_csv(gtf_file, sep='\t', header=None)
    transcript_to_gene_mappings = {}
    lines_processed = 0

    for index, row in gtf_info.iterrows():
        lines_processed += 1
        if lines_processed % 10000 == 0:
            logger.debug("Processed {l} GTF lines.".format(l=lines_processed))

        if row[GTF_FEATURE_COL] != "transcript":
            continue

        attributes_dict = get_gtf_attributes_dict(row[GTF_ATTRIBUTES_COL])
        transcript = attributes_dict[GTF_TRANSCRIPT_ID_ATTRIBUTE]
        if transcript not in transcript_to_gene_mappings:
            transcript_to_gene_mappings[transcript] = \
                attributes_dict[GTF_GENE_ID_ATTRIBUTE]

    return transcript_to_gene_mappings
