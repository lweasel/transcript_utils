import pandas


class GtfRow(object):
    FEATURE_COL = 2
    START_COL = 3
    END_COL = 4
    ATTRIBUTES_COL = 8

    EXON_FEATURE = "exon"

    GENE_ID_ATTRIBUTE = "gene_id"
    TRANSCRIPT_ID_ATTRIBUTE = "transcript_id"

    def __init__(self, row_data):
        self.row_data = row_data

        strip_quotes = lambda x: x.replace('"', '')
        attr_str = self.row_data[GtfRow.ATTRIBUTES_COL]
        self.attr_dict = {attr: strip_quotes(val) for attr, val in
                          [av.split(" ", 1) for av in attr_str.split("; ")]}

    def get_feature(self):
        return self.row_data[GtfRow.FEATURE_COL]

    def get_start(self):
        return self.row_data[GtfRow.START_COL]

    def get_end(self):
        return self.row_data[GtfRow.END_COL]

    def get_gene(self):
        return self.attr_dict[GtfRow.GENE_ID_ATTRIBUTE]

    def get_transcript(self):
        return self.attr_dict[GtfRow.TRANSCRIPT_ID_ATTRIBUTE]

    def is_exon(self):
        return self.get_feature() == GtfRow.EXON_FEATURE


class GtfInfo(object):
    def __init__(self, gtf_file, logger):
        self.gtf_file = gtf_file
        self.data = pandas.read_csv(
            gtf_file, sep="\t", header=None, comment="#")
        self.logger = logger

    def rows(self):
        for index, row in self.data.iterrows():
            yield GtfRow(row)

    def get_transcript_to_gene_mappings(self):
        self.logger.info(
            "Getting transcript->gene mappings from " + self.gtf_file)

        transcript_to_gene_mappings = {}
        lines_processed = 0

        for row in self.rows():
            lines_processed += 1
            if lines_processed % 10000 == 0:
                self.logger.debug("Processed {l} GTF lines.".format(l=lines_processed))

            if row.get_feature() != "transcript":
                continue

            transcript = row.get_transcript()
            if transcript not in transcript_to_gene_mappings:
                transcript_to_gene_mappings[transcript] = row.get_gene()

        return transcript_to_gene_mappings
