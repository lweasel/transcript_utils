BASE_COMPLEMENTS = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


class SingleEndSubsequences(object):
    def __init__(self, read_length):
        self.read_length = read_length

    def get_subsequences(self, sequence):
        subsequences = []
        if len(sequence) == self.read_length:
            subsequences = [sequence]
        elif len(sequence) > self.read_length:
            subsequences = \
                [sequence[i:i + self.read_length]
                 for i in range(len(sequence) - self.read_length + 1)]

        return subsequences

    def print_subsequence(self, gene, subseq_no, subsequence):
        print(">{gene}:{i}".format(gene=gene, i=subseq_no))
        print subsequence


class PairedEndSubsequences(object):
    def __init__(self, read_length, insert_size):
        self.read_length = read_length
        self.insert_size = insert_size

    def get_subsequences(self, sequence):
        insert_sequences = []
        if len(sequence) == self.insert_size:
            insert_sequences = [sequence]
        elif len(sequence) > self.insert_size:
            insert_sequences = \
                [sequence[i:i + self.insert_size]
                 for i in range(len(sequence) - self.insert_size + 1)]

        subsequences = [(ins[:self.read_length] + ins[-self.read_length:])
                        for ins in insert_sequences]

        return subsequences

    def print_subsequence(self, gene, subseq_no, subsequence):
        print(">{gene}:{i}:_1".format(gene=gene, i=subseq_no))
        print(subsequence[:self.read_length])
        print(">{gene}:{i}:_2".format(gene=gene, i=subseq_no))
        print(self.reverse_complement(subsequence[-self.read_length:]))

    def reverse_complement(self, sequence):
        return "".join(
            BASE_COMPLEMENTS.get(base, base) for base in reversed(sequence))
