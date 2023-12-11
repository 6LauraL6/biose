from Bio.Seq import Seq


def info(seq):
    seq = Seq(seq)

    result = str(seq.reverse_complement())

    return result

