"""-------------------------------------------------------------------------------------------------
Excursion analysis can identify regions that have a locally high density region of specific features

24 September 2018   Michael Gribskov

-------------------------------------------------------------------------------------------------"""
import sys
import re
from sequence.fasta import Fasta


class Feature:
    """---------------------------------------------------------------------------------------------
    Features to search for
    ---------------------------------------------------------------------------------------------"""

    def __init__(self):
        self.feature = '.'
        self.window = 0
        self.space = {}

    def count_space(self, fasta):
        """-----------------------------------------------------------------------------------------
        count the spacing between features.  There must be at least two occurences
        :param fasta:
        :return:
        -----------------------------------------------------------------------------------------"""
        result = re.finditer(r'(?=({}))'.format(self.feature), fasta.seq)

        n = 0
        pos = None
        for hit in result:
            if pos:
                n += 1
                ln = hit.start(1) - pos
                print('ln', ln)
                try:
                    self.space['{}'.format(ln)] += 1
                except KeyError:
                    self.space['{}'.format(ln)] = 1

            else:
                pos = hit.start(1)

        return n

    def total(self):
        """-----------------------------------------------------------------------------------------
        return total counts and total length
        usafe
            counts, bases = feature.total()
        -----------------------------------------------------------------------------------------"""
        total_count = 0
        total_len = 0
        for lspace in self.space:
            total_count += self.space[lspace]
            total_len += self.space[lspace] * int(lspace)

        return total_count, total_len


# ==================================================================================================
# main/test
# ==================================================================================================
if __name__ == '__main__':
    fasta = Fasta(file=sys.argv[1])
    feature = Feature()
    feature.feature = "A.C"

    c = 0;
    while fasta.next():
        count = feature.count_space(fasta)
        print(count)
        c += 1
        if c > 100:
            break

    fcount, flen = feature.total()
    print('count:{}  len:{}  avg:{}'.format( fcount, flen, flen/fcount))

exit(0)
