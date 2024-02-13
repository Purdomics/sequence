"""=================================================================================================
Significance of alignment using random simulation

Michael Gribskov     02 August 2020
================================================================================================="""
import sys
from math import log10
import random
from scipy import stats
from sequence.fasta import Fasta
from sequence.score import Score


class Alignment(Score):
    """---------------------------------------------------------------------------------------------

    ---------------------------------------------------------------------------------------------"""

    def __init__(self):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        Score.__init__(self)
        self.s1 = None
        self.s2 = None
        self.i1 = None
        self.i2 = None

    def seqToInt(self):
        """-----------------------------------------------------------------------------------------
        Convert sequence strings to an integer arrays and stores in object.  An integer array is
        more convenient for direct lookups in the scoring table than a string

        :return: int, int length of sequence lists
        -----------------------------------------------------------------------------------------"""
        a2i = self.a2i

        self.i1 = [a2i[c] for c in self.s1.seq]
        self.i2 = [a2i[c] for c in self.s2.seq]

        return len(self.i1), len(self.i2)

    def localScore(self, open, extend):
        """-----------------------------------------------------------------------------------------
        Local alignment score only.
        s1 is the horizontal sequence and s2 is the vertical sequence.  this makes s2 the row
        index and s1 the column index.

        :param open: float, gap opening penalty
        :param extend: float, gap extension penalty
        :return:
        -----------------------------------------------------------------------------------------"""
        cmp = self.table
        i1 = self.i1
        i2 = self.i2
        l1 = len(i1)
        l2 = len(i2)

        edge = self.min + open  # a small value to use for gaps on the edges

        scoremax = 0
        posmax = [0, 0]
        bestrow = edge
        bestcol = [edge for _ in i1]
        score = [0 for _ in i1]
        cell = 0
        diag = 0

        jpos = 0
        for j in i2:
            # print(score)
            diag = 0  # uncomment to print score matrix
            bestrow = edge + extend
            ipos = 0
            for i in i1:
                cell = max(cmp[j][i] + max(diag, bestcol[ipos], bestrow), 0.0)
                if cell > scoremax:
                    scoremax = cell
                    posmax = [ipos, jpos]

                bestrow = max(diag + open, bestrow + extend)
                bestcol[ipos] = max(diag + open, bestcol[ipos] + extend)
                diag = score[ipos]
                score[ipos] = cell
                ipos += 1

            jpos += 1

        # print(score)
        return scoremax, posmax

    def localReverse(self, open, extend, bestscore, bestpos):
        """-----------------------------------------------------------------------------------------
        Start from the endpoint of the best alignment and do an alignment in reverse to find the
        beginning point
        :param open: float, gap opening penalty
        :param extend: float, gap extension penalty
        :param bestscore: float, score of best alignment
        :param bestpos: list of 2 float, position of best alignment
        :return: float (alignment score), list of two floats (position of max)
        -----------------------------------------------------------------------------------------"""
        i, j = bestpos
        i1save = self.i1
        i2save = self.i2

        self.i1 = self.i1[0:i:-1]
        self.i2 = self.i2[0:j:-1]

        score, pos = self.localScore(open, extend)
        print('begin pos: {} = [{},{}]'.format(pos, bestpos[0] - pos[0], bestpos[1] - pos[1]))

        self.i1 = i1save
        self.i2 = i2save

        return score, [bestpos[0] - pos[0], bestpos[1] - pos[1]]


        # def traceback(self, open, extend, bestscore, bestpos):
        #     """-----------------------------------------------------------------------------------------
        #     traceback without score matrix
        #
        #     :param open: float, gap opening penalty
        #     :param extend: float, gap extension penalty
        #     :param bestscore: float, score of best alignment
        #     :param bestpos: list of 2 float, position of best alignment
        #     :return: string, string, the padded aligned sequences
        #     -----------------------------------------------------------------------------------------"""
        #     cmp = self.table
        #     i1 = self.i1
        #     i2 = self.i2
        #     l1 = len(i1)
        #     l2 = len(i2)
        #
        #     score = bestscore
        #     ipos, jpos = bestpos
        #     while score > 0:
        target


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    align = Alignment()

    align.s1 = Fasta(filename=sys.argv[1])
    align.s2 = Fasta(filename=sys.argv[2])
    align.readNCBI('../../dotplot/table/BLOSUM62.matrix')

    # testing
    # align.s1 = Fasta()
    # align.s1.seq = 'ACTGCC'
    # align.s2 = Fasta()
    # align.s2.seq = 'ATGCC'
    # align.readNCBI('../../dotplot/table/NUC4.4.matrix')

    align.seqToInt()
    # random.shuffle(align.i1)          # uncomment to test scores for random alignments
    original_score, bestpos = align.localScore(-10, -1)
    print('original score: {} at {}'.format(original_score, bestpos))
    beginscore, beginpos = align.localReverse(-10,-1, original_score, bestpos)
    exit(1)

    nrandom = 10
    rscore = []
    for i in range(nrandom):
        random.shuffle(align.i2)
        score, pos = align.localScore(-20, -2)
        rscore.append(score)
        # this just gives something to look at
        print('\t{}: {}'.format(i, score))

    # tabulate observed number >= score and peak of score distribution
    i = 0
    N = []
    last = 0
    peak = 0
    peakval = 0
    rscore.sort(reverse=True)
    ro = rscore[0]
    for r in rscore:
        if r < ro:
            N.append([i, ro])
            interval = i - last
            if interval >= peakval:
                peak = ro
                peakval = interval
            last = i
            ro = r

        i += 1
        print('\t{}\t{}\t{:.3f}'.format(i, log10(i), r))
    N.append([i, ro])

    print('peak:{}  count:{}'.format(peak, peakval))

    # fit a straight line from point 2 to the peak
    x = []
    y = []
    i = 0
    for pair in N:
        print('\t{}\t{:.3f}\t{}'.format(pair[0], log10(pair[0]), pair[1]))
        if i < 1 and i < peak:
            y.append(log10(pair[0]))
            x.append(pair[1])

    slope, intercept, r, p, std_err = stats.linregress(x, y)
    print('slope:{:4g}\tint:{:4f}\tr:{:4f}\tp:{:4f}\tstd err:{:4f}'.format(slope, intercept, r, p,
                                                                           std_err))

    logE = slope * original_score + intercept
    E = 10 ** (logE)
    print('E({}) = {:.2g}\tlog(E)={}'.format(original_score, E, logE))
    exit(0)
