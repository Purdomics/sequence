"""=================================================================================================
just a trial of using an object to represent the cell and pointers.  since it is a local alignment
only non'zero cells need to be represented.  maybe this will make it sparse enough to be worthwhile

Michael Gribskov     04 August 2020
================================================================================================="""
import sys
from math import log10
import random
from scipy import stats
from sequence.fasta import Fasta
from sequence.score import Score


class Cell:
    """=============================================================================================

    ============================================================================================="""
    count = 0

    def __init__(self):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        self.n = Cell.count
        Cell.count += 1
        self.score = 0
        self.p = []


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
        self.score = None

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

    def globalBrute(self, open, extend, nogap=False):
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

        edge = Cell()  # a dummy cell for edges
        edge.score = self.min + open + l1 * extend

        score = [[Cell() for i in i1] for j in i2]
        self.score = score
        bestrow = Cell()
        bestcol = [Cell() for i in i1]

        score[0][0].score = cmp[i2[0]][i1[0]]
        gap = open
        for ipos in range(l1):
            bestcol[ipos].score = edge.score
            bestcol[ipos].p = []
            gap += extend

        diag = Cell()

        jpos = 0
        vgap = 0
        diag.score = 0
        for j in i2:
            bestrow.p = []
            bestrow.score = edge.score

            ipos = 0
            for i in i1:
                previous = max(diag.score, bestcol[ipos].score, bestrow.score)
                cell = cmp[j][i] + previous

                for dir in (diag, bestrow, bestcol[ipos]):
                    if dir.score == previous:
                        # set pointers for all directions
                        score[jpos][ipos].p.append(dir)

                # update best row and column values
                if diag.score + open > bestrow.score + extend:
                    # what if scores are equal? some paths missed
                    bestrow.p = diag.p
                    bestrow.score = diag.score + open
                else:
                    bestrow.score += extend

                if diag.score + open > bestcol[ipos].score + extend:
                    bestcol[ipos].p = diag.p
                    bestcol[ipos].score = diag.score + open
                else:
                    bestcol[ipos].score += extend

                # diagonal score for next cell
                if jpos > 0:
                    diag.p = score[jpos - 1][ipos]
                    diag.score = score[jpos - 1][ipos].score
                else:
                    diag.score = edge.score

                score[jpos][ipos].score = cell
                ipos += 1

            # end of loop over columns

            # special case for first cell in each row, best previous is always a column gap
            if jpos:
                vgap += extend
            else:
                vgap = open
            bestcol[0].score = vgap

            diag.score = edge.score
            diag.p = []

            jpos += 1

            # end of loop over rows

        # add the end gap penalties
        if not nogap:
            gap = open
            jpos = len(i2) - 1
            for ipos in range(len(i1) - 2, -1, -1):
                score[jpos][ipos].score += gap
                gap += extend

            gap = open
            ipos = len(i1) - 1
            for jpos in range(len(i2) - 2, -1, -1):
                score[jpos][ipos].score += gap
                gap += extend

        scoremax = 1
        posmax = [1, 1]
        return scoremax, posmax

    def localBrute(self, open, extend):
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

        edge = Cell()  # a dummy cell for the edges
        edge.score = self.min + open

        scoremax = 0
        posmax = [0, 0]
        score = [[Cell() for i in i1] for j in i2]
        self.score = score
        bestrow = Cell()
        # bestrow.p = edge
        # bestrow.score = edge.score
        bestcol = []
        for _ in i1:
            c = Cell()
            c.p = [edge]
            c.score = edge.score
            bestcol.append(c)

        for ipos in range(l1):
            score[0][ipos].score = edge.score
            # score[0][ipos].p = [edge]

        diag = Cell()

        jpos = 0
        for j in i2:
            bestrow.p = [edge]
            bestrow.score = edge.score + extend
            diag.p = [edge]
            diag.score = edge.score

            ipos = 0
            for i in i1:
                previous = max(diag.score, bestcol[ipos].score, bestrow.score, 0.0)
                cell = max(cmp[j][i] + previous, 0.0)
                if cell > 0:
                    # score[jpos][ipos].score = cell
                    if cell > scoremax:
                        scoremax = cell
                        posmax = [jpos, ipos]

                    if previous > 0:
                        # only set pointers if cell score > zero and there is a non-zero best
                        # previous score
                        for dir in (diag, bestrow, bestcol[ipos]):
                            if dir.score == previous:
                                # set pointers for all directions
                                score[jpos][ipos].p.append(dir.p)

                # update best row and column values
                if diag.score + open > bestrow.score + extend:
                    # what if scores are equal? some paths missed
                    bestrow.p = diag.p
                    bestrow.score = diag.score + open
                else:
                    bestrow.score += extend

                if diag.score + open > bestcol[ipos].score + extend:
                    bestcol[ipos].p = diag.p
                    bestcol[ipos].score = diag.score + open
                else:
                    bestcol[ipos].score += extend

                # diagonal score for next cell
                if jpos > 0:
                    diag.p = score[jpos - 1][ipos]
                    diag.score = score[jpos - 1][ipos].score

                score[jpos][ipos].score = cell
                ipos += 1

            jpos += 1

        return scoremax, posmax

    def trace1(self, pos):
        """-----------------------------------------------------------------------------------------
        Trace back one alignment using the first pointer for each cell

        :param pos: list of 2 int, traceback start position
        :return:
        -----------------------------------------------------------------------------------------"""
        s1 = self.s1.seq
        s2 = self.s2.seq
        l1 = len(s1)
        l2 = len(s2)
        cmp = self.table
        a2i = self.a2i

        score = self.score
        current = score[pos[0]][pos[1]]

        a1 = ''
        a2 = ''
        rowold = pos[0]
        colold = pos[1]
        while len(current.p) > 0:
            row, col = Alignment.n2pos(l1, current.n)

            for c in range(colold - 1, col, -1):
                a1 += s1[c]
                a2 += '.'

            for r in range(rowold - 1, row, -1):
                a1 += '.'
                a2 += s2[r]

            a1 += s1[col]
            a2 += s2[row]

            if len(current.p):
                current = current.p[0]

            rowold = row
            colold = col

        m = self.matchString(a1, a2)

        return a1[::-1], a2[::-1], m[::-1]

    def traceAll(self, pos):
        """-----------------------------------------------------------------------------------------
        Trace back one alignment using the first pointer for each cell

        :param pos: list of 2 int, traceback start position
        :return:
        -----------------------------------------------------------------------------------------"""
        s1 = self.s1.seq
        s2 = self.s2.seq
        l1 = len(s1)
        l2 = len(s2)
        cmp = self.table
        a2i = self.a2i

        score = self.score
        stack = []
        a1 = ' ' * (l1 * l2)
        a2 = ' ' * (l1 * l2)
        n = score[pos[0]][pos[1]].n
        nold = score[pos[0]][pos[1]].n
        alen = 0
        stack.append([n, nold, alen])

        save = []
        while stack:
            n, nold, alen = stack.pop()
            row, col = Alignment.n2pos(l1, n)
            rowold, colold = Alignment.n2pos(l1, nold)
            a1 = a1[:alen]
            a2 = a2[:alen]

            for c in range(colold - 1, col, -1):
                a1 += s1[c]
                a2 += '.'
                alen += 1

            for r in range(rowold - 1, row, -1):
                a1 += '.'
                a2 += s2[r]
                alen += 1

            a1 += s1[col]
            a2 += s2[row]
            alen += 1

            # for each path in the pointers of the current cell push on stack
            ptrs = score[row][col].p
            if len(ptrs):
                for p in ptrs:
                    stack.append([p.n, n, alen])

            else:
                # if there are no pointers, it is the end of a path
                save.append([a1[:alen], a2[:alen]])

        return save

    def writeScoreMatrix(self, file, decimal=0, reverse=False, space=2):
        """-----------------------------------------------------------------------------------------
        Write out the score matrix in an aligned table. Could provide scoremax, but then it wouldn't
        work for a sub-table.

        :param file: open filehandle for output, think stdout
        :param decimal:int, number of digits past decimal point
        :return: float, maximum score
        -----------------------------------------------------------------------------------------"""
        # first find the largest value (maximum column width)
        score = self.score
        scoremax = 0.0
        for row in score:
            for col in row:
                scoremax = max(scoremax, col.score)

        fmt = '{{:>{}.{}f}}'.format(len(str(scoremax)) + space, decimal)
        smt = '{{:>{}s}}'.format(len(str(scoremax)) + space)

        if reverse:
            # up, left
            s1 = self.s1.seq
            s2 = self.s2.seq
            i0 = len(s1) - 1
            j0 = len(s2) - 1
            step = -1

            file.write(smt.format(' '))
            i = i0
            while i >= 0:
                file.write(smt.format(s1[i]))
                i += step
            file.write('\n')

            j = j0
            while j >= 0:
                file.write(smt.format(s2[j]))
                i = i0
                while i >= 0:
                    file.write(fmt.format(score[j][i].score))
                    i += step
                j += step

                file.write('\n')
        else:
            # down, right
            file.write(smt.format(' '))
            for c in self.s1.seq:
                file.write(smt.format(c))
            file.write('\n')

            s2 = self.s2.seq
            i = 0
            for row in score:
                file.write(smt.format(s2[i]))
                i += 1
                for col in row:
                    file.write(fmt.format(col.score))

                file.write('\n')

        return scoremax

    def matchString(self, a1, a2):
        """-----------------------------------------------------------------------------------------
        creat a string showing the match between the aligned sequences a1 and a2.

        :param a1: string, aligned sequence 1
        :param a2: string, aligned sequence 2
        :return: string, match string
        -----------------------------------------------------------------------------------------"""
        idchar = '|'
        simchar = ':'
        cmp = self.table
        a2i = self.a2i

        match = ''
        for i in range(len(a1)):
            c1 = a1[i]
            c2 = a2[i]
            if c1 == c2:
                match += idchar
            elif c1 == '.' or c2 == '.':
                match += ' '
            elif cmp[a2i[c1]][a2i[c2]] > 0:
                match += simchar
            else:
                match += ' '

        return match

    @staticmethod
    def n2pos(l1, n):
        """-----------------------------------------------------------------------------------------
        return the row and col corresponding to cell n

        :param l1: int, length of sequence 1 (col)
        :param n: int, cell n
        :return: int, int; row, col
        -----------------------------------------------------------------------------------------"""
        return (n - 1) // l1, (n - 1) % l1


# --------------------------------------------------------------------------------------------------
# testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    align = Alignment()


    align.s1 = Fasta()
    align.s2 = Fasta()

    # align.alphabet = 'ACGT'
    # align.identity(pos=3, neg=-3)

    # align.s1.seq = 'ACTTATCTTAT'
    # align.s1.seq = 'TATTCTATTCA'
    # align.s1.seq = 'TGGTATACTAT'
    # align.s1.seq = 'GATACTATCTA'

    # align.s2.seq = 'AGTATCATATT'
    # align.s2.seq = 'TTATACTATGG'
    # align.s2.seq = 'TACTATTTAGAT'
    # align.s2.seq = 'TTATACTATGA'


    # align.seqToInt()
    # bestscore, bestpos = align.globalBrute(-1, -1, nogap=False)
    # bestscore, bestpos = align.localBrute(-1, -1)
    # print('score: {} at {}\n'.format(bestscore, bestpos))
    # align.writeScoreMatrix(sys.stdout, reverse=True, space=3)
    # alignments = align.traceAll(bestpos)
    # for a in alignments:
    #     m = align.matchString( a[0], a[1])
    #     print('{}\n{}\n{}\n\n'.format(a[0], m, a[1]))

    # align.readNCBI('..//scoring_matrices/alphabet.matrix')
    # align.s1.seq = 'BORROW'
    # align.s2.seq = 'BORABORA'
    align.readNCBI('..//scoring_matrices/dna4-2.matrix')
    # align.s1.seq = "TAGATTTATCAT"
    align.s2.seq = "TAGATTTATCAT"
    # align.s2.seq = "TACTATTTAGAT"
    align.s1.seq = "TGGTATACTAT"

    align.seqToInt()
    # bestscore, bestpos = align.globalBrute(-1, -1, nogap=False)
    bestscore, bestpos = align.localBrute(-1, -1)
    print('score: {} at {}\n'.format(bestscore, bestpos))
    align.writeScoreMatrix(sys.stdout, reverse=False, space=1)

    # testing
    # align.s1 = Fasta()
    # align.s1.seq = 'ACTGCCTTGATC'
    # align.s2 = Fasta()
    # align.s2.seq = 'ATGCCAAAGATC'
    # align.readNCBI('../../dotplot/table/NUC4.4.matrix')

    # align.seqToInt()
    # # random.shuffle(align.i1)          # uncomment to test scores for random alignments
    # original_score, bestpos = align.localBrute(-1, -1)
    # print('original score: {} at {}\n'.format(original_score, bestpos))
    # align.writeScoreMatrix(sys.stdout)
    # a1, a2, m = align.trace1(bestpos)
    # print('\n{}\n{}\n{}'.format(a1, m, a2))
    # alignments = align.traceAll(bestpos)
    # for a in alignments:
    #     m = align.matchString( a[0], a[1])
    #     print('{}\n{}\n{}\n\n'.format(a[0], m, a[1]))

    #
    # align.s1 = Fasta(filename=sys.argv[1])
    # align.s2 = Fasta(filename=sys.argv[2])
    # align.readNCBI('../../dotplot/table/BLOSUM62.matrix')

    exit(0)
