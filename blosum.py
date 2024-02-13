"""=================================================================================================
perform blosum calculation on example data
could be extended to full calculation

Michael Gribskov     09 March 2022
================================================================================================="""
from math import log


class Frequency:
    """=============================================================================================
    Frequency is a vector of counts.  Counts may be either integer counts or probabilities. The
    alphabet is determined by the data

    ============================================================================================="""

    def __init__(self, alphabet=''):
        """-----------------------------------------------------------------------------------------
        alphabet: string, letters in composition
        c: count for each letter
        n: float/int, total count of letters
        count: float/int, original letter count (to allow back transformation
        -----------------------------------------------------------------------------------------"""
        self.c = {}
        self.n = 0
        self.count = 0
        self.alphabet = ''
        if alphabet:
            self.alphabet = alphabet
            self.c = {aa: 0.0 for aa in alphabet}

    @classmethod
    def composition_ref(cls, select, frequency=True):
        """-----------------------------------------------------------------------------------------
        Alternate constructor, return composition of protein database

        Data from:http://Martz.MolviZ.Org,  April, 2020 by Eric Martz,

        Sources

        "creighton": Thomas E. Creighton (EMBO Heidelberg), Proteins, Structures and Molecular Properties, 2nd ed. 1993,
        W. H. Freeman and Co. Page 5:
        Table of amino acid residue masses and frequencies in 1,021 unrelated proteins of known sequence

        "trinquier": Trinquier G, Sanejouand YH. Which effective property of amino acids is best preserved by the genetic
        code?. Protein Eng. 1998;11(3):153–169. doi:10.1093/protein/11.3.153
        https://pubmed.ncbi.nlm.nih.gov/9613840/
        values from this table: http://proteinsandproteomics.org/content/free/tables_1/table08.pdf
        "105,990 sequences in the nonredundant OWL protein database (release 26.0 e) "

        "carugo": Carugo, Oliviero. “Amino acid composition and protein dimension.”
        Protein science 17,12 (2008): 2187-91. doi:10.1110/ps.037762.108
        Values read from Figure 2 for lengths of 200 residues.
        Derived from uniref50 FASTA file from UniProt, sequence identity <50%. About 550,000 sequences, which give several
        thousand sequences at each length.
        mrg: i assume carugo100 is for proteins > 100 residues

        "swiss2022":UniProtKB/Swiss-Prot protein knowledgebase release 2022_01 statistics

        :param select: string, key of composition to select
        :param frequency: boolean, if True convert to frequency
        :return: dict, keys are amino acids (one letter code)
        -----------------------------------------------------------------------------------------"""
        data = {'creighton': {'A': 8.3, 'C': 1.7, 'D': 5.3, 'E': 6.2, 'F': 3.9,
                              'G': 7.2, 'H': 2.2, 'I': 5.2, 'K': 5.7, 'L': 9.0,
                              'M': 2.4, 'N': 4.4, 'P': 5.1, 'Q': 4.0, 'R': 5.7,
                              'S': 6.9, 'T': 5.8, 'V': 6.6, 'W': 1.3, 'Y': 3.2},
                'trinquier': {'A': 7.5, 'C': 1.8, 'D': 5.2, 'E': 6.3, 'F': 3.9,
                              'G': 7.1, 'H': 2.2, 'I': 5.5, 'K': 5.8, 'L': 9.1,
                              'M': 2.8, 'N': 4.6, 'P': 5.1, 'Q': 4.1, 'R': 5.2,
                              'S': 7.4, 'T': 6.0, 'V': 6.5, 'W': 1.3, 'Y': 3.3},
                'carugo200': {'A': 8.7, 'C': 1.5, 'D': 5.1, 'E': 6.2, 'F': 4.0,
                              'G': 6.8, 'H': 2.2, 'I': 5.7, 'K': 5.3, 'L': 9.8,
                              'M': 2.4, 'N': 3.9, 'P': 5.0, 'Q': 3.9, 'R': 6.2,
                              'S': 7.0, 'T': 5.3, 'V': 6.5, 'W': 1.3, 'Y': 2.9},
                'carugo100': {'A': 8.5, 'C': 1.9, 'D': 4.7, 'E': 6.0, 'F': 4.0,
                              'G': 6.5, 'H': 2.3, 'I': 5.8, 'K': 5.6, 'L': 9.6,
                              'M': 2.8, 'N': 3.7, 'P': 4.9, 'Q': 3.9, 'R': 6.7,
                              'S': 7.1, 'T': 5.2, 'V': 6.4, 'W': 1.4, 'Y': 2.7},
                'swiss2022': {'A': 8.25, 'C': 1.38, 'D': 5.46, 'E': 6.72, 'F': 3.86,
                              'G': 7.07, 'H': 2.27, 'I': 5.91, 'K': 5.80, 'L': 9.65,
                              'M': 2.41, 'N': 4.06, 'P': 4.74, 'Q': 3.93, 'R': 5.53,
                              'S': 6.64, 'T': 5.35, 'V': 6.86, 'W': 1.10, 'Y': 2.92, }
                }

        f = cls()

        for aa in data[select]:
            f.c[aa] = data[select][aa]
            f.n += f.c[aa]
            f.alphabet += aa

        f.count = f.n
        if frequency:
            f.divide()

        return f

    def divide(self, factor=None):
        """-----------------------------------------------------------------------------------------
        Divide all values by a factor

        :param factor: float or int, denominator for scaling values
        :return: float or int, total count after normalization
        -----------------------------------------------------------------------------------------"""
        if factor:
            denom = factor
        else:
            denom = self.n

        self.n = 0.0
        for aa in self.c:
            self.c[aa] /= denom
            self.n += self.c[aa]

        return self.n


class Transition:
    """=============================================================================================
    Transition matrices stored as dict of dict. Transition is the upper triangular matrix.

    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        self.alphabet = ''
        self.type = ''
        self.frequency = {}
        self.total = 0.0

    @classmethod
    def from_alphabet(cls, alpha):
        """-----------------------------------------------------------------------------------------
        Alternate constructor: Set up an upper triangular transition matrix with zero values based
        on alpha

        :param alpha: string, possible letters in use
        :return: Transition object
        -----------------------------------------------------------------------------------------"""
        transition = cls()
        al = len(alpha)
        transition.frequency = {alpha[i]: {alpha[j]: 0 for j in range(i, al)} for i in range(al)}
        transition.alphabet = alpha

        return transition

    @classmethod
    def from_composition(cls, expected, n=1.0):
        """-----------------------------------------------------------------------------------------
        Alternate constructor - calculate the expected transition frequencies based on the expected
        count (background).

        t(i,i) = exp(i)**2
        t(i,j) = s * exp(i) * exp(j)

        :param expected: dict, expected frequencies indexed by aa
        :parameter n: int, expected number of transitions
        :return: Transition object
        -----------------------------------------------------------------------------------------"""
        alphabet = expected.alphabet
        back = cls.from_alphabet(alphabet)
        f = expected.c

        transition = back.frequency
        total = 0
        for i in range(len(alphabet)):
            aa_i = alphabet[i]
            for j in range(i, len(alphabet)):
                aa_j = alphabet[j]
                if i == j:
                    transition[aa_i][aa_j] = n * f[aa_i] ** 2
                else:
                    transition[aa_i][aa_j] = n * 2 * f[aa_i] * f[aa_j]

                total += transition[aa_i][aa_j]

        back.total = total
        return back

    @classmethod
    def from_block(cls, block):
        """---------------------------------------------------------------------------------------------
        calculate transition counts across all columns of the block

        :param block: Block object, list of dict with counts at each position
        :return: Transition object
        ---------------------------------------------------------------------------------------------"""
        alphabet = block.alphabet
        transition = cls.from_alphabet(alphabet)
        frequency = transition.frequency

        for pos in range(block.length):
            c = block.composition[pos]
            letter = sorted(list(c.keys()))
            for i in range(len(letter)):
                for j in range(i, len(letter)):
                    ci = letter[i]
                    cj = letter[j]
                    if ci == cj:
                        frequency[ci][ci] += (c[ci] * (c[ci] - 1)) / 2
                    else:
                        frequency[ci][cj] += c[ci] * c[cj]
                        # transition[cj][ci] += c[ci] * c[cj]

        return transition

    def addtransition(self, prior, factor = 1.0, frequency=True):
        """-----------------------------------------------------------------------------------------
        add the frequencies/counts in prior to each position of the composition

        :param prior: Transition object
        :param factor: int/float, multiple prior by this factor
        :param frequency: boolean, if True divide by total counts to get frequencies (probabilities)

        :return: count at last position (all position should be the same)
        -----------------------------------------------------------------------------------------"""
        composition = self.frequency
        bkg = prior.frequency

        self.total = 0.0
        for row in composition.keys():
            for col in composition[row]:
                composition[row][col] += factor * bkg[row][col]
                self.total += composition[row][col]

        return self.total

    def toprobability(self):
        """-----------------------------------------------------------------------------------------
        convert to probability by dividing by self.total
        :return: float, sum of all transitions
        -----------------------------------------------------------------------------------------"""
        composition = self.frequency

        total = self.total
        self.total = 0
        for row in composition.keys():
            for col in composition[row]:
                composition[row][col] /= total
                self.total += composition[row][col]

        return self.total


class Block:
    """=============================================================================================
    for Blocks motifs, or other short gapless alignments. Alignments are a dictionary with the
    sequence names as keys
    ============================================================================================="""

    def __init__(self, block=None, title=''):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        self.title = ''
        self.block = {}
        self.weight = {}
        self.length = 0
        self.composition = []
        self.alphabet = ''  # currently set by composition()
        self.n = 0
        if block:
            self.block = block
            self.setweight()
            self.getn()
            self.getlength()

        if title:
            self.title = title

    # def setalphabet(self):
    #     """-----------------------------------------------------------------------------------------
    #     set the alphabet string based on the current positional composition
    #
    #     :return: string, alphabet used in block
    #     -----------------------------------------------------------------------------------------"""

    def setweight(self, value=1.0):
        """-----------------------------------------------------------------------------------------
        Set all weights to value

        :param value: float, value for weights
        :return: dict, key=name of sequence, value = 1
        -----------------------------------------------------------------------------------------"""
        self.weight = {k: value for k in self.block.keys()}
        return self.weight

    def getn(self):
        """-----------------------------------------------------------------------------------------
        Set the number of sequences in the block and return the number

        :return:
        -----------------------------------------------------------------------------------------"""
        n = 0
        for id in self.block:
            n += 1

        self.n = n
        return n

    def getlength(self):
        """-----------------------------------------------------------------------------------------
        Sets the block length in the object and returns length of sequence in block. Length is the
        minimum sequence length of all the sequences

        :return: int, block length
        -----------------------------------------------------------------------------------------"""
        block = self.block
        seqid = block.keys()
        length = None

        for id in seqid:
            if length:
                length = max(length, len(block[id]))
            else:
                length = len(block[id])

        self.length = length
        return length

    def composition_count(self):
        """---------------------------------------------------------------------------------------------
        Calculate the composition of a block

        :return: dict, overall composition, keys are amino acids
        ---------------------------------------------------------------------------------------------"""
        block = self.block
        weight = self.weight
        composition = self.composition

        overall = {}
        for pos in range(self.length):
            composition.append({})
            for seq in block:
                letter = block[seq][pos]
                if letter not in self.alphabet:
                    self.alphabet += letter

                try:
                    composition[pos][letter] += weight[seq]
                except KeyError:
                    composition[pos][letter] = weight[seq]

                try:
                    overall[letter] += weight[seq]
                except KeyError:
                    overall[letter] = weight[seq]

        return overall


def tabular(table, format):
    """---------------------------------------------------------------------------------------------

    values are printed with f format, labels with s format

    :param table: dict of dict, keys are aa, values assumed to be float
    :param float: string, format for values, e.g., 5.2
    :return:
    ---------------------------------------------------------------------------------------------"""
    out = ''

    if format.find('.'):
        width, precision = format.split('.')
        width = int(width)
    fs = f'{width}s'
    fv = format + 'f'

    # column labels
    column = f'{" ":>{fs}}'
    for aa in table.keys():
        column += f'{aa:>{fs}}'

    out += f'{column}\n'

    for i in table:
        out += f'{i:<{fs}}'
        for j in table:
            try:
                out += f'{table[i][j]:{fv}}'
            except TypeError:
                # capture the Nones
                a = 'NA'
                out += f'{a:>{fs}}'
            except KeyError:
                # KeyError means indices are reversed
                out += f'{table[j][i]:{fv}}'

        out += '\n'

    return out


def logodds(foreground, background, base=2):
    """---------------------------------------------------------------------------------------------
    calculate the symmetric log odds values from two Transition objects.  values that cannot be
    calculated due to numeric errors such as divide by zero or log(0) will be returned as None

    :param foreground: Transition object
    :param background: Transition object
    :return: dict of dict
    ---------------------------------------------------------------------------------------------"""
    alpha = background.alphabet
    fore = foreground.frequency
    back = background.frequency
    logodds = {a1: {a2: None for a2 in alpha} for a1 in alpha}

    for i in range(len(alpha)):
        a1 = alpha[i]
        for j in range(i, len(alpha)):
            a2 = alpha[j]
            try:
                odds = fore[a1][a2] / back[a1][a2]
                if odds == 0:
                    continue
                else:
                    # make symmetric
                    logodds[a1][a2] = log(odds, base)
                    logodds[a2][a1] = log(odds, base)

            except ZeroDivisionError:
                continue

    return logodds


def transition_subset(trans, subset):
    """---------------------------------------------------------------------------------------------
    return a subset of the transition matrix with only the rows and columns in the string subset
    
    :param trans: dict of dict, transition values
    :param subset: string, list of desired keys (amino acids)
    :return: dict of dict, smaller transition matrix
    ---------------------------------------------------------------------------------------------"""
    smaller = {a1: {a2: 0 for a2 in subset} for a1 in subset}

    for a1 in subset:
        for a2 in subset:
            smaller[a1][a2] = trans[a1][a2]

    return smaller


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    # background transition distribution
    exp = Frequency.composition_ref('swiss2022', frequency=True)
    background = Transition.from_composition(exp)

    # METHYLTRANSFERASE BI
    block = {'TCMN_STRGA': 'IADLGGGDGWFLAQILRRHPHATGLLMDLPRVA',
             'TCMO_STRGA': 'FVDLGGARGNLAAHLHRAHPHLRATCFDLPEME',
             'ZRP4_MAIZE': 'LVDVGGGIGAAAQAISKAFPHVKCSVLDLAHVV',
             'CHMT_POPTM': 'LVDVGGGTGAVVNTIVSKYPSIKGINFDLPHVI',
             'COMT_EUCGU': 'VVDVGGGTGAVLSMIVAKYPSMKGINFDLPHVI',
             'COMT_MEDSA': 'LVDVGGGTGAVINTIVSKYPTIKGINFDLPHVI',
             'CRTF_RHOSH': 'LMDVGGGTGAFLAAVGRAYPLMELMLFDLPVVA',
             'OMTA_ASPPA': 'VVDVGGGRGHLSRRVSQKHPHLRFIVQDLPAVI', }

    mat = Block(block, title='Methyltransferase B1')
    mat.alphabet = background.alphabet
    mat.composition_count()
    expected_n = mat.n * (mat.n - 1) * mat.length / 2

    foreground = Transition.from_block(mat)
    foreground.addtransition(background, factor=mat.length, frequency=True)
    foreground.toprobability()
    # trans = transition_count(composition)

    lo = logodds(foreground, background)

    print(tabular(foreground.frequency, '5.1'))
    print(tabular(background.frequency, '5.1'))
    print(tabular(lo, '5.1'))

    # strans = transition_subset(trans,'FILV')
    # print(tabular(strans, '5.1'))
    # sbackground = transition_subset(background, 'FILV')
    # print(tabular(sbackground, '5.1'))
    # slogodds = transition_subset(logodds, 'FILV')
    # print(tabular(slogodds, '5.1'))

    # for column 1
    blocklen = 1
    expected_n = nseq * (nseq - 1) * blocklen / 2
    exp_sub = {}
    for aa in 'FILV':
        exp_sub[aa] = exp[aa]

    background = transition_background(exp_sub, expected_n)
    composition, count = composition_count(block, blocklen, weight)
    trans = transition_count(composition)
    for a1 in background.keys():
        for a2 in background.keys():
            background[a1][a2] *= expected_n

    strans = transition_subset(trans, 'FILV')
    sbackground = transition_subset(background, 'FILV')
    lo = logodds(strans, sbackground)
    slo = transition_subset(lo, 'FILV')

    print('\ntransition')
    print(tabular(strans, '7.3'))
    print('\nbackground')
    print(tabular(sbackground, '7.3'))
    print('\nlog-odds')
    print(tabular(slo, '7.3'))

    exit(0)
