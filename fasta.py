"""-------------------------------------------------------------------------------------------------
Fasta sequence class.  Supports iteration over a multi-fasta file
    filename
    id
    documentation
    sequence

    Synopsis
    from sequence.fasta import Fasta

    fasta = Fasta()
    fasta.open('filename')
    while fasta.next():
        print(fasta.format(linelen=60))

-------------------------------------------------------------------------------------------------"""
import sys


# noinspection PyMethodParameters
class Fasta:
    """=============================================================================================
    Read, write, and modify FastA formatted sequence files. This class is iterable.

    Fasta format:
    >sequence_id documentation_continues_till_newline
    sequence
    sequence
    ...
    >next_sequence_id ...

    ============================================================================================="""
    codon2aa = {"AAA": "K", "AAC": "N", "AAG": "K", "AAT": "N",
                "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
                "AGA": "R", "AGC": "S", "AGG": "R", "AGT": "S",
                "ATA": "I", "ATC": "I", "ATG": "M", "ATT": "I",

                "CAA": "Q", "CAC": "H", "CAG": "Q", "CAT": "H",
                "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
                "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
                "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",

                "GAA": "E", "GAC": "D", "GAG": "E", "GAT": "D",
                "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
                "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
                "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",

                "TAA": "*", "TAC": "Y", "TAG": "*", "TAT": "T",
                "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
                "TGA": "*", "TGC": "C", "TGG": "W", "TGT": "C",
                "TTA": "L", "TTC": "F", "TTG": "L", "TTT": "F"}

    complement = str.maketrans('ACGTUacgtu', 'TGCAAtgcaa')

    def __init__(fasta, filename="", mode='r', fh=None):
        """-----------------------------------------------------------------------------------------
        Fasta class constructor. Attributes
            filename
            fh          filehandle
            id
            doc
            seq
            buffer  (read ahead buffer, only internal) (deprecated, not needed for iterable)

            TODO probably should create both read and write filehandles
        -----------------------------------------------------------------------------------------"""
        fasta.filename = ''
        fasta.id = ''
        fasta.doc = ''
        fasta.seq = ''
        fasta.buffer = ''
        fasta.fh = None

        if fh:
            fasta.fh = fh
        elif filename:
            fasta.fh = fasta.open(filename, mode)
            next(fasta)

    def __next__(self):
        return self

    def get_entry(self):
        """-----------------------------------------------------------------------------------------
        read an entry from the open file handle
        updates self.id, self.doc, self.seq
        -----------------------------------------------------------------------------------------"""
        for line in self.fh:
            if line.isspace():
                continue

            if line.startswith('>'):
                # title line
                field = line.rstrip().split(' ', maxsplit=1)
                self.id = field[0].replace('>', '')
                if len(field) > 1:
                    self.doc = field[1]
                    self.seq = ''

            else:
                # sequence line
                self.seq += line.rstrip()



    def __iter__(self):
        """-----------------------------------------------------------------------------------------
        Generator function that reads the next sequence and updates internal variables
        :return: (string, string, string)   id, doc, and seq of newly read sequence
        -----------------------------------------------------------------------------------------"""
        for line in self.fh:
            if line.isspace():
                continue

            if line.startswith('>'):
                # title line
                if self.seq:
                    yield self

                field = line.rstrip().split(' ', maxsplit=1)
                self.id = field[0].replace('>', '')
                self.doc = ''
                if len(field) > 1:
                    self.doc = field[1]
                self.seq = ''

            else:
                # sequence line
                self.seq += line.rstrip()

        if self.seq:
            # if the loop ends with stored sequence, yield it here
            yield self

        # reset internal values so that the next call won't duplicate
        self.id = ''
        self.doc = ''
        self.seq = ''

        return

    def open(fasta, filename, mode='r'):
        """-----------------------------------------------------------------------------------------
        open a file for reading with error checking. If file can't be opened, exit with status=1.
        sets both self.filename and self.fh

        :param filename: string     path to file
        :param mode: string         'r' or 'w', any valid file open mode
        :return: filehandle
        -----------------------------------------------------------------------------------------"""
        # print('opening:', filename)
        fasta.fh = None
        fasta.filename = filename
        try:
            fasta.fh = open(filename, mode)
        except (OSError, IOError):
            print('Fasta::open - file open error')
            exit(1)

        return fasta.fh

    def next(fasta):
        """-----------------------------------------------------------------------------------------
        return the next entry from an open file into the object
        usage
            while fasta.next():
               ...
        -----------------------------------------------------------------------------------------"""
        a = []
        a = next(fasta)
        return next(fasta)

    def read(fasta):
        """-----------------------------------------------------------------------------------------
        read one sequence from the file, leave the following line in buffer
        usage:

        TODO now that the class is iterable, this isn't really needed, replace with call to __iter__()
        fasta.read()
        -----------------------------------------------------------------------------------------"""

        fasta.id = ''
        fasta.doc = ''
        fasta.seq = ''

        # get the ID and doc
        fasta.getID()

        for line in fasta.fh:
            if line.isspace():
                continue
            try:
                line = line.rstrip('\n\r ')
            except TypeError:
                # in case of a byte string, not an error if line is bytes
                line = line.decode()
                line = line.rstrip('\n\r ')

            if line[0] == '>':
                fasta.buffer = line
                break

            else:
                fasta.seq += line

        if len(fasta.id) > 0 or len(fasta.doc) > 0 or len(fasta.seq) > 0:
            return True

        # fall through to false if nothing can be read
        return False

    def copy(self):
        """-----------------------------------------------------------------------------------------
        Return a new Fasta object containing the current sequence

        :return: Fasta object
        -----------------------------------------------------------------------------------------"""
        new = Fasta()
        new.id = self.id
        new.doc = self.doc
        new.seq = self.seq

        return new

    def getID(fasta):
        """-----------------------------------------------------------------------------------------
        intended to be used internally for sequence reading
        check the buffer, and if not empty read the ID and documentation
        store in id and doc attributes
        id will be stripped of >
        documentation will be and empty string if there is nothing following the ID
        -----------------------------------------------------------------------------------------"""
        # if buffer is empty read a line
        if fasta.buffer:
            line = fasta.buffer
            fasta.buffer = ''
        else:
            line = fasta.fh.readline()
            try:
                line = line.rstrip('\n')
            except TypeError:
                # in case of a byte string
                # line = line.rstrip('\n'.encode())
                line = line.decode()
                line = line.rstrip('\n')

        # get the ID and documentation from the doc line
        try:
            seqid, doc = line.split(" ", 1)
        except ValueError:
            # documentation is missing
            seqid = line
            doc = ''

        fasta.id = seqid.lstrip('>')
        fasta.doc = doc

        return True

    def length(fasta):
        """-----------------------------------------------------------------------------------------
        return the length of the  current sequence
        return 0 if there is none
        usage
            seqlen = fasta.length()
        -----------------------------------------------------------------------------------------"""
        return len(fasta.seq)

    def format(fasta, linelen=50):
        """-----------------------------------------------------------------------------------------
        return a formatted string with the current sequence
        usage
            seq = fasta.format()
        -----------------------------------------------------------------------------------------"""
        string = '>{0} {1}'.format(fasta.id, fasta.doc)
        pos = 0
        while pos < len(fasta.seq):
            string += '\n{0}'.format(fasta.seq[pos:pos + linelen])
            pos += linelen

        return string

    def trimDocByRegex(fasta, target):
        """-----------------------------------------------------------------------------------------
        Shorten documentation by substituting the target regex with nothing
        target must be a compiled regex
        The new documentation string is returned
        usage
            trim = re.compile( 'len=\d+ ' )
            doc = fasta.trimDocAfterMatch( trim )
        -----------------------------------------------------------------------------------------"""
        if target:
            # skip if no regex is given
            fasta.doc = target.sub('', fasta.doc)

        return fasta.doc

    @staticmethod
    def reverseComplement(seq):
        """-----------------------------------------------------------------------------------------
        Return the sequence converted to reverse complement
        :return: string
        -----------------------------------------------------------------------------------------"""
        seq = seq.translate(Fasta.complement)

        return seq[::-1]

    def translate(fasta, frame=0, direction='+'):
        """-----------------------------------------------------------------------------------------
        translate in a nucleic acid sequence in the desired direction (+.-) and frame (0..2)
        incomplete codons at end are not translated
        stop codons are shown as '*'
        codons with ambiguity characters are translated as 'X';
        currently uses standard amino acid code
        TODO use supplied genetic code

        :param frame: integer, offset from beginning of sequence
        :param direction: string, forward (+) or reverse(-)
        :return: Fasta object (new)
        -----------------------------------------------------------------------------------------"""
        if not fasta.isACGT():
            sys.stderr.write('Fasta::translate - sequence must be ACGT')

        # rf is a suffix added to the ID to distinguish the reading frame
        rf = 'f0'
        if direction == '+':
            rf = '{}{}'.format(direction, frame)
        else:
            rf = 'r{}'.format(frame)

        trans = Fasta()
        trans.id = fasta.id + '_{}'.format(rf)
        trans.doc = fasta.id + ' strand={} frame={}'.format(direction, frame)

        seq = fasta.seq
        if direction == '-':
            seq = Fasta.reverseComplement(seq)

        pos = frame
        while pos < len(seq) - 2:
            codon = seq[pos:pos + 3]
            codon = codon.upper()
            # print('{}:{}:{}'.format(pos, codon, Fasta.codon2aa[codon]))
            trans.seq += Fasta.codon2aa[codon]
            pos += 3

        return trans

    def translate_all(self):
        """-----------------------------------------------------------------------------------------
        Return all six translated reading frames as a list of Fasta.  The starting positions of
        the reading frames are:
        0   0
        1   1
        2   2
        3 len-1
        4 len-2
        5 len-3

        :return: list of Fasta, the six reading frames translated
        -----------------------------------------------------------------------------------------"""
        rf = []
        for direction in ('+', '-'):
            for frame in range(3):
                rf.append = self.translate(frame, direction)

        return rf

    def composition(fasta, uppercase=False):
        """-----------------------------------------------------------------------------------------
        Returns a dictionary with the composition of the sequence
        If uppercase is true, characters are converted to uppercase

        :param fasta: self
        :param uppercase: boolean
        :return: dict, int; keys are sequence letters, values are counts
        -----------------------------------------------------------------------------------------"""
        seq = fasta.seq
        if uppercase:
            seq = fasta.seq.upper()

        count = {}
        for ch in seq:
            if ch in count:
                count[ch] += 1
            else:
                count[ch] = 1

        return count

    def frequency(fasta, uppercase=False):
        """-----------------------------------------------------------------------------------------
        Returns a dictionary with the composition of the sequence as frequencies (probabilities).
        Simply calls fasta.composition and divides the values by the character count. If uppercase
        is true, characters are converted to uppercase

        :param uppercase: boolean
        :return: dict, float; keys are sequence letters , values are frequencies
        -----------------------------------------------------------------------------------------"""

        count = fasta.composition(uppercase=uppercase)
        total = 0
        for a in count:
            total += count[a]
        for a in count:
            count[a] /= total

        return count

    def isACGT(fasta, threshold=0.8):
        """-----------------------------------------------------------------------------------------
        Return True if at least threshold fraction of characters in the sequence are ACGT

        :param threshold: float
        :return: Boolean
        -----------------------------------------------------------------------------------------"""
        total = fasta.length()
        if not total:
            return False

        comp = fasta.composition(uppercase=True)
        acgt = 0
        for base in 'ACGT':
            try:
                acgt += comp[base]
            except KeyError:
                # ignore missing bases
                continue

        if acgt / total >= threshold:
            return True

        return False


# --------------------------------------------------------------------------------------------------
# testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    fasta = Fasta()
    fasta.id = 'sample1'
    fasta.doc = '20 each A,C,G,T'
    fasta.seq = 'A' * 20 + 'C' * 20 + 'G' * 20 + 't' * 20

    print(fasta.format(40))

    print('\nComposition')
    comp = fasta.composition()
    for ch in comp:
        print('\t{}\t{}'.format(ch, comp[ch]))
    comp = fasta.composition(uppercase=True)
    print('Uppercase')
    for ch in comp:
        print('\t{}\t{}'.format(ch, comp[ch]))

    print('\nFrequencies (without uppercasing)')
    freq = fasta.frequency()
    for ch in freq:
        print('\t{}\t{}'.format(ch, freq[ch]))

    # test isACGT
    print('\nACGT should be true')
    print(fasta.format(80))
    print('ACGT:{}'.format(fasta.isACGT()))

    print('\nACGT should be true')
    fasta.doc = 'No T'
    fasta.seq = 'A' * 20 + 'C' * 20 + 'G' * 20
    print(fasta.format(80))
    print('ACGT:{}'.format(fasta.isACGT()))

    print('\nACGT should be false')
    fasta.doc = 'ABC - 20 each'
    fasta.seq = 'A' * 20 + 'B' * 20 + 'C' * 20
    print(fasta.format(80))
    print('ACGT:{}'.format(fasta.isACGT()))

    # translation
    print('\nTranslation')
    fasta.id = 'sample1'
    fasta.doc = '20 each A,C,G,T'
    fasta.seq = 'A' * 20 + 'C' * 20 + 'G' * 20 + 't' * 20

    for direction in 'rf':
        for frame in range(3):
            trans = fasta.translate(frame=frame, direction=direction)
            print(trans.format(80))

    exit(0)
