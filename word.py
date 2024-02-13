"""=================================================================================================
word composition of a string

Michael Gribskov     02 June 2021
================================================================================================="""
import copy
from math import log2, log
from itertools import product


class Word:

    def __init__(self, size=6, alphabet='ACGT', sequence=''):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        self.size = size
        self.initialized = False
        self.alphabet = alphabet
        self.sequence = sequence
        self.count = {}
        self.count_total = 0

    def clone(self):
        """-----------------------------------------------------------------------------------------
        clone an existing instance of Word

        :return: instance of Word
        -----------------------------------------------------------------------------------------"""
        return copy.deepcopy(self)

    def counts_get(self, step=0):
        """-----------------------------------------------------------------------------------------
        count the words of length self.size, assumes word dict is initialized with all possible
        words, automatically initializes with zero if not previously done

        :return:int, number of words counted
        -----------------------------------------------------------------------------------------"""
        if not self.initialized:
            self.init_words()

        if step == 0:
            step = self.size

        n = 0
        for pos in range(0, len(self.sequence) - self.size + 1, step):
            # print(self.sequence[pos: pos + self.size])
            word = self.sequence[pos: pos + self.size]
            # if 'TAA' in word:
            #     print('oops')
            self.count[word] += 1
            self.count_total += 1
            n += 1

        return n

    def counts_get_sparse(self):
        """-----------------------------------------------------------------------------------------
        count the words of length self.size. works with an uninitialized list of words
        :return:
        -----------------------------------------------------------------------------------------"""

        n = 0
        for pos in range(0, len(self.sequence) - self.size + 1, self.size):
            # print(self.sequence[pos: pos + self.size])
            word = self.sequence[pos: pos + self.size]
            if word in self.count:
                self.count[word] += 1
            else:
                self.count[word] = 1
            self.count_total += 1

            n += 1

        return n

    def init_words(self, initval=0):
        """-----------------------------------------------------------------------------------------
        initialize the word list as the number of permutations of the alphabet

        :param initval: int, initial value for counts, e.g., a1 for a plus one prior
        :return: int, total count of words
        -----------------------------------------------------------------------------------------"""
        n = len(self.alphabet) ** self.size * initval
        for w in product(self.alphabet, repeat=self.size):
            self.count[''.join(w)] = initval
        self.count_total = n

        self.initialized = True

        return n

    def remove_stops_inframe(self):
        """-----------------------------------------------------------------------------------------
        Remove in frame stop codons.  stops are in frame if their position in the word  %3 == 0

        :return: int, remaining words
        -----------------------------------------------------------------------------------------"""
        stopcodons = ['TAA', 'TAG', 'TGA']
        n = 0
        delkeys = []
        for word in self.count:
            is_stop = False
            for stop in stopcodons:
                if word.find(stop) % 3 == 0:
                    delkeys.append(word)
                    is_stop = True
                    break
            if is_stop:
                continue

            n += 1

        for word in delkeys:
            self.count_total -= self.count[word]
            del self.count[word]

        return n

    def count_random(self, composition):
        """-----------------------------------------------------------------------------------------
        recalculate counts  based on iid model and letter composition.
        Counts are only calculated for existing keys in self.count and the total count of the
        random count is equal to the original total

        :param composition: dict, keys are letters in alphabet plus 'total'
        :return: int, total count
        -----------------------------------------------------------------------------------------"""
        new_total = 0

        # convert letter composition to frequency
        freq = {}
        for letter in composition:
            if letter == 'total':
                continue

            freq[letter] = composition[letter] / composition['total']

        for word in self.count:
            p = 1.0
            for letter in word:
                p *= freq[letter]
            self.count[word] = p
            new_total += self.count[word]

        self.count_total = new_total
        return new_total

    def composition_from_count(self):
        """-----------------------------------------------------------------------------------------
        Calculate the letter composition from the word count, include the total count with the
        key 'total'

        :return: dict, count or sum of frequencies of letters
        -----------------------------------------------------------------------------------------"""
        letter = {c: 0 for c in self.alphabet}
        letter['total'] = 0

        for word in self.count:
            for w in word:
                letter[w] += self.count[word]
                letter['total'] += self.count[word]

        return letter

    def to_frequency(self):
        """-----------------------------------------------------------------------------------------
        Convert counts to frequency by dividing by self.count_total

        :return: int, new total_count
        -----------------------------------------------------------------------------------------"""
        model = self.clone()
        model.count_total = 0
        for word in model.count:
            model.count[word] = model.count[word]/self.count_total
            model.count_total += model.count[word]

        return model

    def to_logfrequency(self):
        """-----------------------------------------------------------------------------------------
        Convert counts to frequency by dividing by self.count_total, the take log

        :return: Word object
        -----------------------------------------------------------------------------------------"""
        model = self.clone()
        model.count_total = 0
        for word in model.count:
            model.count[word] = log(model.count[word]/self.count_total)
            model.count_total += model.count[word]

        return model


    def add(self, count):
        """-----------------------------------------------------------------------------------------
        Add the counts in count, a Word object, to this count.

        :param count: Word object
        :return: int, total count
        -----------------------------------------------------------------------------------------"""
        for word in self.count:
            self.count[word] += count.count[word]
            self.count_total += count.count[word]

        return self.count_total

    def add_dict(self, countdict):
        """-----------------------------------------------------------------------------------------
        Add the counts in count, a dict indexed by sequence name, to this count.

        :param count: Word object
        :return: int, total count
        -----------------------------------------------------------------------------------------"""
        for word in self.count:
            self.count[word] += countdict[word]
            self.count_total += countdict[word]

        return self.count_total

# --------------------------------------------------------------------------------------------------
# testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    seq = 'GATGTCAAGGAAGATGAAGATGGACATGC'
    hexamer = Word(sequence=seq)
    c = hexamer.init_words(1)
    nwords = hexamer.remove_stops_inframe()
    c += hexamer.counts_get()
    print(f'total count:{hexamer.count_total}')
    for word in sorted(hexamer.count, key=lambda k: hexamer.count[k]):
        print(f'{word}\t{hexamer.count[word]}')
    hexamer_freq = hexamer.to_frequency()

    random = hexamer.clone()
    random.count_random(hexamer.composition_from_count())
    print(f'total count:{random.count_total}')
    for word in sorted(random.count, key=lambda k: random.count[k]):
        print(f'{word}\t{random.count[word]}')
    random_freq = random.to_frequency()

    logodds = {}
    for word in hexamer.count:
        logodds[word] = log2( hexamer_freq[word] / random_freq[word])

    for word in sorted(logodds, key=lambda k: logodds[k]):
        print(f'{word}\t{logodds[word]:.2f}')

    print(seq)

    exit(0)
