from sequence.fasta import Fasta
import re

# regular expression for parsing documentation
lenre  = re.compile('len=(\d+)')
pathre = re.compile('.*path=\[([^\]]*)\]')
idre   = re.compile(r'>*TRINITY_([^_]+)_c(\d+)_g(\d+)_i(\d+)')

class Trinity(Fasta):
    """-----------------------------------------------------------------------------------------------------------------
    Trinity
    Read multiple trinity transcript file in fasta format

    usage
        trinity = Trinity()
        trinity.fh = open(file, 'r')
        while trinity.next():
            ...
    -----------------------------------------------------------------------------------------------------------------"""

    def __init__(self):
        """'
        class constructor.
        Trinity is a subclass of Fasta
        """
        super().__init__()
        self.len = 0
        self.path = []
        self.cluster = ''
        self.component = 0
        self.gene = 0
        self.isoform = 0
        self.shortid = ''

        return None

    def getLen(self):
        """-----------------------------------------------------------------------------------------------------------------
        get the sequence length from the documentation
        :return: length fread from length field in documentation
        -----------------------------------------------------------------------------------------------------------------"""
        self.len = int(lenre.match(self.doc).group(1))
        return self.len

    def getPath(self):
        """-----------------------------------------------------------------------------------------------------------------
        The path describes how the predicted transcript is built from segments
        :return: list of path components from documentation
        -----------------------------------------------------------------------------------------------------------------"""
        path = pathre.match(self.doc).group(1)
        plist = path.split(' ')
        self.path = plist

        return plist


    def getIDParts(self):
        """-----------------------------------------------------------------------------------------------------------------
        Breakdown the trinity ID string to give the separate parts of the ID
        Cluster,  component, gene and isoform
        :return: cluster,  component, gene, isoform
         ----------------------------------------------------------------------------------------------------------------"""
        cluster, component, gene, isoform = idre.match(self.id).groups()
        self.cluster   = cluster
        self.component = int(component)
        self.gene      = int(gene)
        self.isoform   = int(isoform)
        self.shortid   = '{cl}.{co}.{g}.{i}'.format(cl=cluster, co=component, g=gene, i=isoform)

        return self.shortid

    def next(self):
        """-------------------------------------------------------------------------------------------------------------
        Overrides fasta method to add additional attributes of Trinity class.
        :return: True/False
        -------------------------------------------------------------------------------------------------------------"""
        if super().next():
            self.getLen()
            self.getIDParts()
            self.getPath()
            return True
        else:
            return False


if __name__ == '__main__':

    trinity = Trinity()

    #file = r'C:\Users\gribs\Dropbox\rice\Trinity.fasta'
    file = r'A:\mrg\repos\biocomputing\data\Trinity.fasta'
    trinity.fh = open(file, 'r')

    nseq = 0
    while trinity.next():
        nseq += 1

        trinity.doc = 'len={}'.format(trinity.len)
        trinity.id = trinity.shortid
        print( trinity.format())

        if nseq > 10: break

