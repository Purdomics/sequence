"""=================================================================================================
extract protein coding sequences from genbank nucleotide file.  Designed for mitochondrial genome

Michael Gribskov     15 October 2022
================================================================================================="""
import sys


def get_id(file, taxa_names):
    """---------------------------------------------------------------------------------------------
    get a unique ID from the description line

    :param taxa_names: list     names of taxa
    :param file: filehandle     source genbank file
    :return: logical            True if id could be read, false at EOF
    ---------------------------------------------------------------------------------------------"""
    remove = 'isolate breed voucher strain: mitochondrion complete partial genome almost mitochondrial Hybrid horse historical'
    # punc = {'.': '', ',': '', '\n':''}
    punc = ''.maketrans('\\/:_', '___ ', '.,())')

    desc = ''
    desc_found = None
    for line in file:
        if line.startswith('DEFINITION'):
            # including capitals in the remove list conflicts with haplogroups
            line = line.replace('DEFINITION ', '').rstrip()
            line = line.replace('DNA', '')
            # print(f'\n{line.rstrip()}')
            desc = line
            desc_found = True
            break
        desc_found = False

    filtered = []
    if desc_found:
        # desc = desc.replace('_', ' ')
        words = desc.translate(punc).split()
        for w in words:
            if w not in remove:
                filtered.append(w)

    id = '_'.join(filtered)
    id = id.replace('haplogroup_', 'hap')
    id = id.replace('haplotype_', 'hap')

    if id in taxa_names:
        base, n = id.rsplit('_', maxsplit=1)
        if n.isdigit():
            id = f'{base}_{int(n) + 1}'
        else:
            id = f'{id}_1'

    # print(f'{id}\t\t{line}')
    taxa_names.append(id)

    return desc_found


def get_cds(file, taxa_names, gene_count, gene_seq):
    """---------------------------------------------------------------------------------------------
    
    :param gb: 
    :param taxa_names: 
    :param gene_count: 
    :param gene_seq: 
    :return: 
    ---------------------------------------------------------------------------------------------"""
    gene_synonyms = {'COXI':    'COX1',
                     'COXII':   'COX2',
                     'COXIII':  'COX3',
                     'COI':     'COX1',
                     'COII':    'COX2',
                     'COIII':   'COX3',
                     'cox1':    'COX1',
                     'NADH1':   'ND1',
                     'NADH2':   'ND2',
                     'ATPase8': 'ATP8',
                     'ATPase6': 'ATP6',
                     'NADH3':   'ND3',
                     'NADH4L':  'ND4L',
                     'NADH4':   'ND4',
                     'NADH5':   'ND5',
                     'NADH6':   'ND6',
                     'cytb':    'CYTB',
                     'CTYB':    'CYTB'
                     }

    cdsfound = False
    transfound = False
    cds = ''
    gene = ''
    id = taxa_names[-1]
    for line in file:
        line = line.rstrip()
        if line.startswith('ORIGIN'):
            # end of feature table
            break

        if line.startswith('     CDS'):
            cdsfound = True

        if cdsfound:
            if line.startswith('                     /gene='):
                # /gene="CYTB"
                gene = line.replace('                     /gene=', '').replace('"', '').rstrip()
                if gene in gene_synonyms:
                    gene = gene_synonyms[gene]

                continue

            if line.startswith('                     /translation'):
                transfound = True

            if transfound:
                cds += line.strip()
                if line.endswith('"'):
                    # end of translation
                    cds = cds.replace('/translation=', '')
                    cds = cds.replace('"', '')
                    if 'X' not in cds:
                        if gene in gene_count:
                            gene_count[gene] += 1
                            gene_seq[gene][id] = cds
                        else:
                            gene_count[gene] = 1
                            gene_seq[gene] = {id:cds}

                    transfound = False
                    cdsfound = False
                    cds = ''
                    gene = ''

    return


def get_dna(gb):
    """---------------------------------------------------------------------------------------------
    read the nucleotide sequence and calculate checksum (sum of characters)

    :param gb: filehandle   open multiple genbank file
    :return: int            checksum
    ---------------------------------------------------------------------------------------------"""
    seq = ''
    for line in gb:
        if line.startswith('//'):
            break

        field = line.strip().split()
        seq += field[1]

    # calculate checksum
    checksum = 0
    for a in seq:
        checksum += ord(a)

    return checksum


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    gb = open(sys.argv[1], 'r')
    taxa_names = []
    gene_count = {}
    gene_seq = {}
    checksums = []
    duplicate = []

    unique = 0
    while get_id(gb, taxa_names):
        while get_cds(gb, taxa_names, gene_count, gene_seq):
            pass
        checksum = get_dna(gb)
        # print(f'checksum: {checksum}')

        if checksum in checksums:
            # duplicate sequence
            duplicate.append(taxa_names[-1])
            print(f'{taxa_names[-1]} is a duplicate')
        else:
            checksums.append(checksum)
            unique += 1

    print(f'{unique} unique')

    for c in sorted(gene_count, key=lambda x: gene_count[x], reverse=True):
        print(f':{c}:\t{gene_count[c]}')

    longcount = 0
    long = open('../data/long.fa', 'w')
    for t in taxa_names:
        longseq=''
        missing = False
        for g in gene_seq:
            # print(f'{t}\t{gene_seq[g][t]}')
            if t in gene_seq[g]:
                longseq += gene_seq[g][t]
            else:
                missing = True
                break
        if missing:
            continue
        else:
            longcount += 1
            long.write(f'>{t}\n{longseq}\n')

    print( f'{longcount} long sequences')
    long.close()

    exit(0)
