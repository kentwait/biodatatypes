# Naming
STANDARD_AA = [
    'Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Glu', 'Gln', 'Gly', 'His', 'Ile',
    'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val'
]
# One-letter lookup tables
FROM_AMINO_ACID_ONE_LETTER_TOKEN = {
    'A': 'Ala',
    'R': 'Arg',
    'N': 'Asn',
    'D': 'Asp',
    'C': 'Cys',
    'E': 'Glu',
    'Q': 'Gln',
    'G': 'Gly',
    'H': 'His',
    'I': 'Ile',
    'L': 'Leu',
    'K': 'Lys',
    'M': 'Met',
    'F': 'Phe',
    'P': 'Pro',
    'S': 'Ser',
    'T': 'Thr',
    'W': 'Trp',
    'Y': 'Tyr',
    'V': 'Val',
    'U': 'Sec',
    'O': 'Pyl',
    '#' : 'Mask',
    '*': 'Ter',
    '@': 'Special',
    '-' : 'Gap',
    'B': 'Asx',
    'Z': 'Glx',
    'J': 'Xle',
    'X': 'Xaa',
}
TO_AMINO_ACID_ONE_LETTER_TOKEN = {
    'Ala': 'A',
    'Arg': 'R',
    'Asn': 'N',
    'Asp': 'D',
    'Cys': 'C',
    'Glu': 'E',
    'Gln': 'Q',
    'Gly': 'G',
    'His': 'H',
    'Ile': 'I',
    'Leu': 'L',
    'Lys': 'K',
    'Met': 'M',
    'Phe': 'F',
    'Pro': 'P',
    'Ser': 'S',
    'Thr': 'T',
    'Trp': 'W',
    'Tyr': 'Y',
    'Val': 'V',
    'Sec': 'U',
    'Pyl': 'O',
    'Mask': '#',
    'Special': '@',
    'Gap': '-',
    'Ter': '*',
    'Asx': 'B',
    'Glx': 'Z',
    'Xle': 'J',
    'Xaa': 'X',
}
# Three-letter lookup tables
FROM_AMINO_ACID_THREE_LETTER_TOKEN = {
    'Ala': 'Ala',
    'Arg': 'Arg',
    'Asn': 'Asn',
    'Asp': 'Asp',
    'Cys': 'Cys',
    'Glu': 'Glu',
    'Gln': 'Gln',
    'Gly': 'Gly',
    'His': 'His',
    'Ile': 'Ile',
    'Leu': 'Leu',
    'Lys': 'Lys',
    'Met': 'Met',
    'Phe': 'Phe',
    'Pro': 'Pro',
    'Ser': 'Ser',
    'Thr': 'Thr',
    'Trp': 'Trp',
    'Tyr': 'Tyr',
    'Val': 'Val',
    'Sec': 'Sec',
    'Pyl': 'Pyl',
    'Msk': 'Mask',
    'Spl': 'Special',
    'Gap': 'Gap',
    'Ter': 'Ter',
    'Asx': 'Asx',
    'Glx': 'Glx',
    'Xle': 'Xle',
    'Xaa': 'Xaa',
}
TO_AMINO_ACID_THREE_LETTER_TOKEN = {
    'Ala': 'Ala',
    'Arg': 'Arg',
    'Asn': 'Asn',
    'Asp': 'Asp',
    'Cys': 'Cys',
    'Glu': 'Glu',
    'Gln': 'Gln',
    'Gly': 'Gly',
    'His': 'His',
    'Ile': 'Ile',
    'Leu': 'Leu',
    'Lys': 'Lys',
    'Met': 'Met',
    'Phe': 'Phe',
    'Pro': 'Pro',
    'Ser': 'Ser',
    'Thr': 'Thr',
    'Trp': 'Trp',
    'Tyr': 'Tyr',
    'Val': 'Val',
    'Sec': 'Sec',
    'Pyl': 'Pyl',
    'Mask': 'Msk',
    'Special': 'Spl',
    'Gap': 'Gap',
    'Ter': 'Ter',
    'Asx': 'Asx',
    'Glx': 'Glx',
    'Xle': 'Xle',
    'Xaa': 'Xaa',
}

# Properties
# Polar and nonpolar are disjoint
POLAR_AA = ['Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'His', 'Lys', 'Ser', 'Thr', 'Tyr']
NONPOLAR_AA = ['Ala', 'Gly', 'Ile', 'Leu', 'Met', 'Phe', 'Pro', 'Trp', 'Val']
# Acidic and basic are subsets of polar
ACIDIC_AA = ['Asp', 'Glu']
NEGATIVE_AA = ['Asp', 'Glu']
BASIC_AA = ['Arg', 'His', 'Lys']
POSITIVE_AA = ['Arg', 'His', 'Lys']
# Aromatic is a subset of nonpolar
AROMATIC_AA = ['His', 'Phe', 'Trp', 'Tyr']
# Aliphatic is a subset of nonpolar, but not aromatic
ALIPHATIC_AA = ['Ala', 'Gly', 'Ile', 'Leu', 'Pro', 'Val']
# Hydrophobic and hydrophilic are disjoint
HYDROPHOBIC_AA = ['Ala', 'Ile', 'Leu', 'Met', 'Phe', 'Pro', 'Trp', 'Val']
HYDROPHILIC_AA = ['Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'His', 'Lys', 'Ser', 'Thr', 'Tyr']
# Contains sulfur
SULFUR_AA = ['Cys', 'Met']
# Contains amide
AMIDE_AA = ['Asn', 'Gln']
# Contains hydroxyl
HYDROXYL_AA = ['Ser', 'Thr', 'Tyr']
