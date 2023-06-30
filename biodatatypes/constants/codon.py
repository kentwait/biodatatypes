TRANSLATION_TABLE = {
    'AAA': 'Lys',
    'AAC': 'Asn',
    'AAG': 'Lys',
    'AAT': 'Asn',
    'ACA': 'Thr',
    'ACC': 'Thr',
    'ACG': 'Thr',
    'ACT': 'Thr',
    'AGA': 'Arg',
    'AGC': 'Ser',
    'AGG': 'Arg',
    'AGT': 'Ser',
    'ATA': 'Ile',
    'ATC': 'Ile',
    'ATG': 'Met',
    'ATT': 'Ile',
    'CAA': 'Gln',
    'CAC': 'His',
    'CAG': 'Gln',
    'CAT': 'His',
    'CCA': 'Pro',
    'CCC': 'Pro',
    'CCG': 'Pro',
    'CCT': 'Pro',
    'CGA': 'Arg',
    'CGC': 'Arg',
    'CGG': 'Arg',
    'CGT': 'Arg',
    'CTA': 'Leu',
    'CTC': 'Leu',
    'CTG': 'Leu',
    'CTT': 'Leu',
    'GAA': 'Glu',
    'GAC': 'Asp',
    'GAG': 'Glu',
    'GAT': 'Asp',
    'GCA': 'Ala',
    'GCC': 'Ala',
    'GCG': 'Ala',
    'GCT': 'Ala',
    'GGA': 'Gly',
    'GGC': 'Gly',
    'GGG': 'Gly',
    'GGT': 'Gly',
    'GTA': 'Val',
    'GTC': 'Val',
    'GTG': 'Val',
    'GTT': 'Val',
    'TAA': 'Ter',
    'TAC': 'Tyr',
    'TAG': 'Ter',
    'TAT': 'Tyr',
    'TCA': 'Ser',
    'TCC': 'Ser',
    'TCG': 'Ser',
    'TCT': 'Ser',
    'TGA': 'Ter',
    'TGC': 'Cys',
    'TGG': 'Trp',
    'TGT': 'Cys',
    'TTA': 'Leu',
    'TTC': 'Phe',
    'TTG': 'Leu',
    'TTT': 'Phe',
    # Special triplet tokens
    'Mask': 'Mask',
    'Gap': 'Gap',
    'Special': 'Special',
    # Ambiguous triplet tokens
    'NNN': 'Xaa'   
}
FROM_TRIPLET_TOKEN = {
    'AAA': 'AAA',
    'AAC': 'AAC',
    'AAG': 'AAG',
    'AAT': 'AAT',
    'ACA': 'ACA',
    'ACC': 'ACC',
    'ACG': 'ACG',
    'ACT': 'ACT',
    'AGA': 'AGA',
    'AGC': 'AGC',
    'AGG': 'AGG',
    'AGT': 'AGT',
    'ATA': 'ATA',
    'ATC': 'ATC',
    'ATG': 'ATG',
    'ATT': 'ATT',
    'CAA': 'CAA',
    'CAC': 'CAC',
    'CAG': 'CAG',
    'CAT': 'CAT',
    'CCA': 'CCA',
    'CCC': 'CCC',
    'CCG': 'CCG',
    'CCT': 'CCT',
    'CGA': 'CGA',
    'CGC': 'CGC',
    'CGG': 'CGG',
    'CGT': 'CGT',
    'CTA': 'CTA',
    'CTC': 'CTC',
    'CTG': 'CTG',
    'CTT': 'CTT',
    'GAA': 'GAA',
    'GAC': 'GAC',
    'GAG': 'GAG',
    'GAT': 'GAT',
    'GCA': 'GCA',
    'GCC': 'GCC',
    'GCG': 'GCG',
    'GCT': 'GCT',
    'GGA': 'GGA',
    'GGC': 'GGC',
    'GGG': 'GGG',
    'GGT': 'GGT',
    'GTA': 'GTA',
    'GTC': 'GTC',
    'GTG': 'GTG',
    'GTT': 'GTT',
    'TAA': 'TAA',
    'TAC': 'TAC',
    'TAG': 'TAG',
    'TAT': 'TAT',
    'TCA': 'TCA',
    'TCC': 'TCC',
    'TCG': 'TCG',
    'TCT': 'TCT',
    'TGA': 'TGA',
    'TGC': 'TGC',
    'TGG': 'TGG',
    'TGT': 'TGT',
    'TTA': 'TTA',
    'TTC': 'TTC',
    'TTG': 'TTG',
    'TTT': 'TTT',
    '---': 'Gap',
    '###': 'Mask',
    '@@@': 'Special',
    'NNN': 'NNN',
}
STOP_CODONS = ['TAA', 'TAG', 'TGA']

# Degenerate codons
# Unique codons
ONEFOLD_DEGENERATE_CODONS = [
    'ATG',  # Met
    'TGG',  # Trp
]
# 2 codons code for the same amino acid
TWOFOLD_DEGENERATE_CODONS = [
    'AAT', 'AAC',  # Asn
    'GAT', 'GAC',  # Asp
    'TGT', 'TGC',  # Cys
    'CAA', 'CAG',  # Gln
    'GAA', 'GAG',  # Glu
    'CAT', 'CAC',  # His
    'AAA', 'AAG',  # Lys
    'TTT', 'TTC',  # Phe
    'TAT', 'TAC',  # Tyr
]
# 3 codons code for the same amino acid
THREEFOLD_DEGENERATE_CODONS = [
    'ATT', 'ATC', 'ATA',  # Ile
]
# 4 codons code for the same amino acid
FOURFOLD_DEGENERATE_CODONS = [
    'GCT', 'GCC', 'GCA', 'GCG',  # Ala
    'GGT', 'GGC', 'GGA', 'GGG',  # Gly
    'CCT', 'CCC', 'CCA', 'CCG',  # Pro
    'ACT', 'ACC', 'ACA', 'ACG',  # Thr
    'GTT', 'GTC', 'GTA', 'GTG',  # Val
]
# 6 codons code for the same amino acid
SIXFOLD_DEGENERATE_CODONS = [
    'CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG',  # Leu
    'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG',  # Arg
    'TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC',  # Ser
]
