"""
Transforms genome sequences for HDR.

target_seq contains codon
target_seq contains insert loc
target_seq contains guide_seq
guide_seq contains pam
guide_seq contains cut loc
hdr_template contains hdr_seq
mutated is copy of hdr_template
mutated contains mutated guide_seq
hdr_primer contains hdr_template

codon_stop is reversed operation of codon_start

# TODO (gdingle): MIT score on mutation
"""
# For doctest
try:
    from lib.validators import validate_seq
except ImportError:
    from validators import validate_seq  # type: ignore


def get_hdr_template(target_seq: str, hdr_seq: str, hdr_tag: str='start_codon') -> str:
    """
    Inserts HDR sequence in correct position in target sequence.
    Based on https://czi.quip.com/YbAhAbOV4aXi/.

    >>> get_hdr_template('ATGTCCCAGCCGGGAAT', 'NNN')
    'ATGnnnTCCCAGCCGGGAAT'
    >>> get_hdr_template('TCCCAGCCGGGTGA', 'NNN', 'stop_codon')
    'TCCCAGCCGGGnnnTGA'
    """

    # TODO (gdingle): - [ ] need to create a Donor/hdr_template,  200bp, 90bp
    # insert (hdr_seq), 55 bp on either side of homology arm around insertion
    # point, needs to exist within primer product

    validate_seq(target_seq)
    validate_seq(hdr_seq)
    if hdr_tag == 'stop_codon':
        codon = target_seq[-3:]
        assert codon in ['TAG', 'TGA', 'TAA']
        return target_seq[0:-3] + hdr_seq.lower() + codon
    else:
        codon = target_seq[0:3]
        assert codon in ['ATG']
        return codon + hdr_seq.lower() + target_seq[3:]


def get_hdr_primer(
        primer_product: str,
        hdr_template: str,
        hdr_tag: str='start_codon') -> str:
    """
    Locates target codon in primer product then inserts HDR sequence.
    >>> get_hdr_primer('ATGTCCCAGCCGGGAAT', 'ATGnnn')
    'ATGnnnTCCCAGCCGGGAAT'
    >>> get_hdr_primer('AACAAGTGAATAAA', 'nnnTGA', 'stop_codon')
    'AACAAGTGAAnnnTAAAAA'
    """
    validate_seq(primer_product)
    validate_seq(hdr_template)
    assert hdr_tag in ('start_codon', 'stop_codon')
    if hdr_tag == 'start_codon':
        assert hdr_template[0:3] == 'ATG'
        codon_index = primer_product.find('ATG')
        if codon_index == -1:
            # TODO (gdingle): good return value?
            return 'start_codon not found'
        assert primer_product[codon_index:codon_index + 3] == 'ATG'
        hdr_seq = hdr_template[3:].lower()
        return primer_product[:codon_index] + \
            get_hdr_template(primer_product[codon_index:], hdr_seq, hdr_tag)

    # TODO (gdingle): this is really uglly and should be refactored.
    elif hdr_tag == 'stop_codon':
        stop_codons = ['TAG', 'TGA', 'TAA']
        assert hdr_template[-3:] in stop_codons, hdr_template
        # TODO (gdingle): does this work? waht about triplets of amino acids?
        stops = [primer_product.rfind(stop) for stop in stop_codons]
        if all(s == -1 for s in stops):
            return 'stop_codon not found'
        # TODO (gdingle): is this biologically correct?
        codon_index = max(stops)
        assert primer_product[codon_index:codon_index + 3] in stop_codons
        hdr_seq = hdr_template[:-3].lower()
        return get_hdr_template(primer_product[:codon_index + 3], hdr_seq, hdr_tag) + \
            primer_product[-3:]

    assert False


def mutate(row: dict, hdr_seq: str) -> str:
    """
    mutate from closest to PAM outwards, apply MIT score, halt at some threshold
    mutate including PAM
    adjust MIT score for 23bp length
    """
    hdr_template = row['hdr_template']
    assert hdr_seq in hdr_template.upper()
    # TODO (gdingle): handle guides on reverse as well
    cut = hdr_template.find(hdr_seq) + len(hdr_seq) - row['hdr_dist']
    # Align mutation region to nearest codon
    assert len(hdr_seq) % 3 == 0
    guide_right = cut - (cut % 3)
    # A max of 17 contiguous bp in the protospacer may survive HDR,
    # so we only touch 5 codons in that sequence
    assert len(row['guide_seq']) == 20
    guide_left = guide_right - 15
    mutation_target_seq = hdr_template[guide_left:guide_right]
    assert len(mutation_target_seq) == 15
    new_hdr_template = (
        hdr_template[:guide_left] +
        mutate_guide_seq(mutation_target_seq) +
        hdr_template[guide_right:]
    )
    assert len(new_hdr_template) == len(hdr_template)
    return new_hdr_template


def mutate_guide_seq(guide_seq: str) -> str:
    """
    Silently mutates input sequence by substituing a different
    codon that encodes the same amino acid wherever possible.

    The input is assumed to start with a codon of 3bp.

    Based on https://czi.quip.com/YbAhAbOV4aXi/.

    Data from http://biopython.org/DIST/docs/api/Bio.SeqUtils.CodonUsage-pysrc.html

    >>> mutate_guide_seq('TGTTGCGATGAC')
    'tgctgtgacgat'
    """
    synonymous = {
        # TODO (gdingle): comment out rare codons
        # see https://www.genscript.com/tools/codon-frequency-table
        # TODO (gdingle): consider offtarget analysis of mutated
        'CYS': ['TGT', 'TGC'],
        'ASP': ['GAT', 'GAC'],
        'SER': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'],
        'GLN': ['CAA', 'CAG'],
        'MET': ['ATG'],
        'ASN': ['AAC', 'AAT'],
        'PRO': ['CCT', 'CCG', 'CCA', 'CCC'],
        'LYS': ['AAG', 'AAA'],
        'STOP': ['TAG', 'TGA', 'TAA'],
        'THR': ['ACC', 'ACA', 'ACG', 'ACT'],
        'PHE': ['TTT', 'TTC'],
        'ALA': ['GCA', 'GCC', 'GCG', 'GCT'],
        'GLY': ['GGT', 'GGG', 'GGA', 'GGC'],
        'ILE': ['ATC', 'ATA', 'ATT'],
        'LEU': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'],
        'HIS': ['CAT', 'CAC'],
        'ARG': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'],
        'TRP': ['TGG'],
        'VAL': ['GTA', 'GTC', 'GTG', 'GTT'],
        'GLU': ['GAG', 'GAA'],
        'TYR': ['TAT', 'TAC'],
    }
    synonymous_index = dict(
        (codon, aa)
        for aa, codons in synonymous.items()
        for codon in codons
    )
    validate_seq(guide_seq)
    new_guide = ''
    for i in range(0, len(guide_seq), 3):
        codon = guide_seq[i:i + 3].upper()
        syns = set(synonymous[synonymous_index[codon]])
        if len(syns) > 1:
            syns.remove(codon)
        # TODO (gdingle): better to choose random syn?
        new_guide += syns.pop().lower()
    assert len(new_guide) == len(guide_seq)
    assert new_guide != guide_seq
    return new_guide


if __name__ == '__main__':
    import doctest
    doctest.testmod()
