"""
Transformations of genome sequences for HDR.
"""


class HDR:
    """
    Encapsulates all the data and operations of a sequence for homolgous
    directed repair targeted at either a start codon or stop codon.

    The target sequence should be in the direction of the gene. That is,
    it have either a ATG or one of TAG, TGA, or TAA.
    """

    def __init__(
            self,
            target_seq: str,
            hdr_seq: str,
            hdr_tag: str = 'start_codon',
            hdr_dist: int = 0) -> None:

        _validate_seq(target_seq)
        self.target_seq = target_seq
        _validate_seq(hdr_seq)
        self.hdr_seq = hdr_seq

        assert hdr_tag in ('start_codon', 'stop_codon')
        self.hdr_tag = hdr_tag
        if hdr_tag == 'start_codon':
            self.valid_codons = set(['ATG'])
            self.insert_at = self._target_codon() + 3
        else:
            self.valid_codons = set(['TAG', 'TGA', 'TAA'])
            self.insert_at = self._target_codon()
        assert any(c in target_seq for c in self.valid_codons)

        assert hdr_dist >= 0
        self.hdr_dist = hdr_dist
        # TODO (gdingle): infer whether forward or reverse guide

    @property
    def guide_direction(self):
        """
        Based on https://czi.quip.com/YbAhAbOV4aXi/ .

        See get_guide_cut_to_insert.

        # TODO (gdingle): fix ambiguity
        We need to try both directions because hdr_dist is an absolute value.
        There is a small chance that there could be PAMs in equidistant in both
        directions.

        >>> hdr = HDR('CCATGGCTGAGCTGGATCCGTTCGGC', 'NNN', hdr_dist=14)
        >>> hdr.guide_direction
        '+'
        >>> hdr = HDR('CCATGGCTGAGCTGGATCCGTTCGGC', 'NNN', hdr_dist=1)
        >>> hdr.guide_direction
        '-'
        """
        cut_at = self.insert_at + self.hdr_dist
        pam1 = self.target_seq[cut_at + 3:cut_at + 6]
        pam2 = self.target_seq[cut_at - 6:cut_at - 3]
        is_for = pam1.endswith('GG')
        is_rev = pam2.startswith('CC')
        assert is_for or is_rev, (pam1, pam2)
        assert not (is_for and is_rev)
        return '+' if is_for else '-'

    def _target_codon(self) -> int:
        for codon in self.valid_codons:
            start = self.target_seq.find(codon)
            if start >= 0:
                return start
        assert False

    @property
    def hdr_template(self) -> str:
        return self._hdr_template(False)

    @property
    def hdr_template_mutated(self) -> str:
        return self._hdr_template(True)

    def _hdr_template(self, mutate: bool = False) -> str:
        """
        >>> HDR('ATGTCCCAGCCGGGAAT', 'NNN')._hdr_template()
        'ATGnnnTCCCAGCCGGGAAT'
        >>> HDR('TCCCAGCCGGGTGA', 'NNN', 'stop_codon')._hdr_template()
        'TCCCAGCCGGGnnnTGA'
        """
        target_seq = self.mutated if mutate else self.target_seq
        return (
            target_seq[:self.insert_at] +
            self.hdr_seq.lower() +
            target_seq[self.insert_at:])

    def _split_out_guide(self) -> tuple:
        """
        Based on example at https://czi.quip.com/YbAhAbOV4aXi/.

        >>> hdr = HDR('ATGGCTGAGCTGGATCCGTTCGGC', 'NNN', 'start_codon', 14)
        >>> hdr._split_out_guide()
        ('', 'ATGGCTGAGCTGGAT', 'CCG', 'TTCGGC', '')
        """
        assert len(self.target_seq) >= 24, 'Need to include guide plus PAM plus one for codon'
        # align to codons
        last_codon = self.hdr_dist - (self.hdr_dist % 3)
        mutate_to = self.insert_at + last_codon
        # skip cut codon
        after_cut = mutate_to + 3
        # TODO (gdingle): need to align mutate_to with codons

        # A max of 17 contiguous bp in the protospacer may survive HDR,
        # so we only touch 5 codons in that sequence. 5 * 3 == 15.
        # TODO (gdingle): stop_codons
        assert self.hdr_tag == 'start_codon', 'not implemented'
        if self.hdr_tag == 'start_codon':
            return (
                self.target_seq[:mutate_to - 15],
                self.target_seq[mutate_to - 15:mutate_to],  # part before cut
                self.target_seq[mutate_to:after_cut],
                self.target_seq[after_cut:after_cut + 6],  # part after cut
                self.target_seq[after_cut + 6:])
        else:
            return (
                self.target_seq[:mutate_to],
                self.target_seq[mutate_to:mutate_to + 15],  # part before cut
                self.target_seq[mutate_to + 15:after_cut],
                self.target_seq[after_cut:after_cut + 6],  # part after cut
                self.target_seq[after_cut:])

    @property
    def mutated(self) -> str:
        """
        >>> hdr = HDR('ATGGCTGAGCTGGATCCGTTCGGC', 'NNN', 'start_codon', 14)
        >>> hdr.mutated
        'atggcggaactagacCCGtttgga'
        """

        # TODO (gdingle): optimize mutation with mutated_score
        before, guide_left, remainder, guide_right, after = self._split_out_guide()
        return ''.join((
            before,
            mutate_silently(guide_left),
            remainder,
            mutate_silently(guide_right),
            after,
        ))

    @property
    def mutated_score(self) -> float:
        """
        >>> hdr = HDR('ATGGCTGAGCTGGATCCGTTCGGC', 'NNN', 'start_codon', 14)
        >>> hdr.mutated_score
        7.775213510406977e-06
        """
        # TODO (gdingle): do for stop_codon
        cut_at = self.insert_at + self.hdr_dist
        # guide is shifted to include PAM because we mutate the PAM
        guide_right = cut_at + 6
        guide_left = cut_at - 14

        target_guide = self.target_seq[guide_left:guide_right]
        # TODO (gdingle): make it work for reverse guide
        assert target_guide.endswith('GG')

        return mit_hit_score(
            self.mutated[guide_left:guide_right],
            target_guide)


def mutate_silently(guide_seq: str) -> str:
    """
    Silently mutates input sequence by substituing a different
    codon that encodes the same amino acid wherever possible.

    The input is assumed to start with a codon of 3bp.

    Based on https://czi.quip.com/YbAhAbOV4aXi/.

    Data from http://biopython.org/DIST/docs/api/Bio.SeqUtils.CodonUsage-pysrc.html

    >>> mutate_silently('TGTTGCGATGAC')
    'tgctgtgacgat'
    """
    synonymous = {
        # TODO (gdingle): comment out rare codons
        # see https://www.genscript.com/tools/codon-frequency-table
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
    _validate_seq(guide_seq)
    new_guide = ''
    for i in range(0, len(guide_seq), 3):
        codon = guide_seq[i:i + 3].upper()
        syns = list(synonymous[synonymous_index[codon]])
        if len(syns) > 1:
            syns.remove(codon)
        # TODO (gdingle): better to choose random syn?
        new_guide += syns.pop().lower()
    assert len(new_guide) == len(guide_seq)
    assert new_guide != guide_seq
    return new_guide


def _validate_seq(seq: str):
    assert all(b.upper() in 'AGCTN' for b in seq)
    assert len(seq)


# TODO (gdingle): we may not need this at all because we can create from hdr_template above
# def get_hdr_primer(
#         primer_product: str,
#         hdr_template: str,
#         hdr_tag: str='start_codon') -> str:
#     """
#     Locates target codon in primer product then inserts HDR sequence.
#     >>> get_hdr_primer('ATGTCCCAGCCGGGAAT', 'ATGnnn')
#     'ATGnnnTCCCAGCCGGGAAT'
#     >>> get_hdr_primer('AACAAGTGAATAAA', 'nnnTGA', 'stop_codon')
#     'AACAAGTGAAnnnTAAAAA'
#     """
#     _validate_seq(primer_product)
#     _validate_seq(hdr_template)
#     assert hdr_tag in ('start_codon', 'stop_codon')
#     if hdr_tag == 'start_codon':
#         assert hdr_template[0:3] == 'ATG'
#         codon_index = primer_product.find('ATG')
#         if codon_index == -1:
#             # TODO (gdingle): good return value?
#             return 'start_codon not found'
#         assert primer_product[codon_index:codon_index + 3] == 'ATG'
#         hdr_seq = hdr_template[3:].lower()
#         return primer_product[:codon_index] + \
#             get_hdr_template(primer_product[codon_index:], hdr_seq, hdr_tag)

#     # TODO (gdingle): this is really uglly and should be refactored.
#     elif hdr_tag == 'stop_codon':
#         stop_codons = ['TAG', 'TGA', 'TAA']
#         assert hdr_template[-3:] in stop_codons, hdr_template
#         # TODO (gdingle): does this work? waht about triplets of amino acids?
#         stops = [primer_product.rfind(stop) for stop in stop_codons]
#         if all(s == -1 for s in stops):
#             return 'stop_codon not found'
#         # TODO (gdingle): is this biologically correct?
#         codon_index = max(stops)
#         assert primer_product[codon_index:codon_index + 3] in stop_codons
#         hdr_seq = hdr_template[:-3].lower()
#         return get_hdr_template(primer_product[:codon_index + 3], hdr_seq, hdr_tag) + \
#             primer_product[-3:]

#     assert False


def mit_hit_score(seq1: str, seq2: str) -> float:
    """Compute MIT mismatch score between two 20-mers

    See 'Scores of single hits' on http://crispr.mit.edu/about
    See calcHitScore in
    https://github.com/maximilianh/crisporWebsite/blob/master/crispor.py

    Parameters
    ----------
    seq1, seq2 : sequence
        two 20-mers to compare

    Returns
    -------
    float
        MIT mismatch score between the two sequences

    >>> mit_hit_score('AAGGCCAACCGGCGCCGCGC', 'GCGCGGCGCCGGTTGGCCTT')
    6.039504885480631e-06
    >>> mit_hit_score('GAAGGCCAACCGGCGCCGCG', 'CGCGGCGCCGGTTGGCCTTC')
    1.6703747039472636e-05
    """
    # aka Matrix "M"
    hit_score_m = [0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445, 0.508,
                   0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583]

    assert(len(seq1) == 20)
    max_dist = 19

    assert(len(seq1) == len(seq2)), (seq1, seq2)

    dists = []  # distances between mismatches, for part 2
    mm_count = 0  # number of mismatches, for part 3
    last_mm_pos = None  # position of last mismatch, used to calculate distance

    score1 = 1.0
    for pos in range(0, len(seq1)):
        if seq1[pos] != seq2[pos]:
            mm_count += 1
            if last_mm_pos != None:
                dists.append(pos - last_mm_pos)  # type: ignore
            score1 *= 1 - hit_score_m[pos]
            last_mm_pos = pos
    # 2nd part of the score
    if mm_count < 2:  # special case, not shown in the paper
        score2 = 1.0
    else:
        avg_dist = sum(dists) / len(dists)
        score2 = 1.0 / (((max_dist - avg_dist) / float(max_dist)) * 4 + 1)
    # 3rd part of the score
    if mm_count == 0:  # special case, not shown in the paper
        score3 = 1.0
    else:
        score3 = 1.0 / (mm_count**2)

    return score1 * score2 * score3 * 100


if __name__ == '__main__':
    import doctest
    doctest.testmod()
