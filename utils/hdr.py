"""
Transformations of genome sequences for HDR.
"""
from typing import Iterator


# TODO (gdingle): rename HdrSeq
class HDR:
    """
    Encapsulates all the data and operations of a sequence for homolgous
    directed repair targeted at either a start codon or stop codon.

    The target sequence should be in the direction of the gene. That is,
    it have either a ATG or one of TAG, TGA, or TAA.

    target_mutation_score is the minimum MIT score needed to stop silent mutation.
    """

    def __init__(
            self,
            target_seq: str,
            hdr_seq: str,
            hdr_tag: str = 'start_codon',
            hdr_dist: int = 0,
            target_mutation_score: int = 50) -> None:

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

        assert target_mutation_score < 100 and target_mutation_score > 0
        self.target_mutation_score = target_mutation_score

    @property
    def guide_direction(self):
        """
        Get guide direction by expected PAM locations.

        # TODO (gdingle): fix ambiguity ... take guide_direction in __init__?
        We try both directions because we don't know guide direction yet.
        There is a small chance that there could be PAMs equidistant in both
        directions.

        See get_guide_cut_to_insert.

        >>> hdr = HDR('CCATGGCTGAGCTGGATCCGTTCGGC', 'NNN', hdr_dist=14)
        >>> hdr.guide_direction
        '+'
        >>> hdr = HDR('CCATGGCTGAGCTGGATCCGTTCGGC', 'NNN', hdr_dist=1)
        >>> hdr.guide_direction
        '-'
        """
        cut_at = self.cut_at
        pam1 = self.target_seq[cut_at + 3:cut_at + 6]
        pam2 = self.target_seq[cut_at - 6:cut_at - 3]
        is_for = pam1.endswith('GG')
        is_rev = pam2.startswith('CC')
        assert is_for or is_rev, (pam1, pam2)
        assert not (is_for and is_rev)
        return '+' if is_for else '-'

    @property
    def cut_at(self):
        return self.insert_at + self.hdr_dist

    @property
    def guide_seq(self):
        """
        Returns 23bp guide sequence that includes PAM.

        >>> hdr = HDR('CCATGGCTGAGCTGGATCCGTTCGGC', 'NNN', hdr_dist=14)
        >>> hdr.guide_seq
        'ATGGCTGAGCTGGATCCGTTCGG'

        >>> hdr = HDR('CCATGGCTGAGCTGGATCCGTTCGGC', 'NNN', hdr_dist=1)
        >>> hdr.guide_seq
        'CCATGGCTGAGCTGGATCCGTTC'
        """
        cut_at = self.cut_at
        if self.guide_direction == '+':
            return self.target_seq[cut_at - 17:cut_at + 6]
        else:
            return self.target_seq[cut_at - 6:cut_at + 17]

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
        """
        # TODO (gdingle): test me
        """
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

    @property
    def mutated(self) -> str:
        """
        >>> hdr = HDR('CCATGGCTGAGCTGGATCCGTTCGGC', 'NNN', hdr_dist=14)
        >>> hdr.mutated
        'CCATGGCTGAGCTGGATCCGtttGGC'
        """
        start = self.target_seq.index(self.guide_seq_aligned)
        mutated = self.guide_mutated
        return ''.join((
            self.target_seq[:start],
            mutated,
            self.target_seq[start + len(mutated):],
        ))

    @property
    def guide_seq_aligned(self) -> str:
        """
        Subset of guide sequence aligned to codons.

        >>> hdr = HDR('CCATGGCTGAGCTGGATCCGTTCGGC', 'NNN', hdr_dist=14)
        >>> hdr.guide_seq
        'ATGGCTGAGCTGGATCCGTTCGG'
        >>> hdr.guide_seq_aligned
        'ATGGCTGAGCTGGATCCGTTC'

        >>> hdr = HDR('CCATGGCTGAGCTGGATCCGTTCGGC', 'NNN', hdr_dist=1)
        >>> hdr.guide_seq
        'CCATGGCTGAGCTGGATCCGTTC'
        >>> hdr.guide_seq_aligned
        'ATGGCTGAGCTGGATCCGTTC'
        """

        # TODO (gdingle): do we want to extend to include PAM?

        if self.guide_direction == '+':
            codon_offset = abs(self.hdr_dist % 3)
            return self.guide_seq[:-codon_offset]
        else:
            codon_offset = 3 - abs(self.hdr_dist % 3)
            return self.guide_seq[codon_offset:]

    @property
    def guide_mutated(self) -> str:
        """
        Silently mutates codons in the guide sequence, going from the PAM side inwards.

        hdr_dist is relative to insert_at, which is at a codon boundary, so make
        the start of mutation a multiple of three distant.

        >>> hdr = HDR('CCATGGCTGAGCTGGATCCGTTCGGC', 'NNN', hdr_dist=14)
        >>> hdr.guide_mutated
        'ATGGCTGAGCTGGATCCGttt'
        >>> hdr.target_mutation_score = 1
        >>> hdr.guide_mutated
        'ATGGCTGAGCTGGATcccttt'

        >>> hdr = HDR('CCATGGCTGAGCTGGATCCGTTCGGC', 'NNN', hdr_dist=1)
        >>> hdr.guide_mutated
        'atggcgGAGCTGGATCCGTTC'
        """

        # TODO (gdingle): is it okay to use mit_hit_score on sequence that does not end precisely
        # in 3bp PAM? should we try to align the hit_score_m? lols

        for mutated in mutate_silently(self.guide_seq_aligned, self.guide_direction):
            score = mit_hit_score(
                mutated.upper(),
                self.guide_seq_aligned.upper(),
                self.guide_direction,
            )
            if score <= self.target_mutation_score:
                break

        return mutated

    @property
    def mutated_score(self) -> float:
        """
        >> hdr = HDR('ATGGCTGAGCTGGATCCGTTCGGC', 'NNN', 'start_codon', 14)
        >> hdr.mutated_score
        8.609700038185587e-08
        """
        # TODO (gdingle): how to get proper direction to score?
        # that includes PAM?
        return mit_hit_score(
            self.guide_mutated,
            self.guide_seq,
            self.guide_direction,)


def mutate_silently(guide_seq: str, guide_direction: str = '-') -> Iterator[str]:
    """
    Generator that silently mutates input sequence by substituing a different
    codon that encodes the same amino acid. Changes one codon per iteration.
    Direction is from PAM inwards.

    The input is assumed to start with a codon of 3bp.

    Based on https://czi.quip.com/YbAhAbOV4aXi/.

    Data from http://biopython.org/DIST/docs/api/Bio.SeqUtils.CodonUsage-pysrc.html

    >>> it = mutate_silently('TGTTGCGATGAC')
    >>> next(it)
    'tgcTGCGATGAC'
    >>> next(it)
    'tgctgtGATGAC'

    No possible synonyms.
    >>> next(mutate_silently('ATG'))
    'atg'

    Incomplete codon.
    >>> next(mutate_silently('AT'))
    'AT'
    >>> list(mutate_silently('ATGAT'))
    ['atgAT', 'atgAT']

    Right to left.
    >>> it = mutate_silently('TGTTGCGATGAC', '+')
    >>> next(it)
    'TGTTGCGATgat'
    >>> next(it)
    'TGTTGCgacgat'
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

    # TODO (gdingle):
    if guide_direction == '+':
        codons = _right_to_left_codons(guide_seq)
    else:
        codons = _left_to_right_codons(guide_seq)

    new_guide = []
    for codon in codons:

        if len(codon) < 3:
            # Exit loop on remaining bp
            new_guide.append(codon)
        else:
            # Make copy and remove current codon
            syns = list(synonymous[synonymous_index[codon]])
            syns.remove(codon)

            if len(syns):
                # TODO (gdingle): better to choose random syn?
                new_guide.append(syns.pop().lower())
            else:
                new_guide.append(codon.lower())

        if guide_direction == '+':
            new_guide_str = ''.join(new_guide[::-1])
            combined = guide_seq[:-len(new_guide_str)] + new_guide_str
        else:
            new_guide_str = ''.join(new_guide)
            combined = new_guide_str + guide_seq[len(new_guide_str):]

        assert len(combined) == len(guide_seq), (combined, guide_seq)
        yield combined


def _validate_seq(seq: str):
    assert all(b.upper() in 'AGCTN' for b in seq), seq
    assert len(seq), seq


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


def mit_hit_score(seq1: str, seq2: str, guide_direction='+') -> float:
    """Compute MIT mismatch score between two 20-mers

    See 'Scores of single hits' on http://crispr.mit.edu/about
    See calcHitScore in
    https://github.com/maximilianh/crisporWebsite/blob/master/crispor.py

    Parameters
    ----------
    seq1, seq2 : sequence
        two 20-mers to compare

    guide_direction : optional direction for starting with PAM

    Returns
    -------
    float
        MIT mismatch score between the two sequences

    Extremes.
    >>> mit_hit_score('AAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAA')
    100.0
    >>> mit_hit_score('AAAAAAAAAAAAAAAAAAAA', 'GAAAAAAAAAAAAAAAAAAA')
    100.0
    >>> mit_hit_score('AAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAG')
    41.7
    >>> mit_hit_score('ZZZZZZZZZZZZZZZZZZZZ', 'AAAAAAAAAAAAAAAAAAAA')
    8.609700038185587e-08

    Realistic.
    >>> mit_hit_score('AAGGCCAACCGGCGCCGCGC', 'GCGCGGCGCCGGTTGGCCTT')
    6.039504885480631e-06
    >>> mit_hit_score('GAAGGCCAACCGGCGCCGCG', 'CGCGGCGCCGGTTGGCCTTC')
    1.6703747039472636e-05

    Other direction.
    >>> mit_hit_score('AAAAAAAAAAAAAAAAAAAA', 'GAAAAAAAAAAAAAAAAAAA', '-')
    41.7
    >>> mit_hit_score('AAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAG', '-')
    100.0
    """
    # aka Matrix "M"
    hit_score_m = [0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445, 0.508,
                   0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583]

    # Go towards PAM
    if guide_direction == '-':
        seq1 = seq1[::-1]
        seq2 = seq2[::-1]

    assert(len(seq1) == len(seq2)), (seq1, seq2)

    # Use most important 20bp only
    seq1 = seq1[-20:]
    seq2 = seq2[-20:]

    assert(len(seq1) == 20)
    max_dist = 19

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


def _right_to_left_codons(seq: str) -> Iterator[str]:
    """
    >>> it = _right_to_left_codons('TGTTGCGATGAC')
    >>> next(it)
    'GAC'
    >>> next(it)
    'GAT'
    >>> next(it)
    'TGC'
    >>> next(it)
    'TGT'
    >>> next(it)
    Traceback (most recent call last):
    ...
    StopIteration
    """
    for i in range(len(seq), 0, -3):
        codon = seq[i - 3:i]
        yield codon


def _left_to_right_codons(seq: str) -> Iterator[str]:
    """
    >>> it = _left_to_right_codons('TGTTGCGATGAC')
    >>> next(it)
    'TGT'
    >>> next(it)
    'TGC'
    """
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        yield codon


if __name__ == '__main__':
    import doctest
    doctest.testmod()
