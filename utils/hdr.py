"""
Transformations of genome sequences for HDR.

# TODO (gdingle):
• what default MIT score threshold should we set for mutating a sequence from the original?
• what synonymous codon to use when there are multiple alternatives?
• is it okay to use MIT score on a sequence that does not end precisely in 3bp PAM as it was designed for?
"""
from typing import Iterator


class HDR:
    """
    Encapsulates all the HDR transformations of sequences described in
    https://czi.quip.com/YbAhAbOV4aXi/ . Get a mutated HDR template that varies depending on start
    or stop codon, the cut-to-insert distance, the strandedness of the guide, and the amount of
    mutation desired.

    The target sequence should be in the direction of the gene. That is,
    it has either a ATG or one of TAG, TGA, or TAA.

    The target sequence must be codon aligned!

    target_mutation_score is the minimum MIT score needed to stop silent mutation.
    """

    def __init__(
            self,
            target_seq: str,
            hdr_seq: str,
            hdr_tag: str = 'start_codon',
            hdr_dist: int = 0,
            guide_direction: str = '',
            target_mutation_score: float = 50.0) -> None:

        _validate_seq(target_seq)
        # TODO (gdingle): too strict?
        # assert len(target_seq) % 3 == 0, 'must be codon aligned'
        self.target_seq = target_seq
        _validate_seq(hdr_seq)
        self.hdr_seq = hdr_seq

        assert hdr_tag in ('start_codon', 'stop_codon')
        self.hdr_tag = hdr_tag

        assert abs(hdr_dist) < len(target_seq)
        self.hdr_dist = hdr_dist

        assert target_mutation_score < 100 and target_mutation_score > 0
        self.target_mutation_score = target_mutation_score

        if guide_direction:
            assert guide_direction in ('+', '-')
            self.guide_direction = guide_direction
            # Run inference to double check... this is failing for some!
            # # TODO (gdingle): fix me IMPORTANT!
            # self._guide_direction()
        else:
            self.guide_direction = self._guide_direction()

        if hdr_tag == 'start_codon':
            self.boundary_codons = set(['ATG'])
            # just after start codon
            self.insert_at = self._target_codon_at() + 3
        else:
            self.boundary_codons = set(['TAG', 'TGA', 'TAA'])
            # just before stop codon
            self.insert_at = self._target_codon_at()
        assert any(c in target_seq for c in self.boundary_codons)

    def __repr__(self):
        return "HDR('{}', '{}', '{}', {}, '{}', {})".format(
            self.target_seq,
            self.hdr_seq,
            self.hdr_tag,
            self.hdr_dist,
            self.guide_direction,
            self.target_mutation_score,
        )

    def _guide_direction(self) -> str:
        """
        Infer guide direction by expected PAM locations.

        We try both directions because we don't know guide direction yet.
        There is a small chance that there could be PAMs equidistant in both
        directions.

        See get_guide_cut_to_insert.

        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGC', 'NNN', hdr_dist=14)
        >>> hdr._guide_direction()
        '+'
        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGC', 'NNN', hdr_dist=1)
        >>> hdr._guide_direction()
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
        cut_at = self.insert_at + self.hdr_dist
        assert cut_at > 0
        return cut_at

    @property
    def guide_seq(self):
        """
        Returns 23bp guide sequence that includes PAM.

        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGC', 'NNN', hdr_dist=14)
        >>> hdr.guide_seq
        'ATGGCTGAGCTGGATCCGTTCGG'

        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGC', 'NNN', hdr_dist=1)
        >>> hdr.guide_seq
        'CCATGGCTGAGCTGGATCCGTTC'
        """
        cut_at = self.cut_at
        if self.guide_direction == '+':
            guide_seq = self.target_seq[cut_at - 17:cut_at + 6]
        else:
            guide_seq = self.target_seq[cut_at - 6:cut_at + 17]
        assert len(guide_seq) == 23
        return guide_seq

    def _target_codon_at(self) -> int:
        for i, codon in enumerate(_left_to_right_codons(self.target_seq)):
            if codon in self.boundary_codons:
                return i * 3
        assert False

    @property
    def template(self) -> str:
        """
        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGC', 'NNN', hdr_dist=14)
        >>> hdr.template
        'GCCATGnnnGCTGAGCTGGATCCGTTCGGC'
        """
        return self._template(False)

    @property
    def template_mutated(self) -> str:
        """
        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGC', 'NNN', hdr_dist=14)
        >>> hdr.template_mutated
        'GCCATGnnnGCTGAGCTGGATCCGtttGGC'
        """
        return self._template(True)

    def _template(self, mutate: bool = False) -> str:
        target_seq = self.mutated if mutate else self.target_seq
        return (
            target_seq[:self.insert_at] +
            self.hdr_seq.lower() +
            target_seq[self.insert_at:])

    @property
    def mutated(self) -> str:
        """
        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGC', 'NNN', hdr_dist=14)
        >>> hdr.mutated
        'GCCATGGCTGAGCTGGATCCGtttGGC'
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
        Returns 21bp subset of guide sequence aligned to codons.

        Extra base pairs are removed from the PAM side, because that is
        where we want to mutate whole codons.

        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGC', 'NNN', hdr_dist=14)
        >>> hdr.guide_seq
        'ATGGCTGAGCTGGATCCGTTCGG'
        >>> hdr.guide_seq_aligned
        'ATGGCTGAGCTGGATCCGTTC'

        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGC', 'NNN', hdr_dist=1)
        >>> hdr.guide_seq
        'CCATGGCTGAGCTGGATCCGTTC'
        >>> hdr.guide_seq_aligned
        'ATGGCTGAGCTGGATCCGTTC'

        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGG', 'NNN', hdr_dist=15)
        >>> hdr.guide_seq
        'TGGCTGAGCTGGATCCGTTCGGG'
        >>> hdr.guide_seq_aligned
        'GCTGAGCTGGATCCGTTCGGG'
        """

        # TODO (gdingle): do we want to extend to always include entire PAM?

        codon_offset = abs(self.hdr_dist % 3)
        if self.guide_direction == '+':
            aligned = self.guide_seq[:-codon_offset] if codon_offset else self.guide_seq
            return aligned[-21:]
        else:
            aligned = self.guide_seq[3 - codon_offset:] if codon_offset else self.guide_seq
            return aligned[:21]

    @property
    def guide_mutated(self) -> str:
        """
        Silently mutates codons in the guide sequence, going from the PAM side inwards.

        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGC', 'NNN', hdr_dist=14)
        >>> hdr.guide_mutated
        'ATGGCTGAGCTGGATCCGttt'

        Varying target score.
        >>> hdr.target_mutation_score = 1
        >>> hdr.guide_mutated
        'ATGGCTGAGCTGGATcccttt'
        >>> hdr.target_mutation_score = 0.1
        >>> hdr.guide_mutated
        'ATGGCTGAGCTGgaccccttt'

        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGC', 'NNN', hdr_dist=1)
        >>> hdr.guide_mutated
        'atggcgGAGCTGGATCCGTTC'
        """

        # TODO (gdingle): is it okay to use mit_hit_score on sequence that does not end precisely
        # in 3bp PAM? should we try to align to the hit_score_m? lols

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
        >>> hdr = HDR('ATGGCTGAGCTGGATCCGTTCGGC', 'NNN', 'start_codon', 14)
        >>> hdr.mutated_score
        0.05972723076923077
        """
        return mit_hit_score(
            self.guide_mutated,
            self.guide_seq_aligned,
            self.guide_direction)


def mutate_silently(guide_seq: str, guide_direction: str = '-') -> Iterator[str]:
    """
    Generator that silently mutates input sequence by substituing a different
    codon that encodes the same amino acid. Changes one codon per iteration.
    Direction is from PAM inwards.

    The input is assumed to a multiple of 3bp codons.

    Data from http://biopython.org/DIST/docs/api/Bio.SeqUtils.CodonUsage-pysrc.html

    >>> it = mutate_silently('TGTTGCGATGAC')
    >>> next(it)
    'tgcTGCGATGAC'
    >>> next(it)
    'tgctgtGATGAC'

    No possible synonyms.
    >>> next(mutate_silently('ATG'))
    'atg'

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

    if guide_direction == '+':
        codons = _right_to_left_codons(guide_seq)
    else:
        codons = _left_to_right_codons(guide_seq)

    new_guide = []
    for codon in codons:
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
    # cut_at = hdr.cut_at
    # print(hdr.insert_at, cut_at, cut_at - 6, cut_at + 17)
    # s = hdr.target_seq[cut_at - 6:cut_at + 17]
    # print(s)
