"""
Transformations of genome sequences for HDR.
"""
import itertools

from typing import Iterator

import math


class HDR:
    """
    Encapsulates all the HDR transformations of sequences described in
    https://czi.quip.com/YbAhAbOV4aXi/ . Get a mutated HDR inserted that varies
    depending on start or stop codon, the cut-to-insert distance, the
    strandedness of the guide, and the amount of mutation desired.

    The target sequence should be in the direction of the gene. Reading from
    left to right, it should have either a ATG or one of TAG, TGA, or TAA.

    The target sequence must be codon aligned so the target codon can be found!

    The CDS sequence is needed to avoid the intron/exon junctions. It must also
    contain the target codon like the target sequence.

    target_mutation_score is the minimum MIT score needed to stop silent mutation.

    guide_strand_same refers to strand of target_seq.
    """

    # When mutating, compare mutated 20-mer guide to all 20-mer sequences in
    # target_seq.
    compare_all_positions = True
    # When mutating, stop mutating at the first score that is below the
    # target_mutation_score. This may leave out some closer matches.
    stop_mutating_at_first_success = True
    # Try mutate all codons, not just linearly.
    mutate_all_permutations = False

    def __init__(
            self,
            target_seq: str,
            hdr_seq: str = '',
            hdr_tag: str = 'start_codon',
            hdr_dist: int = 0,
            guide_strand_same: bool = None,
            cds_seq: str = '',
            codon_at: int = -1,
            # Default based on analysis of
            # https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1012-2
            target_mutation_score: float = 0.1) -> None:

        _validate_seq(target_seq)
        self.target_seq = target_seq
        _validate_seq(hdr_seq)
        self.hdr_seq = hdr_seq
        _validate_seq(cds_seq)
        self.cds_seq = cds_seq

        assert hdr_tag in ('start_codon', 'stop_codon')
        self.hdr_tag = hdr_tag

        assert abs(hdr_dist) < len(target_seq)
        self.hdr_dist = hdr_dist

        assert target_mutation_score < 100 and target_mutation_score > 0
        self.target_mutation_score = target_mutation_score

        # TODO (gdingle): refactor _target_codon_at
        self._codon_at = codon_at
        if hdr_tag == 'start_codon':
            self.boundary_codons = set(['ATG'])
            # just after start codon
            self.insert_at = self._target_codon_at() + 3
        else:
            self.boundary_codons = set(['TAG', 'TGA', 'TAA'])
            # just before stop codon
            self.insert_at = self._target_codon_at()

        if guide_strand_same is not None:
            assert guide_strand_same in (True, False)
            self.guide_strand_same = guide_strand_same
            # TODO (gdingle): Run inference to double check? or remove _guide_strand_same?
            # self._guide_strand_same()
        else:
            self.guide_strand_same = self._guide_strand_same()

    def __repr__(self):
        return "HDR('{}', '{}', '{}', {}, '{}', '{}' {})".format(
            self.target_seq,
            self.hdr_seq,
            self.hdr_tag,
            self.hdr_dist,
            self.guide_strand_same,
            self.cds_seq,
            self.target_mutation_score,
        )

    def _guide_strand_same(self) -> bool:
        """
        Infer guide direction by expected PAM locations.

        We try both directions because we don't know guide direction yet.
        There is a small chance that there could be PAMs equidistant in both
        directions.

        See get_guide_cut_to_insert.

        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGC', hdr_dist=14)
        >>> hdr._guide_strand_same()
        True
        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGC', hdr_dist=1)
        >>> hdr._guide_strand_same()
        False
        """

        cut_at = self.cut_at
        pam1 = self.target_seq[cut_at + 3:cut_at + 6]
        pam2 = self.target_seq[cut_at - 6:cut_at - 3]
        is_for = pam1.endswith('GG')
        is_rev = pam2.startswith('CC')
        assert is_for or is_rev, (pam1, pam2)
        assert not (is_for and is_rev)
        return True if is_for else False

    def _target_codon_at(self) -> int:
        if self._codon_at != -1:
            return self._codon_at
        # TODO (gdingle): sometimes there is an extra stop codon that is picked up first
        # ... outside CDS? how to fix?
        for i, codon in enumerate(_left_to_right_codons(self.target_seq)):
            if codon in self.boundary_codons:
                return i * 3

        assert False

    @property
    def cut_at(self):
        cut_at = self.insert_at + self.hdr_dist
        assert cut_at >= 0, (self.insert_at, self.hdr_dist)
        return cut_at

    @property
    def junction(self) -> tuple:
        """
        Returns the intron/exon junction range on the intron side relative
        to the target sequence. Exclusive range: (start, end].

        Jason Li says: The exon-proximal side of the junction is relatively
        unimportant to the intron-proximal (exon-distal) side of each
        junction, meaning we only want to avoid the 3 nts on the intron side
        of things.

        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGC', hdr_dist=14,
        ... cds_seq='ATGGCTGAGCTGGATCC')
        >>> hdr.junction
        (20, 23)

        >>> hdr = HDR('NNNNNNTAANNNNNN', hdr_dist=0, hdr_tag='stop_codon',
        ... cds_seq='TAA', guide_strand_same=True)
        >>> hdr.junction
        (3, 6)

        >>> hdr = HDR('ATGNGG', cds_seq='ATGNNNNNN', hdr_dist=-3)
        >>> hdr.junction
        ()
        """
        index = self.target_seq.find(self.cds_seq)
        if index == -1:
            # assume target region does not go outside CDS
            return tuple([])
        if self.hdr_tag == 'start_codon':
            # Assumes junction is towards middle of gene
            index = index + len(self.cds_seq)
            return (index, index + 3)
        else:
            return (index - 3, index)

    @property
    def cut_in_junction(self) -> bool:
        """
        Determines whether the cut location is inside an intron/exon junction.

        Cut in junction.
        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGC', hdr_dist=14,
        ... cds_seq='ATGGCTGAGCTGGATCC')
        >>> (hdr.cut_at, hdr.junction, hdr.cut_in_junction)
        (20, (20, 23), True)

        Cut just after junction.
        >>> hdr = HDR('CCNNNNTAANNNNNN', hdr_dist=0, hdr_tag='stop_codon',
        ... cds_seq='TAA', guide_strand_same=True)
        >>> (hdr.cut_at, hdr.junction, hdr.cut_in_junction)
        (6, (3, 6), False)

        No junction to cut.
        >>> hdr = HDR('ATGNGG', cds_seq='ATGNNNNNN', hdr_dist=-3)
        >>> hdr.cut_in_junction
        False
        """
        junction = self.junction
        if not junction:
            return False
        return self.cut_at >= junction[0] and self.cut_at < junction[1]

    @property
    def guide_seq(self):
        """
        Returns 23bp guide sequence that includes PAM.

        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGC', hdr_dist=14)
        >>> hdr.guide_seq
        'ATGGCTGAGCTGGATCCGTTCGG'

        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGC', hdr_dist=1)
        >>> hdr.guide_seq
        'CCATGGCTGAGCTGGATCCGTTC'
        """
        cut_at = self.cut_at
        if self.guide_strand_same:
            guide_seq = self.target_seq[cut_at - 17:cut_at + 6]
        else:
            guide_seq = self.target_seq[cut_at - 6:cut_at + 17]
        assert len(guide_seq) == 23
        return guide_seq

    @property
    def guide_seq_aligned(self) -> str:
        """
        Returns subset of guide sequence aligned to codons.

        Extra base pairs are removed from the PAM side, because that is
        where we want to mutate whole codons.

        GCCATG|GCTGAGCTGGATCC|GTT|CGGC
        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGC', hdr_dist=14)
        >>> hdr.guide_seq
        'ATGGCTGAGCTGGATCCGTTCGG'
        >>> hdr.guide_seq_aligned
        'ATGGCTGAGCTGGATCCGTTC'

        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGC', hdr_dist=1)
        >>> hdr.guide_seq
        'CCATGGCTGAGCTGGATCCGTTC'
        >>> hdr.guide_seq_aligned
        'ATGGCTGAGCTGGATCCGTTC'

        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGG', hdr_dist=15)
        >>> hdr.guide_seq
        'TGGCTGAGCTGGATCCGTTCGGG'
        >>> hdr.guide_seq_aligned
        'GCTGAGCTGGATCCGTTCGGG'
        """

        # TODO (gdingle): do we want to extend to always include entire PAM?

        codon_offset = abs(self.hdr_dist % 3)
        if self.guide_strand_same == True:
            aligned = self.guide_seq[:-codon_offset] if codon_offset else self.guide_seq
            return aligned[-21:]
        else:
            aligned = self.guide_seq[3 - codon_offset:] if codon_offset else self.guide_seq
            return aligned[:21]

    @property
    def inserted(self) -> str:
        """
        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGC', 'NNN', hdr_dist=14)
        >>> hdr.inserted
        'GCCATGnnnGCTGAGCTGGATCCGTTCGGC'
        """
        return self._inserted(False)

    @property
    def inserted_mutated(self) -> str:
        """
        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGC', 'NNN', hdr_dist=14, target_mutation_score=50.0)
        >>> hdr.inserted_mutated
        'GCCATGnnnGCTGAGCTGGATCCGTTtGGC'
        """
        return self._inserted(True)

    def _inserted(self, mutate: bool = False) -> str:
        target_seq = self.mutated if mutate else self.target_seq
        return (
            target_seq[:self.insert_at] +
            self.hdr_seq.lower() +
            target_seq[self.insert_at:])

    @property
    def mutated(self) -> str:
        """
        Mutates target sequence. If the guide PAM is outside the coding region,
        the PAM is mutated in place. Otherwise, some codons in the guide are
        mutated silently.

        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGC', hdr_dist=14, target_mutation_score=50.0)
        >>> hdr.mutated
        'GCCATGGCTGAGCTGGATCCGTTtGGC'

        PAM is outside.
        >>> hdr = HDR('CCTTGGCTGATGTGGATCCGTTCGGC', hdr_dist=-12)
        >>> hdr.mutated
        'CCTTccCTGATGTGGATCCGTTCGGC'
        """
        if self.pam_outside_cds:
            # Skip other kinds of mutations because PAM mutation is enough
            return self._pam_mutated

        start = self.target_seq.index(self.guide_seq_aligned)
        mutated = self.guide_mutated
        return ''.join((
            self.target_seq[:start],
            mutated,
            self.target_seq[start + len(mutated):],
        ))

    @property
    def _pam_mutated(self) -> str:
        """
        Target seq with 3bp PAM mutated inside it.

        >>> hdr = HDR('CCTTGGCTGATGTGGATCCGTTCGGC', hdr_dist=-12)
        >>> hdr._pam_mutated
        'CCTTccCTGATGTGGATCCGTTCGGC'

        >>> hdr = HDR('ATGCCTTGGCTGATATGGATCCGT', hdr_dist=6, guide_strand_same=False)
        >>> hdr._pam_mutated
        'ATGggTTGGCTGATATGGATCCGT'
        """
        before, pam, after = (
            self.target_seq[:self.pam_at],
            self.target_seq[self.pam_at:self.pam_at + 3],
            self.target_seq[self.pam_at + 3:]
        )
        assert len(pam) == 3
        if self.guide_strand_same:
            assert 'GG' in pam, pam
            pam_mutated = pam.replace('GG', 'cc')
        else:
            assert 'CC' in pam, pam
            pam_mutated = pam.replace('CC', 'gg')
        combined = before + pam_mutated + after
        assert len(combined) == len(self.target_seq)
        return combined

    @property
    def guide_mutated(self) -> str:
        """
        Silently mutates codons in the guide sequence, going from the PAM side inwards.

        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGC', hdr_dist=14, target_mutation_score=50.0)
        >>> hdr.guide_mutated
        'ATGGCTGAGCTGGATCCGTTt'

        Varying target score.
        >>> hdr.target_mutation_score = 1
        >>> hdr.guide_mutated
        'ATGGCTGAGCTGGATCCcTTt'
        >>> hdr.target_mutation_score = 0.1
        >>> hdr.guide_mutated
        'ATGGCTGAGCTGGAcCCcTTt'
        >>> hdr.target_mutation_score = 0.01
        >>> hdr.guide_mutated
        'ATGGCcGAaCTcGAcCCcTTt'

        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGC', hdr_dist=1, target_mutation_score=50.0)
        >>> hdr.guide_mutated
        'ATGGCcGAGCTGGATCCGTTC'

        Varying target score.
        >>> hdr.target_mutation_score = 1
        >>> hdr.guide_mutated
        'ATGGCcGAaCTGGATCCGTTC'
        >>> hdr.target_mutation_score = 0.1
        >>> hdr.guide_mutated
        'ATGGCcGAaCTcGAcCCGTTC'
        >>> hdr.target_mutation_score = 0.01
        >>> hdr.guide_mutated
        'ATGGCcGAaCTcGAcCCcTTt'

        Permutations.
        >>> hdr.stop_mutating_at_first_success = False
        >>> hdr.mutate_all_permutations = True
        >>> hdr.target_mutation_score = 1
        >>> hdr.guide_mutated
        'ATGGCTGAGCTcGAcCCcTTt'
        >>> hdr.target_mutation_score = 0.1
        >>> hdr.guide_mutated
        'ATGGCcGAaCTcGAcCCGTTt'
        >>> hdr.target_mutation_score = 0.01
        >>> hdr.guide_mutated
        'ATGGCcGAaCTcGAcCCcTTt'
        """

        # TODO (gdingle): is it okay to use mit_hit_score on sequence that does not end precisely
        # in 3bp PAM? should we try to align to the hit_score_m? lols
        candidates = []
        for mutated in mutate_silently(
            self.guide_seq_aligned,
            self.guide_strand_same,
            all_permutations=self.mutate_all_permutations
        ):
            if self.compare_all_positions:
                scores = []
                assert len(self.target_seq) >= len(mutated)
                for i in range(0, len(self.target_seq) - len(mutated) + 1):
                    test_seq = self.target_seq[i:i + len(mutated)]
                    scores.append(mit_hit_score(
                        mutated.upper(),
                        test_seq.upper(),
                        self.guide_strand_same,
                    ))
                score = max(scores)
            else:
                score = mit_hit_score(
                    mutated.upper(),
                    self.guide_seq_aligned.upper(),
                    self.guide_strand_same)

            if self.stop_mutating_at_first_success:
                if score <= self.target_mutation_score:
                    return mutated
            else:
                candidates.append((score, mutated))

        if candidates:
            filtered = [(s, c) for s, c in candidates
                        if s <= self.target_mutation_score]
            if filtered:
                return max(filtered)[1]

        return mutated

    @property
    def mutated_score(self) -> float:
        """
        >>> hdr = HDR('ATGGCTGAGCTGGATCCGTTCGGC', hdr_dist=14, target_mutation_score=50.0)
        >>> hdr.guide_mutated
        'ATGGCTGAGCTGGATCCGTTt'
        >>> hdr.mutated_score
        41.7

        Verify compare_all_positions returns same for normal input.
        >>> hdr.compare_all_positions = True
        >>> hdr.mutated_score
        41.7
        >>> hdr.compare_all_positions = False
        >>> hdr.mutated_score
        41.7

        Artifical compare_all_positions example. The mutated seq was copied into target seq.
        compare_all_positions then causes more mutation.
        >>> hdr = HDR('ATGAAAAAAAAAAAAAAAAAAGG' + 'ATGAAAAAAAAAAAAAAgAAg', hdr_dist=14)
        >>> hdr.compare_all_positions = False
        >>> hdr.guide_mutated
        'ATGAAAAAAAAAAAgAAgAAg'
        >>> hdr.mutated_score
        0.06084376104417672
        >>> hdr.compare_all_positions = True
        >>> hdr.guide_mutated
        'ATGAAgAAgAAgAAgAAgAAg'
        >>> hdr.mutated_score
        0.00844207184487952
        """
        return mit_hit_score(
            self.guide_mutated,
            self.guide_seq_aligned,
            self.guide_strand_same)

    @property
    def mutation_in_junction(self) -> bool:
        """
        Determines whether there is a mutation inside an intron/exon junction.

        Mutation just inside 3 bp window.
        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGC', hdr_dist=14,
        ... cds_seq='ATGGCTGAGCTGGATCCG', target_mutation_score=50.0)
        >>> (hdr.mutated, hdr.junction, hdr.mutation_in_junction)
        ('GCCATGGCTGAGCTGGATCCGTTtGGC', (21, 24), True)

        Mutation just outside 3 bp window.
        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGC', hdr_dist=14,
        ... cds_seq='ATGGCTGAGCTGGATCC', target_mutation_score=50.0)
        >>> (hdr.mutated, hdr.junction, hdr.mutation_in_junction)
        ('GCCATGGCTGAGCTGGATCCGTTtGGC', (20, 23), False)

        >>> hdr = HDR('ATGNGG', cds_seq='ATGNNNNNN', hdr_dist=-3)
        >>> hdr.mutation_in_junction
        False
        """
        junction = self.junction
        if not junction:
            return False
        junction_seq = self.mutated[junction[0]:junction[1]]
        assert len(junction_seq) <= 3
        # Lowercase means mutated
        if any(c.lower() == c for c in junction_seq):
            return True
        else:
            return False

    @property
    def should_mutate(self) -> bool:
        """
        Determines whether a guide should be mutated depending on the cut to
        insert distance and the guide orientation. The rule is: mutate if more
        than 14bp of the PAM-side of protospacer will be intact after insertion.

        1a. 14bp or more intact on PAM side, positive guide.
        GCC|ATG|GCTGAGCTGGATCC|GTT|CGG|C
            codon              cut pam
        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGC', hdr_dist=14)
        >>> hdr.should_mutate
        True

        1b. Less than 14bp intact on PAM side, positive guide.
        GCC|ATG|GCTGAGCT|GTT|CGG|C
            codon           cut pam
        >>> hdr = HDR('GCCATGGCTGAGCTGTTCGGC', hdr_dist=8)
        >>> hdr.should_mutate
        False

        2a. 14bp or more intact on PAM side, negative guide.
        |CCA|CGA|GCGGCGGCGGCG|ATG|
         pam cut              codon
        >>> hdr = HDR('CCACGAGCGGCGGCGGCGATG', hdr_dist=-15, guide_strand_same=False)
        >>> hdr.should_mutate
        True

        2b. Less than 14bp intact on PAM side, negative guide.
        |CCA|CGA|GCG|ATG|GCTGAGCTGGATCCG
         pam cut     codon
        >>> hdr = HDR('CCACGAGCGATGGCTGAGCTGGATCCG', hdr_dist=-6, guide_strand_same=False)
        >>> hdr.should_mutate
        False

        3a. Insert is outside of guide, positive guide.
        |CCT|TGG|CTG|ATG|TGGATCCGTTCGGC
         cut pam     codon
        >>> hdr = HDR('CCTTGGCTGATGTGGATCCGTTCGGC', hdr_dist=-12)
        >>> hdr.should_mutate
        True

        3b. Insert is outside of guide, negative guide.
        |ATG|CCT|TGG|CTGATATGGATCCGT
         cod pam cut
        >>> hdr = HDR('ATGCCTTGGCTGATATGGATCCGT', hdr_dist=6, guide_strand_same=False)
        >>> hdr.should_mutate
        True
        """
        if self.guide_strand_same:
            guide_right = self.cut_at + 3
            intact = guide_right - self.insert_at
        else:
            guide_left = self.cut_at - 3
            intact = self.insert_at - guide_left

        # intact <= -3 means the insert is outside the guide + pam
        return intact <= -3 or intact >= 14

    @property
    def pam_at(self) -> int:
        """
        >>> hdr = HDR('ATGCCTTGGCTGATATGGATCCGT', hdr_dist=6, guide_strand_same=False)
        >>> hdr.pam_at
        3

        >>> hdr = HDR('TCTTGGCTGATGTGGATCCGTTCGGC', hdr_dist=-12)
        >>> hdr.pam_at
        3
        """
        if self.guide_strand_same:
            return self.cut_at + 3
        else:
            return self.cut_at - 6

    @property
    def pam_outside_cds(self) -> bool:
        """
        >>> hdr = HDR('ATGCCTTGGCTGATATGGATCCGT', hdr_dist=6, guide_strand_same=False)
        >>> hdr.pam_outside_cds
        False

        >>> hdr = HDR('CCTTGGCTGATGTGGATCCGTTCGGC', hdr_dist=-12)
        >>> hdr.pam_outside_cds
        True
        """
        if self.hdr_tag == 'start_codon':
            return self.pam_at <= self.insert_at - 6
        else:
            return self.pam_at >= self.insert_at + 6


def mutate_silently(
        guide_seq: str,
        guide_strand_same: bool = False,
        skip_stop_codon: bool = True,
        all_permutations: bool = False) -> Iterator[str]:
    """
    Generator that silently mutates input sequence by substituing a different
    codon that encodes the same amino acid. Changes one codon per iteration.
    Direction is from PAM inwards, unless all_permutations is True. The new
    codon is the selected by frequency in the human genome.

    Data from http://biopython.org/DIST/docs/api/Bio.SeqUtils.CodonUsage-pysrc.html

    The input is assumed to a multiple of 3bp codons.

    By default, does not mutate stop codons, because such mutations are not
    always silent.

    If all_permutations is True, all possible orderings of mutations will be
    returned.

    >>> it = mutate_silently('TGTTGCGATGAC')
    >>> next(it)
    'TGcTGCGATGAC'
    >>> next(it)
    'TGcTGtGATGAC'

    No possible synonyms.
    >>> next(mutate_silently('ATG'))
    'ATG'

    Right to left.
    >>> it = mutate_silently('TGTTGCGATGAC', True)
    >>> next(it)
    'TGTTGCGATGAt'
    >>> next(it)
    'TGTTGCGAcGAt'

    Skip stop codon.
    >>> it = mutate_silently('TAG')
    >>> next(it)
    'TAG'
    >>> it = mutate_silently('TAG', skip_stop_codon=False)
    >>> next(it)
    'Tga'

    all_permutations
    >>> it = mutate_silently('TGTTGCGATGAC', all_permutations=True)
    >>> next(it)
    'TGTTGCGATGAC'
    >>> next(it)
    'TGTTGCGATGAt'
    >>> next(it)
    'TGTTGCGAcGAC'

    """
    synonymous = {
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
    # Fraction of occurences among synonyms in Human genome.
    # see https://www.genscript.com/tools/codon-frequency-table
    syn_fractions = {
        'ATG': 1, 'TGG': 1, 'CAG': 0.75, 'CAC': 0.59, 'AAG': 0.58, 'GAG': 0.58,
        'TAC': 0.57, 'TTC': 0.55, 'TGC': 0.55, 'AAC': 0.54, 'GAC': 0.54, 'TGA':
        0.52, 'ATC': 0.48, 'GTG': 0.47, 'AAT': 0.46, 'GAT': 0.46, 'TTT': 0.45,
        'TGT': 0.45, 'TAT': 0.43, 'AAA': 0.42, 'GAA': 0.42, 'CTG': 0.41, 'CAT':
        0.41, 'GCC': 0.4, 'ATT': 0.36, 'ACC': 0.36, 'GGC': 0.34, 'CCC': 0.33,
        'TAA': 0.28, 'CCT': 0.28, 'ACA': 0.28, 'CCA': 0.27, 'GCT': 0.26, 'CAA':
        0.25, 'GGA': 0.25, 'GGG': 0.25, 'GTC': 0.24, 'ACT': 0.24, 'AGC': 0.24,
        'GCA': 0.23, 'TCC': 0.22, 'CGG': 0.21, 'TAG': 0.2, 'CTC': 0.2, 'AGA':
        0.2, 'AGG': 0.2, 'CGC': 0.19, 'GTT': 0.18, 'TCT': 0.18, 'ATA': 0.16,
        'GGT': 0.16, 'TCA': 0.15, 'AGT': 0.15, 'TTG': 0.13, 'CTT': 0.13, 'ACG':
        0.12, 'GTA': 0.11, 'CCG': 0.11, 'CGA': 0.11, 'GCG': 0.11, 'CGT': 0.08,
        'TTA': 0.07, 'CTA': 0.07, 'TCG': 0.06,
    }
    synonymous_index = dict(
        (codon, aa)
        for aa, codons in synonymous.items()
        for codon in codons
    )
    _validate_seq(guide_seq)

    def _mutate_codons(codons) -> Iterator:
        mutated_codons = []
        for codon in codons:
            mutated_codons.append(_mutate_codon(codon))
            yield mutated_codons

    def _mutate_codon(codon: str) -> str:
        # Make copy and remove current codon
        syns = list(synonymous[synonymous_index[codon]])
        syns.remove(codon)

        if skip_stop_codon and codon in ['TAG', 'TGA', 'TAA']:
            return codon
        elif len(syns):
            fractions = tuple((syn_fractions[syn], syn) for syn in syns)
            top = max(fractions)[1]
            lowered = ''.join([
                c if c == top[i] else top[i].lower()
                for i, c in enumerate(codon)
            ])
            return lowered
        else:
            return codon

    def _all_permutations(guide_seq) -> Iterator:
        """This will return increasing numbers of mutations, from right to left"""
        codons = list(_left_to_right_codons(guide_seq))
        masks = itertools.product([False, True], repeat=len(codons))
        for mask in masks:
            assert len(mask) == len(codons)
            new_codons = []
            for i, do_mutate in enumerate(mask):
                if do_mutate:
                    new_codons.append(_mutate_codon(codons[i]))
                else:
                    new_codons.append(codons[i])
            assert len(new_codons) == len(codons)
            yield ''.join(new_codons)

    def _pam_inwards(guide_seq) -> Iterator:
        if guide_strand_same:
            codons = _right_to_left_codons(guide_seq)
        else:
            codons = _left_to_right_codons(guide_seq)

        for mutated_codons in _mutate_codons(codons):
            if guide_strand_same:
                new_guide_str = ''.join(mutated_codons[::-1])
                combined = guide_seq[:-len(new_guide_str)] + new_guide_str
            else:
                new_guide_str = ''.join(mutated_codons)
                combined = new_guide_str + guide_seq[len(new_guide_str):]

            assert len(combined) == len(guide_seq), (combined, guide_seq)
            yield combined

    if all_permutations:
        yield from _all_permutations(guide_seq)

    else:
        yield from _pam_inwards(guide_seq)


def _validate_seq(seq: str):
    assert all(b.upper() in 'AGCTN' for b in seq), seq
    if seq != '':
        assert len(seq) >= 3, seq


def mit_hit_score(
    seq1: str,
    seq2: str,
    guide_strand_same=True,
    include_pam=False
) -> float:
    """Compute MIT mismatch score between two 20-mers or 23-mers.

    See 'Scores of single hits' on http://crispr.mit.edu/about
    See calcHitScore in
    https://github.com/maximilianh/crisporWebsite/blob/master/crispor.py

    Parameters
    ----------
    seq1, seq2 : sequence
        two 20-mers to compare

    guide_strand_same : optional direction for starting with PAM

    include_pam : optional include extra 3bp for PAM.

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
    >>> mit_hit_score('AAAAAAAAAAAAAAAAAAAA', 'GAAAAAAAAAAAAAAAAAAA', False)
    41.7
    >>> mit_hit_score('AAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAG', False)
    100.0

    Real case.
    >>> seq1 = list(reversed('CTAAGAGCATTTACACAATACA'))
    >>> seq2 = list(reversed('ctgAGAGCATTTACACAATACA'))
    >>> mit_hit_score(seq1, seq2)
    0.05972723076923077

    Include PAM.
    >>> mit_hit_score('AAAAAAAAAAAAAAAAAAAAGGG', 'AAAAAAAAAAAAAAAAAAAATGG', include_pam=True)
    100.0
    >>> mit_hit_score('AAAAAAAAAAAAAAAAAAAAAGG', 'AAAAAAAAAAAAAAAAAAAAATT', include_pam=True)
    0.20754716981132063
    """
    # aka Matrix "M"
    hit_score_m = [0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445, 0.508,
                   0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583]
    if include_pam:
        # Add some high values, determined intuitively.
        hit_score_m += [0, 0.8, 0.8]

    # Go towards PAM
    if guide_strand_same is False:
        seq1 = seq1[::-1]
        seq2 = seq2[::-1]

    assert(len(seq1) == len(seq2)), (seq1, seq2)

    # Use most important 20bp only
    if include_pam:
        seq1 = seq1[-23:]
        seq2 = seq2[-23:]
        assert (
            seq1[-3:].endswith(('GG', 'CC')) or
            seq2[-3:].endswith(('GG', 'CC'))
        ), (seq1[-3:], seq2[-3:])
        assert len(seq1) == 23
        max_dist = 22
    else:
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


def manu_score(specificity_score: float, hdr_dist: int) -> float:
    """
    Composite score to optimize guide selection by Manuel Leonetti.
    See https://goo.gl/KsUYJa.

    specificity_score as reported by Crispor.

    hdr_dist is cut-to-insert distance.

    >>> manu_score(100, 0)
    1.0
    >>> manu_score(0, 0)
    0.0
    >>> manu_score(0, 100)
    0.0
    >>> manu_score(60, 2)
    0.7232171842235433
    """
    score = _specificity_weight(specificity_score) * _dist_weight(hdr_dist)
    assert score >= 0 and score <= 1
    return score


def _specificity_weight(specificity_score: float):
    """
    >>> _specificity_weight(20)
    0
    >>> _specificity_weight(60)
    0.75
    >>> _specificity_weight(80)
    1
    """
    assert specificity_score >= 0 and specificity_score <= 100
    if specificity_score <= 45:
        return 0
    elif specificity_score >= 65:
        return 1
    else:
        return 0.05 * (specificity_score - 45)


def _dist_weight(hdr_dist: int, variance=55) -> float:
    """
    >>> _dist_weight(0)
    1.0

    >>> _dist_weight(5)
    0.7967034698934616
    >>> _dist_weight(10)
    0.402890321529133
    >>> _dist_weight(-20)
    0.026347980814448734
    """
    hdr_dist = abs(hdr_dist)  # make symmetric
    assert hdr_dist >= 0 and hdr_dist <= 100  # 100 is resonable upper bound

    # Returns a gaussian
    weight = math.exp((-1 * hdr_dist**2) / (2 * variance))
    assert weight >= 0 and weight <= 1
    return weight


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
