"""
Derived from crispr_primer.py. Provides some additional checks for primer
quality not covered by primer3.

See:
https://github.com/czbiohub/packer-images/blob/master/assets/crispr-primer/crispr_primer.py
"""

# TODO (gdingle): how are tags determined? are they always the same?
LEFT_ADAPTER_TAG = 'CTCTTTCCCTACACGACGCTCTTCCGATCT'
RIGHT_ADAPTER_TAG = 'CTGGAGTTCAGACGTGTGCTCTTCCGATCT'


def is_left_right_binding(left_primer: str, right_primer: str) -> bool:
    """
    >>> left_primer = 'GCAC'
    >>> right_primer = 'ACAT'
    >>> is_left_right_binding(left_primer, right_primer)
    False
    """
    last_4_left_com = complementary_sequence(left_primer[-4:])
    last_4_right_rcom = complementary_sequence(right_primer[::-1][:4])
    if last_4_left_com in right_primer[::-1]:
        return True
    if last_4_right_rcom in left_primer:
        return True
    return False


def is_same_side_binding(left_primer: str, right_primer: str) -> bool:
    """
    >>> left_primer = 'GCAC'
    >>> right_primer = 'ACAT'
    >>> is_same_side_binding(left_primer, right_primer)
    False
    """
    last_4_left_com = complementary_sequence(left_primer[-4:])
    last_4_right_rcom = complementary_sequence(right_primer[::-1][:4])
    if last_4_left_com in left_primer:
        return True
    if last_4_right_rcom in right_primer[::-1]:
        return True
    return False


def is_self_binding_with_tags(left_primer: str, right_primer: str) -> bool:
    """
    >>> left_primer = 'GCAC'
    >>> right_primer = 'ACAT'
    >>> is_self_binding_with_tags(left_primer, right_primer)
    True
    """
    last_4_left_com = complementary_sequence(left_primer[-4:])
    last_4_right_rcom = complementary_sequence(right_primer[::-1][:4])
    if last_4_left_com in (RIGHT_ADAPTER_TAG + right_primer)[::-1]:
        return True
    if last_4_right_rcom in (LEFT_ADAPTER_TAG + left_primer):
        return True
    return False


def complementary_sequence(seq: str) -> str:
    """
    >>> complementary_sequence('ATCG')
    'TAGC'
    """
    # TODO (gdingle): refactor with fastqs reverse_complement
    seq_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join([seq_map[c] for c in seq.upper()])


if __name__ == '__main__':
    import doctest
    doctest.testmod()
