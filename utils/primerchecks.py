"""
Derived from crispr_primer.py. Provides some additional checks for primer
quality not covered by primer3.

See:
https://github.com/czbiohub/packer-images/blob/master/assets/crispr-primer/crispr_primer.py
"""

# TODO (gdingle): how are tags determined? are they always the same?
# TODO (gdingle): this should correspond to the user selected adapters, no?
LEFT_ADAPTER_TAG = 'CTCTTTCCCTACACGACGCTCTTCCGATCT'
RIGHT_ADAPTER_TAG = 'CTGGAGTTCAGACGTGTGCTCTTCCGATCT'


def is_self_binding(left_primer: str, right_primer: str) -> bool:
    """
    Do primers bind on 3' to 3' ends?
    See https://en.wikipedia.org/wiki/Primer_dimer#Mechanism_of_formation

    No binding, because no complement.
    >>> left_primer = 'GCAC'
    >>> right_primer = left_primer
    >>> is_self_binding(left_primer, right_primer)
    False

    Left binds to right.
    >>> right_primer = complementary_sequence(left_primer[::-1])
    >>> is_self_binding(left_primer, right_primer)
    True

    Left binds to left.
    >>> right_primer = left_primer + complementary_sequence(left_primer)
    >>> is_self_binding(left_primer, right_primer)
    True
    """
    last_4_left_com = complementary_sequence(left_primer[-4:])
    last_4_right_rcom = complementary_sequence(right_primer[::-1][:4])

    # left to right
    if last_4_left_com in right_primer[::-1]:
        return True
    if last_4_right_rcom in left_primer:
        return True

    # left to left and right to right
    if last_4_left_com in left_primer:
        return True
    if last_4_right_rcom in right_primer[::-1]:
        return True

    return False


def is_self_binding_with_adapters(left_primer: str, right_primer: str) -> bool:
    """
    >>> left_primer = 'ACAT'
    >>> right_primer = 'ACAT'
    >>> is_self_binding_with_adapters(left_primer, right_primer)
    False

    # left in right adapter
    >>> left_primer = 'GCAC'
    >>> right_primer = 'ACAT'
    >>> complementary_sequence(left_primer[::-1][:4]) in RIGHT_ADAPTER_TAG
    True
    >>> is_self_binding_with_adapters(left_primer, right_primer)
    True

    # true in combination with tag
    >>> left_primer = 'ACAT'
    >>> right_primer = 'TGT'
    >>> is_self_binding_with_adapters(left_primer, right_primer)
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
