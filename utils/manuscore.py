import math

# TODO (gdingle): re-enable after calibration
# from functools import lru_cache

_specificity_weight_low = 45
_specificity_weight_high = 65
_dist_weight_variance = 55
_before_start_codon_penalty = 0.2


# TODO (gdingle): re-enable after calibration
# @lru_cache(maxsize=1024 * 1024)
def manu_score(specificity_score: float, hdr_dist: int, codon='stop_codon') -> float:
    """
    Composite score to optimize guide selection by Manuel Leonetti.
    See https://goo.gl/KsUYJa.

    specificity_score as reported by Crispor.

    hdr_dist is cut-to-insert distance.

    codon was added to avoid regulatory region before or inside of start codon.

    >>> manu_score(100, 0)
    1.0
    >>> manu_score(0, 0)
    0.0
    >>> manu_score(0, 100)
    0.0
    >>> manu_score(60, 2)
    0.7232171842235433
    >>> manu_score(100, 0, 'start_codon')
    1.0
    >>> manu_score(100, -1, 'start_codon')
    0.19819005765761588
    >>> manu_score(60, -2, 'start_codon')
    0.14464343684470868
    """
    assert codon in ('start_codon', 'stop_codon'), codon
    score = _specificity_weight(specificity_score) * _dist_weight(hdr_dist)
    if codon == 'start_codon' and hdr_dist < -1:
        score *= _before_start_codon_penalty
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
    low = _specificity_weight_low
    high = _specificity_weight_high

    assert specificity_score >= 0 and specificity_score <= 100
    if specificity_score <= low:
        return 0
    elif specificity_score >= high:
        return 1
    else:
        return 1 / (high - low) * (specificity_score - low)


def _dist_weight(hdr_dist: int) -> float:
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
    variance = _dist_weight_variance

    hdr_dist = abs(hdr_dist)  # make symmetric
    assert hdr_dist >= 0 and hdr_dist <= 100  # 100 is resonable upper bound

    # Returns a gaussian
    weight = math.exp((-1 * hdr_dist**2) / (2 * variance))
    assert weight >= 0 and weight <= 1
    return weight


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.FAIL_FAST | doctest.ELLIPSIS)
