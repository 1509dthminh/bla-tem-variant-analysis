import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from global_alignment import (alignment_score, build_score_matrix, check_invariants,
                   global_align, remove_gaps)

def rand_dna(rnd, n):
    return "".join(rnd.choice("ACGT") for _ in range(n))

def test_edges():
    assert check_invariants("", "")[0] == 0
    assert check_invariants("A", "")[0] == -2
    assert check_invariants("", "A")[0] == -2
def test_substitutations(): # đột biến thay thế
    score, ax, ay = check_invariants("ACGT","AAGT")
    assert score ==5
    assert "-" not in ax and "-" not in ay
def test_deletion(): 
    score, ax, ay=check_invariants("ATGCCGTAAGCT", "ATGCCGAAGCT")
    assert score == 20                
    assert ax.count("-") == 0 and ay.count("-") == 1
def test_insertion():
    score, ax, ay = check_invariants("ATGCCGAAGCT", "ATGCCGTAAGCT")
    assert score == 20
    assert ax.count("-") == 1 and ay.count("-") == 0

def test_hand_matrix_2x2():
    # x="AC", y="AG"; match=2, mismatch=-1, gap=-2
    assert build_score_matrix("AC", "AG") == [
        [0, -2, -4],
        [-2, 2, 0],
        [-4, 0, 1],
    ]


def test_matrix_borders_cumulative_gap():
    m = build_score_matrix("ACGT", "AC")
    assert [row[0] for row in m] == [0, -2, -4, -6, -8]
    assert m[0] == [0, -2, -4]

def test_matrix_shape():
    m = build_score_matrix("ACGT", "AC")
    assert len(m) == 5 and all(len(r) == 3 for r in m)

def test_symmetry():
    for a, b in [("ATGCCGTAAGCT", "ATGCCGAAGCT"), ("ACGT", "AGT"), ("", "ACG")]:
        assert global_align(a, b)[0] == global_align(b, a)[0]


def test_random_pairs_invariants():
    rnd = random.Random(0)
    for _ in range(300):
        a = rand_dna(rnd, rnd.randint(0, 30))
        b = rand_dna(rnd, rnd.randint(0, 30))
        score, _, _ = check_invariants(a, b)
        assert score == build_score_matrix(a, b)[len(a)][len(b)]   # = ô cuối ma trận
        assert score == global_align(b, a)[0]                       # đối xứng


def test_score_upper_bound():
    rnd = random.Random(1)
    for _ in range(100):
        a = rand_dna(rnd, rnd.randint(1, 25))
        b = rand_dna(rnd, rnd.randint(1, 25))
        assert global_align(a, b)[0] <= 2 * min(len(a), len(b))

def test_deterministic():
    a, b = "ATGCCGTAAGCT", "ATGCCGAAGCT"
    assert global_align(a, b) == global_align(a, b)

def test_tie_break_diag_first():
    # "AAAA" vs "AA" có nhiều alignment cùng điểm 0; quy ước DIAG > UP > LEFT
    score, ax, ay = global_align("AAAA", "AA")
    assert score == 0
    assert ax == "AAAA"
    assert ay == "--AA"

def test_known_substitutions_recovered(k):
    rnd = random.Random(100 + k)
    ref = rand_dna(rnd, 300)
    positions = rnd.sample(range(300), k)
    var = list(ref)
    for p in positions:
        var[p] = rnd.choice([c for c in "ACGT" if c != ref[p]])
    var = "".join(var)
    score, ax, ay = check_invariants(ref, var)
    assert count_gaps(ax, ay) == 0
    assert count_mismatches(ax, ay) == k
    assert score == 2 * (300 - k) - k


def test_known_deletion_position_recovered():
    rnd = random.Random(7)
    ref = rand_dna(rnd, 200)
    var = ref[:100] + ref[101:]          # xóa base thứ 101 (chỉ số 100)
    score, ax, ay = check_invariants(ref, var)
    assert score == 2 * 199 - 2
    assert ax.count("-") == 0 and ay.count("-") == 1
    assert remove_gaps(ay) == var

def test_codon_indel_keeps_frame():
    # xóa liền 3 nt (một codon): 3 gap liên tiếp, phần còn lại khớp hết
    rnd = random.Random(5)
    ref = rand_dna(rnd, 120)
    var = ref[:60] + ref[63:]
    score, ax, ay = check_invariants(ref, var)
    assert score == 2 * 117 - 2 * 3
    assert count_gaps(ax, ay) == 3