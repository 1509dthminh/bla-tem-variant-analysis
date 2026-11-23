import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from global_alignment import global_align, remove_gaps
from call_mutations import call_mutations, label_columns, validate_alignment

def mut(position, ref, alt, type_, length=1, variant="V"):
    return {"variant": variant, "position": position, "ref": ref,
            "alt": alt, "type": type_, "length": length}

def test_no_mutation():
    assert call_mutations("ATGC", "ATGC", "V") == []
    assert call_mutations("", "", "V") == []

def test_three_events_example():
    got = call_mutations("ATGCCG-TAAGCT", "ATGTCGGTA-GCT", "V")
    assert got == [mut(4, "C", "T", "substitution"),
                   mut(6, "-", "G", "insertion"),
                   mut(9, "A", "-", "deletion")]
def test_insertion_at_start_has_position_zero():
    assert call_mutations("-ATG", "GATG", "V") == [mut(0, "-", "G", "insertion")]


def test_consecutive_gaps_merged_into_one_event():
    assert call_mutations("ATGCTAA", "AT---AA", "V") == [mut(3, "GCT", "-", "deletion", 3)]
    assert call_mutations("AT--GC", "ATCCGC", "V") == [mut(2, "-", "CC", "insertion", 2)]


def test_deletion_then_insertion_are_two_events():
    got = call_mutations("ATG-C", "AT-AC", "V")
    assert [m["type"] for m in got] == ["deletion", "insertion"]


def test_invalid_alignment_raises():
    with pytest.raises(ValueError):
        call_mutations("ATG", "AT", "V")           # khác độ dài
    with pytest.raises(ValueError):
        call_mutations("A-C", "A-C", "V")          # cột '-' đối '-'
    with pytest.raises(ValueError):
        validate_alignment("A-", "A-")


def test_label_columns():
    labels = label_columns("ATGCCG-TAAGCT", "ATGTCGGTA-GCT")
    assert len(labels) == 13
    assert labels.count("substitution") == 1
    assert labels.count("insertion") == 1
    assert labels.count("deletion") == 1
    assert labels.count("match") == 10


def test_known_substitutions_recovered_through_global_align():
    rnd = random.Random(0)
    ref = "".join(rnd.choice("ACGT") for _ in range(200))
    positions = [10, 90, 199]                       # vị trí 0-based
    var = list(ref)
    for p in positions:
        var[p] = rnd.choice([c for c in "ACGT" if c != ref[p]])
    _, ax, ay = global_align(ref, "".join(var))
    got = call_mutations(ax, ay, "V")
    assert [m["position"] for m in got] == [p + 1 for p in positions]   # đổi sang 1-based
    assert all(m["type"] == "substitution" for m in got)


def test_known_deletion_recovered_through_global_align():
    # Gap tuyến tính không phân biệt gap liền hay tách rời, nên cùng một lần xóa 3 base
    # có thể được báo thành nhiều sự kiện (cùng điểm). Vì vậy kiểm tra TỔNG độ dài xóa
    # và chuỗi biến thể dựng lại được, không kiểm tra từng vị trí.
    rnd = random.Random(1)
    ref = "".join(rnd.choice("ACGT") for _ in range(150))
    var = ref[:60] + ref[63:]                       # xóa 3 base
    _, ax, ay = global_align(ref, var)
    got = call_mutations(ax, ay, "V")
    assert got and all(m["type"] == "deletion" for m in got)
    assert sum(m["length"] for m in got) == 3
    assert ay.replace("-", "") == var


def test_single_deletion_in_unique_context_is_one_event():
    # xóa 1 base trong ngữ cảnh không lặp: luôn đúng một sự kiện, đúng vị trí
    ref = "ATGCCGTAAGCT"
    _, ax, ay = global_align(ref, "ATGCCGAAGCT")
    assert len(call_mutations(ax, ay, "V")) == 1