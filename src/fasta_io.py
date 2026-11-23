import csv
from pathlib import Path


def read_fasta(path) -> tuple[str, str]:
    with open(path, 'r', encoding='utf-8') as file:
        line_0 = file.readline()
        if not line_0:
            raise ValueError(f'{path}: FASTA file must not be empty')
        if not line_0.startswith('>'):
            raise ValueError(f'{path}: header must start with >')

        header = line_0[1:].strip()
        if not header:
            raise ValueError(f'{path}: header must not be blank')

        the_left = file.readlines()

        if any(line.lstrip().startswith('>') for line in the_left):
            raise ValueError(f'{path}: FASTA must contain exactly one record')

        temp_the_left = [i.strip() for i in the_left]
        sequence = "".join(temp_the_left)

        if not sequence:
            raise ValueError(f'{path}: sequence must not be blank')
    return (header, sequence)

def normalize(seq) -> str:
    # Giữ ký tự ngoài ASCII để không làm lệch độ dài và vị trí lỗi (ß -> SS).
    return "".join(char.upper() if char.isascii() else char
                   for char in seq if not char.isspace())

def validate(seq: str) -> None:
    alphabet = ['A', 'T', 'G', 'C']
    if not seq:
        raise ValueError('Sequence is not exist')

    error = []
    for i in range(len(seq)):
        if seq[i]  not in alphabet:
            error.append(i)
    if len(error) != 0:
        error_char = [seq[i] for i in error]
        error_print = dict(zip(error, error_char))
        raise ValueError(f'DNA must strict DNA, detect {error_print}')

def n_invalid(seq) -> int:
    alphabet = ['A', 'T', 'G', 'C']
    count = 0
    for i in range(len(seq)):
        if seq[i]  not in alphabet:
            count += 1
    return count


def gc_content(seq) -> float:
    validate(seq)  # Sai thì hàm này tự ném ValueError.
    return 100 * (seq.count("G") + seq.count("C")) / len(seq)

def write_fasta(path, header, seq) -> None:
    with open(path, "w", encoding="utf-8") as file:
        file.write(f">{header}\n{seq}\n")



def main(project_dir: Path | None = None) -> None:
    if project_dir is None:
        project_dir = Path(__file__).resolve().parent.parent

    raw_dir = project_dir / "data" / "raw"
    clean_dir = project_dir / "data" / "clean"
    results_dir = project_dir / "results"

    input_paths = sorted(raw_dir.glob("TEM-*_cds.fasta"))
    if not input_paths:
        raise FileNotFoundError(f"Không tìm thấy file CDS trong {raw_dir}")

    # Đọc hết trước khi ghi để lỗi cấu trúc không làm thay đổi đầu ra cũ.
    records = []
    stats = []
    for input_path in input_paths:
        _, seq = read_fasta(input_path)
        seq = normalize(seq)
        if not seq:
            raise ValueError(f"{input_path}: sequence must not be blank")
        name = input_path.stem.removesuffix("_cds")
        error = None
        try:
            validate(seq)
        except ValueError as exc:
            error = str(exc)

        stats.append({
            "name": name,
            "length": len(seq),
            "gc_percent": f"{gc_content(seq):.4f}" if error is None else "",
            "n_invalid": n_invalid(seq),
        })
        records.append((name, seq, error))

    clean_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)
    for name, seq, error in records:
        output_path = clean_dir / f"{name}.fasta"
        if error is not None:
            # Không để FASTA từ lần chạy trước bị dùng nhầm cho chuỗi lỗi.
            output_path.unlink(missing_ok=True)
            print(f"{name}: skipped - {error}")
            continue
        write_fasta(output_path, name, seq)
        print(f"{name}: wrote {output_path.name}, {len(seq)} bp")

    stats_path = results_dir / "seq_stats.csv"
    with stats_path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=["name", "length", "gc_percent", "n_invalid"])
        writer.writeheader()
        writer.writerows(stats)
    print(f"Wrote {stats_path}")


if __name__ == "__main__":
    main()
