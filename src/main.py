import csv
from pathlib import Path
from fasta_io import read_fasta, normalize
from global_alignment import global_align
from call_mutations import call_mutations

def generate_mutations_dna_csv(project_dir: Path | None = None) -> Path:
    if project_dir is None:
        project_dir = Path(__file__).resolve().parent.parent

    clean_dir = project_dir / "data" / "clean"
    results_dir = project_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    # 1. Đọc chuỗi tham chiếu gốc TEM-1
    ref_path = clean_dir / "TEM-1.fasta"
    if not ref_path.exists():
        raise FileNotFoundError(f"Không tìm thấy file gene gốc: {ref_path}")
    _, ref_seq = read_fasta(ref_path)
    ref_seq = normalize(ref_seq)

    # 2. Tìm tất cả các biến thể còn lại (TEM-2, TEM-3, TEM-4,...)
    variant_files = sorted(clean_dir.glob("TEM-*.fasta"))
    all_mutations = []

    for var_path in variant_files:
        var_name = var_path.stem
        if var_name == "TEM-1":
            continue  # Bỏ qua gene gốc

        _, var_seq = read_fasta(var_path)
        var_seq = normalize(var_seq)

        # 3. Gióng hàng toàn cục giữa TEM-1 và biến thể
        _, aln_ref, aln_var = global_align(ref_seq, var_seq)

        # 4. Gọi đột biến
        muts = call_mutations(aln_ref, aln_var, var_name)
        all_mutations.extend(muts)

    # 5. Ghi kết quả ra results/mutations_dna.csv
    output_path = results_dir / "mutations_dna.csv"
    fieldnames = ["variant", "position", "ref", "alt", "type", "length"]

    with output_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_mutations)

    print(f"Đã lưu kết quả đột biến DNA vào: {output_path}")
    return output_path

if __name__ == "__main__":
    generate_mutations_dna_csv()