from pathlib import Path 
from Bio import SeqIO

RAW_DATA = Path("data/raw")
STOPS = {"TAA","TAG","TGA"}
BACTERIAL_GENETIC_CODE = 11  # bảng mã di truyền chuẩn cho vi khuẩn (NCBI table 11)

genebank = sorted(RAW_DATA.glob("TEM-?.gb"))
if not genebank:
    print(f"Warning: Không thấy file có dạng TEM-?.gb trong{RAW_DATA}")

# check chất lg file
for gb in genebank:
    name = gb.stem # loại bỏ phần đuôi mở rộng
    fasta_filepath = RAW_DATA / f"{name}_cds.fasta"

    #check xem file fasta đang xét có tồn tại không
    if not fasta_filepath.exists():
        print(f"Warning: không thấy file này")
        continue
    # đọc dữ liệu từ 2  định dạng file
    gb_rec = SeqIO.read(gb,"genbank")
    fasta_rec = SeqIO.read(fasta_filepath,"fasta")
    # tìm chú thích phần mã hóa (cds) trong file gb
    cds_features = [feat for feat in gb_rec.features if feat.type == "CDS"]
    if not cds_features:
        print(f"Warning: bỏ qua file {name} do không có vùng cds")
        continue

    primary_cds = cds_features[0]

    # trích vùng dna
    gb_cds_seq = primary_cds.location.extract(gb_rec.seq)
    gb_cds_dna = str(gb_cds_seq).upper()
    fasta_cds_dna = str(fasta_rec.seq).upper()

    # trích + translate protein
    # trích trong gb
    gb_annotated_protein = primary_cds.qualifiers["translation"][0]
    # translate từ dna vừa trích sang protein
    translate_protein = str(gb_cds_seq.translate(
        table=BACTERIAL_GENETIC_CODE, 
        to_stop=True  #dừng trans khi gặp stop, không print '*'
    ))

    # check logic
    dna_matching = (fasta_cds_dna==gb_cds_dna)
    valid_start = (fasta_cds_dna.startswith("ATG"))
    valid_stop = (fasta_cds_dna[-3:] in STOPS)
    divisible_by3 = (len(fasta_cds_dna)%3==0)
    protein_matching = (translate_protein==gb_annotated_protein)

    # print result 
    # chuyển tọa độ python sang tọa độ sinh học
    bio_start_pos = int(primary_cds.location.start)+1
    bio_end_pos = int(primary_cds.location.end)
    strand_direction = "+" if primary_cds.location.strand == 1 else "-"

    print(f"Biến thể: {name}")
    print(f"Tổng quan: chiều dài {len(gb_rec)} bp | "
        f"CDS: {bio_start_pos}..{bio_end_pos} (Mạch {strand_direction}) | "
        f"độ dài CDS: {len(fasta_cds_dna)} bp")
    print("-"*65)
    print(f"DNA file FASTA khớp với GenBank CDS: {dna_matching}")
    print(f"Bắt đầu bằng codon mở đầu (ATG): {valid_start}")
    print(f"Kết thúc bằng codon dừng hợp lệ: {valid_stop}")
    print(f"Độ dài chia hết cho 3 (đúng khung): {divisible_by3}")
    print(f"Protein tự dịch khớp với GenBank: {protein_matching}")