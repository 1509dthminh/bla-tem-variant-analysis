# Phân tích biến thể gene blaTEM kháng kháng sinh

Project cuối kì MAT3390-Nhập môn Tin sinh học (kì 1 năm 2026)

# 1. Mục tiêu:
Gene *bla*TEM-1 mã hóa enzyme B-lactamase khiến vi khuẩn kháng penicillin. Các biến thể (TEM-2,TEM-3,..) khác biến thể gốc TEM-1 ở vài nucleotide. Project cố gắng trả lời các câu hỏi sau:
1. Các biến thể của TEM-1 là đột biến gì ? (thay thế, chèn hay xóa?)
2. Protein từ đột biến gene có thay đổi nhiều không ?
3. Protein từ các biến thể có cùng họ B-lactamase TEM không ?
# 2. Dữ liệu:
Source: NCBI, lưu ở data/raw/
Mỗi biến thể có 3 file: 
- '.gb': record genbank đầy đủ 
- '_full.fasta': toàn bộ record (gồm cả vùng 2 bên gene)
- '_cds.fasta': đoạn mã hóa, input của pipeline
 
# 3. Cấu trúc thư mục:
```
bla-tem-variant-analysis/
├── README.md
├── requirements.txt
├── accessions.txt
├── data/
│   ├── raw/          dữ liệu gốc
│   └── clean/        dữ liệu normalized
├── src/              
├── tests/            
├── notebooks/        
├── results/          
├── docs/             
└── slides/           
```
# 4. Pipeline:
Tiền xử lí dữ liệu -> ORF, translate -> global alignment DNA -> gọi đột biến DNA -> phân loại theo codon -> gióng hàng protein -> đo mức thay đổi của protein -> đột biến mô phỏng-> phân tích họ protein -> kiểm thử + đánh giá
# 5. Ví dụ cụ thể
Chạy tiền xử lí từ thư mục dự án (Python 3.10 trở lên):

```powershell
python src/fasta_io.py
python tests/test_fasta_io.py
```

Module FASTA và bộ test này chỉ dùng thư viện chuẩn Python. Nếu `python` chưa có
trong PATH, dùng đường dẫn tới `python.exe` của môi trường đang sử dụng.

- Đọc `data/raw/TEM-*_cds.fasta`, chuyển chữ hoa và bỏ khoảng trắng.
- Ghi từng trình tự hợp lệ vào `data/clean/TEM-n.fasta`, header `>TEM-n`.
- Ghi `results/seq_stats.csv` với các cột `name,length,gc_percent,n_invalid`.
  GC tính theo phần trăm và hiển thị bốn chữ số thập phân.
- Chuỗi có ký tự ngoài `ACGT` vẫn có dòng thống kê; GC để trống, không xuất FASTA
  và xóa file clean cùng tên từ lần chạy trước nếu có. Vị trí ký tự lỗi bắt đầu từ 0.
- File rỗng, sai cấu trúc hoặc có nhiều record làm chương trình dừng trước khi ghi đầu ra.
  Chạy lại sẽ ghi đè kết quả, không nối thêm dòng; dữ liệu raw được giữ nguyên.
# 6. Khai báo sử dụng AI
có. chi tiết thêm sau
# 7. Tài liệu tham khảo
- Giáo trình môn học
- Các paper tham khảo sẽ bổ sung sau 
