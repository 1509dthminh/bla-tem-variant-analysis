# Phân tích biến thể gene blaTEM kháng kháng sinh

Project cuối kì MAT3390-Nhập môn Tin sinh học (kì 1 năm 2026)

# 1. Mục tiêu:
Gene *bla*TEM-1 mã hóa enzyme B-lactamase khiến vi khuẩn kháng penicillin. Các biến thể (TEM-2,TEM-3,..) khác biến thể gốc TEM-1 ở vài nucleotide. Project sẽ:
1. Normalize trình tự từng biến thể 
2. ORF và translate ra protein
3. Global alignment (Needleman-Wunsch) giữa TEM-1 và từng biến thể
4. Xác định đột biến 
5. Xác định protein bằng BLASTp
6. Kiểm định kết quả

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
thêm sau
# 5. Ví dụ cụ thể
thêm sau
# 6. Khai báo sử dụng AI
có. chi tiết thêm sau
# 7. Tài liệu tham khảo
- Giáo trình môn học
- Các paper tham khảo sẽ bổ sung sau 
