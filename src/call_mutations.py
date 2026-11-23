# check sai sót ở align
def validate_alignment(refer_aln,vari_aln):
    if len(refer_aln)!=len(vari_aln): raise ValueError("2 seq phải cùng length")
    for i in range(len(refer_aln)):
        if refer_aln[i]== "-" and vari_aln[i]=="-": raise ValueError("2 seq đều có gap cùng 1 locus")

# gán nhãn cho mỗi vị trí
def label_columns(refer_aln,vari_aln):
    validate_alignment(refer_aln,vari_aln)
    list_label = []
    for i in range(len(refer_aln)):
        if refer_aln[i]==vari_aln[i]: list_label.append("match")
        elif refer_aln[i]=="-": list_label.append("insertion")
        elif vari_aln[i]=="-": list_label.append("deletion")
        else: list_label.append("substitution")
    return list_label

def call_mutations(refer_aln,vari_aln,name_vari):
    validate_alignment(refer_aln,vari_aln)
    mutations=[]
    ref_pos = 0 # số base đã qua trên refer
    current = None # lưu đột biến chèn/xóa liên tiếp đang mở để gộp thành 1 record duy nhất
    for i in range(len(refer_aln)):
        # xử lí đột biến chèn
        if refer_aln[i]=="-": 
            # trường hợp đang có chuỗi chèn mở sẵn -> gom nhóm
            if current is not None and current["type"]=="insertion":
                current["alt"]+=vari_aln[i]
                current["length"]+=1
            else: 
                # nếu trước đó có đột biến khác chưa được lưu
                if current is not None:
                    mutations.append(current)
                current = {"variant": name_vari,"position":ref_pos,"ref":"-","alt":vari_aln[i],"type":"insertion","length":1}
            continue
        ref_pos+=1
        # đột biến xóa
        if vari_aln[i]=="-":
            # trường hợp đang có chuỗi chèn mở sẵn -> gom nhóm
            if current is not None and current["type"]=="deletion":
                current["ref"]+=refer_aln[i]
                current["length"]+=1
            else: 
                # nếu trước đó có đột biến khác chưa được lưu
                if current is not None:
                    mutations.append(current)
                current = {"variant": name_vari,"position":ref_pos,"ref":refer_aln[i],"alt":"-","type":"deletion","length":1}
            continue
        # đóng cột gap đang mở 
        if current is not None:
            mutations.append(current)
            current = None
        # đột biến thay thế
        if refer_aln[i]!=vari_aln[i]:
            mutations.append({"variant": name_vari,"position":ref_pos,"ref":refer_aln[i],"alt":vari_aln[i],"type":"substitution","length":1})
    # nếu gap ở cuối đuôi
    if current is not None: mutations.append(current)
    return mutations