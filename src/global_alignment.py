from typing import List, Tuple

def build_score_matrix(
    x: str, y: str, match: int = 2, mismatch: int = -1, gap: int = -2
) -> List[List[int]]:
    # x và y normalized

    # TODO 1: tạo ma trận (m+1) x (n+1)
    # TODO 2: khởi tạo hàng 0 và cột 0 bằng gap tích lũy
    # TODO 3: điền các ô theo recurrence Needleman–Wunsch
    matrix = [[0]*(len(y)+1) for i in range(len(x)+1)]
    for i in range(1,len(x)+1):
      matrix[i][0]=matrix[i-1][0]+gap
    for j in range(1,len(y)+1):
      matrix[0][j]=matrix[0][j-1]+gap
    for i in range(1,len(x)+1):
      for j in range(1,len(y)+1):
        cost = 0
        if y[j-1]!=x[i-1]: cost = mismatch
        else: cost = match
        matrix[i][j]=max(matrix[i][j-1]+gap,matrix[i-1][j]+gap,matrix[i-1][j-1]+cost)
    return matrix

#khôi phục một optimal alignment từ score matrix.
def traceback_global(
    x: str, y: str, score_matrix: List[List[int]],
    match: int = 2, mismatch: int = -1, gap: int = -2
) -> Tuple[str, str]:
    # x và y normalized
    # TODO 1: bắt đầu tại i=len(x), j=len(y)
    # TODO 2: kiểm tra predecessor theo thứ tự DIAG > UP > LEFT
    # TODO 3: tích lũy ký tự/gap
    # TODO 4: reverse và trả về hai chuỗi aligned
    i = len(x)
    j = len(y)
    aligned_x = []
    aligned_y = []
    while(i>0 or j>0):
      if (i>0 and j>0)and(score_matrix[i][j]==score_matrix[i-1][j-1]+(mismatch if x[i-1]!=y[j-1] else match)):
        aligned_x.append(x[i-1])
        aligned_y.append(y[j-1])
        i-=1
        j-=1
      else:
        if i>0 and (j==0 or score_matrix[i][j]==score_matrix[i-1][j]+gap):
          aligned_x.append(x[i-1])
          aligned_y.append("-")
          i-=1
        else:
          aligned_y.append(y[j-1])
          aligned_x.append("-")
          j-=1
    return "".join(reversed(aligned_x)), "".join(reversed(aligned_y))

def alignment_score(
    aligned_x: str, aligned_y: str,
    match: int = 2, mismatch: int = -1, gap: int = -2
) -> int:
    # TODO: kiểm tra hai chuỗi cùng độ dài; tính score từng cột.
    # Không chấp nhận cột ('-', '-').
    if len(aligned_x)!= len(aligned_y): raise ValueError()
    score = 0
    for i in range(len(aligned_x)):
      if aligned_x[i]=="-" and aligned_y[i]=="-": raise ValueError()
      elif aligned_x[i]=="-" or aligned_y[i]=="-": score+=gap
      elif aligned_x[i]!=aligned_y[i]: score+=mismatch
      else: score+=match
    return score

def global_align(
    xx: str, yy: str, match: int = 2, mismatch: int = -1, gap: int = -2
) -> Tuple[int, str, str]:
    # xx, yy normalized
    # TODO: normalize -> build_score_matrix -> traceback_global -> return
    matrix = build_score_matrix(xx,yy,match=match,mismatch=mismatch,gap=gap)
    trace_x,trace_y = traceback_global(xx,yy,matrix,match=match,mismatch=mismatch,gap=gap)
    score = alignment_score(trace_x,trace_y,match=match,mismatch=mismatch,gap=gap)
    return score,trace_x,trace_y

def remove_gaps(seq: str) -> str:
    return seq.replace("-", "")

def check_invariants(x, y, match=2, mismatch=-1, gap=-2):
    score, ax, ay = global_align(x, y, match, mismatch, gap)
    # x,y normalized
    assert len(ax) == len(ay)
    assert remove_gaps(ax) == x
    assert remove_gaps(ay) == y
    assert alignment_score(ax, ay, match, mismatch, gap) == score
    return score, ax, ay


    