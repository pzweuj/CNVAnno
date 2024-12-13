# coding=utf-8
# pzw
# 20241211
# 基础函数
# 区域1完全包含区域2
def complete_contain(chrom1, start1, end1, chrom2, start2, end2):
    if chrom1 != chrom2:
        return False
    if start1 <= start2 and end1 >= end2:
        return True
    return False

# 区域1交集区域2
def overlap_intersect(chrom1, start1, end1, chrom2, start2, end2):
    if chrom1 != chrom2:
        return False
    if start1 < end2 and end1 > start2:
        return True
    return False

# 重叠程度
## 注，区域2作为数据库记录区域
def overlap_degree(chrom1, start1, end1, chrom2, start2, end2):
    if chrom1 != chrom2:
        return 0
    if complete_contain(chrom1, start1, end1, chrom2, start2, end2):
        return 1
    if overlap_intersect(chrom1, start1, end1, chrom2, start2, end2):
        coverage = (min(end1, end2) - max(start1, start2)) / (end2 - start2)
        return coverage
    return 0

# 重合相似度
def overlap_similarity(chrom1, start1, end1, chrom2, start2, end2):
    if chrom1 != chrom2:
        return 0
    if overlap_intersect(chrom1, start1, end1, chrom2, start2, end2):
        similarity = (min(end1, end2) - max(start1, start2)) / (max(end1, end2) - min(start1, start2))
        return similarity
    return 0

# end
