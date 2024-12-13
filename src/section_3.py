# coding=utf-8
# pzw
# 20241210

from src import common

# 处理函数
def coding_gene_intersect(row, coding_gene_df, mode='LOSS'):
    # 查询
    coding_genes = coding_gene_df[coding_gene_df.apply(lambda x: common.overlap_intersect(row["#chrom"], row["start"], row["end"], x["# chrom"], x["start"], x["end"]), axis=1)]

    # 统计数目
    refseq_coding_gene_amount = len(list(set(coding_genes["gene"].tolist())))

    # 证据计算
    mode_change = {
        "LOSS": {
            "3B": 25,
            "3C": 35
        },
        "GAIN": {
            "3B": 35,
            "3C": 50
        }
    }

    change_dict = mode_change[mode]
    if refseq_coding_gene_amount >= change_dict["3C"]:
        row["3C"] = 1
    elif refseq_coding_gene_amount >= change_dict["3B"]:
        row["3B"] = 1
    else:
        row["3A"] = 1
    return row

# 根据RefSeq Coding基因数目评价，最简单的逻辑
def section_3_cal(input_df, gene_info_df, mode='LOSS'):
    print("[Section3] Start!")
    df = input_df.copy()
    all_gene_df_no_ensembl = gene_info_df[~gene_info_df["gene"].str.startswith("ENSG")]
    coding_gene_df = all_gene_df_no_ensembl[all_gene_df_no_ensembl["gene_type"] == "protein_coding"]
    
    # 初始化
    df.loc[:, "3A"] = 0
    df.loc[:, "3B"] = 0
    df.loc[:, "3C"] = 0

    # 处理
    df = df.apply(lambda x: coding_gene_intersect(x, coding_gene_df, mode), axis=1)
    print("[Section3] Done!")
    return df

# end
