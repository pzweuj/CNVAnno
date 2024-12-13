# coding=utf-8
# pzw
# 20241204

"""
1A. 包含编码基因或其他已知重要元件      || 0
1B. 不包含编码基因或其他已知重要元件    || -0.6
"""

from src import common

# 处理函数
def coding_gene_intersect(row, coding_gene_df, all_gene_df_no_ensembl, important_elements_df):
    # 查询
    coding_genes = coding_gene_df[coding_gene_df.apply(lambda x: common.complete_contain(row["#chrom"], row["start"], row["end"], x["# chrom"], x["start"], x["end"]), axis=1)]
    all_genes = all_gene_df_no_ensembl[all_gene_df_no_ensembl.apply(lambda x: common.complete_contain(row["#chrom"], row["start"], row["end"], x["# chrom"], x["start"], x["end"]), axis=1)]
    important_elements = important_elements_df[important_elements_df.apply(lambda x: common.complete_contain(row["#chrom"], row["start"], row["end"], x["# chrom"], x["start"], x["end"]), axis=1)]
    
    # 形成结果，需要去重
    row["CodingGenes"] = list(set(coding_genes["gene"].tolist()))
    row["AllGenes"] = list(set(all_genes["gene"].tolist()))
    row["tmp_ImportEles"] = list(set(important_elements["gene"].tolist()))

    # 统计数目
    row["CodingGenes_Num"] = len(row["CodingGenes"])
    row["AllGenes_Num"] = len(row["AllGenes"])

    # 证据计算
    row["1A"] = row["CodingGenes_Num"] + len(row["tmp_ImportEles"])
    row["1B"] = 1 if row["1A"] == 0 else 0

    return row

# 输入df，查询是否包含编码基因或其他已知重要元件
def section_1_cal(input_df, gene_info_df):
    print("[Section1] Start!")
    df = input_df.copy()
    print("[Section1] input regions num:", len(df))

    # 重要元件列表
    important_elements_list = [
        "antisense", "lincRNA", "lncRNA", "miRNA", "misc_RNA", "Mt_rRNA", "Mt_tRNA", "ribozyme", "rRNA", "snoRNA", "snRNA", "vault_RNA",
        "IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene"
    ]

    # 拆分数据库
    # 不要ENSG基因
    all_gene_df_no_ensembl = gene_info_df[~gene_info_df["gene"].str.startswith("ENSG")]
    coding_gene_df = all_gene_df_no_ensembl[all_gene_df_no_ensembl["gene_type"] == "protein_coding"]
    important_elements_df = all_gene_df_no_ensembl[all_gene_df_no_ensembl["gene_type"].isin(important_elements_list)]
    print("[Section1] Databases Ready!")

    # 查询是否包含编码基因或其他已知重要元件
    # 给input_df初始化新列
    df.loc[:, "CodingGenes"] = df.apply(lambda row: [], axis=1)
    df.loc[:, "AllGenes"] = df.apply(lambda row: [], axis=1)
    df.loc[:, "CodingGenes_Num"] = 0
    df.loc[:, "AllGenes_Num"] = 0

    # 临时列，用来储存其他元件
    df.loc[:, "tmp_ImportEles"] = df.apply(lambda row: [], axis=1)

    # 证据列
    df.loc[:, "1A"] = 0
    df.loc[:, "1B"] = 1

    # 处理df
    df = df.apply(lambda row: coding_gene_intersect(row, coding_gene_df, all_gene_df_no_ensembl, important_elements_df), axis=1)

    # 清理临时列
    df.drop(columns=["tmp_ImportEles"], inplace=True)
    print("[Section1] Done!")

    # 输出结果
    return df

# end
