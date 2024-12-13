# coding=utf-8
# pzw
# 20241204

from src import common
import pandas as pd

# 通用分析
def section_2_2a(row, gene_hi_df, region_hi_df):
    complete_overlap_gene = gene_hi_df[gene_hi_df.apply(lambda x: common.overlap_similarity(row["#chrom"], row["start"], row["end"], x["# chrom"], x["start"], x["end"]) > 0.75, axis=1)]
    complete_overlap_region = region_hi_df[region_hi_df.apply(lambda x: common.overlap_similarity(row["#chrom"], row["start"], row["end"], x["# chrom"], x["start"], x["end"]) > 0.75, axis=1)]
    # 2A 完全覆盖已建立的HI区
    row["2A"] = len(complete_overlap_gene) + len(complete_overlap_region)
    return row

def section_2_2b(row, gene_hi_vus_df, region_hi_vus_df):
    # 首先要部分重叠
    partial_overlap_gene = gene_hi_vus_df[gene_hi_vus_df.apply(lambda x: common.overlap_intersect(row["#chrom"], row["start"], row["end"], x["# chrom"], x["start"], x["end"]), axis=1)]
    partial_overlap_region = region_hi_vus_df[region_hi_vus_df.apply(lambda x: common.overlap_intersect(row["#chrom"], row["start"], row["end"], x["# chrom"], x["start"], x["end"]), axis=1)]

    # 不用纠结是否与已知致病区的重叠，因为2A已经处理了，在最后矛盾处理优先2A即可
    row["2B"] = len(partial_overlap_gene) + len(partial_overlap_region)
    return row

# LOSS模式分析
def section_2_loss_other_evi(row, cds_hi_df, utr_hi_df, exon_hi_df):
    # 2C分析，首先确认是否与第一个外显子重叠
    first_exon_df = exon_hi_df[(exon_hi_df["First_Exon"] == "Y")]
    first_exon_overlap = first_exon_df[first_exon_df.apply(lambda x: common.overlap_intersect(row["#chrom"], row["start"], row["end"], x["# chrom"], x["start"], x["end"]), axis=1)]
    current_gene_list_2c = first_exon_overlap["gene"].tolist()

    ## 2C-1分析，是否包含CDS
    first_cds_df = cds_hi_df[(cds_hi_df["First_Exon"] == "Y") & (cds_hi_df["gene"].isin(current_gene_list_2c))]
    partial_overlap_first_cds = first_cds_df[first_cds_df.apply(lambda x: common.overlap_intersect(row["#chrom"], row["start"], row["end"], x["# chrom"], x["start"], x["end"]), axis=1)]
    row["2C-1"] = len(partial_overlap_first_cds)

    ## 2C-2分析，是否只与UTR5重叠，即分析它不与CDS重叠的区域即可
    utr5_df = utr_hi_df[(utr_hi_df["location"] == "UTR5") & (utr_hi_df["gene"].isin(current_gene_list_2c))]
    partial_overlap_utr5 = utr5_df[utr5_df.apply(lambda x: common.overlap_intersect(row["#chrom"], row["start"], row["end"], x["# chrom"], x["start"], x["end"]), axis=1)]
    row["2C-2"] = 0 if row["2C-1"] > 0 else len(partial_overlap_utr5)

    # 2D分析，首先确认是否与最后一个外显子重叠
    last_exon_df = exon_hi_df[(exon_hi_df["Last_Exon"] == "Y")]
    last_exon_overlap = last_exon_df[last_exon_df.apply(lambda x: common.overlap_intersect(row["#chrom"], row["start"], row["end"], x["# chrom"], x["start"], x["end"]), axis=1)]
    current_gene_list_2d = last_exon_overlap["gene"].tolist()

    ## 是否与多个CDS重叠
    partial_overlap_cds = cds_hi_df[(cds_hi_df["gene"].isin(current_gene_list_2d))]
    ### 统计基因覆盖CDS个数
    partial_overlap_cds_gene_count = {}
    for gene in current_gene_list_2d:
        partial_overlap_cds_gene_count[gene] = len(partial_overlap_cds[partial_overlap_cds["gene"] == gene])
    
    # 任意基因大于1，修正2D-4
    if any(value > 1 for value in partial_overlap_cds_gene_count.values()):
        row["2D-4"] = 1
    elif all(value == 0 for value in partial_overlap_cds_gene_count.values()):
        row["2D-1"] = 1
    elif all(value == 1 for value in partial_overlap_cds_gene_count.values()):
        # 此时再判断是否有致病位点检出
        last_exon_overlap_clinvar = last_exon_overlap[last_exon_overlap["clinvar_count"] > 0]
        if len(last_exon_overlap_clinvar) > 0:
            row["2D-2"] = 1
        else:
            row["2D-3"] = 1
    return row

# 2E分析 || 需要判断PVS1
## ClassifyCNV的方案是采用一个现成的PVS1列表进行交集判断
## AutoCNV的方案是使用AutoPVS1来重新判断
## 我既不想部署AutoPVS1，也不想使用VEP注释，还不想用ClassifyCNV的数据库，因此直接判断基因
def section_2_loss_2e(row, gene_info_df, gene_moi_dict):
    # 2E分析，区域被一个基因完全覆盖
    overlap_df_2e = gene_info_df[gene_info_df.apply(lambda x: common.complete_contain(x["# chrom"], x["start"], x["end"], row["#chrom"], row["start"], row["end"]), axis=1)]
    pvs1 = "UK"
    if len(overlap_df_2e) == 1:
        # 判断PVS1逻辑, 首先去获得基因
        gene = overlap_df_2e['gene'].tolist()[0]
        moi = gene_moi_dict.get(gene, "UK")
        if moi == "AD":
            pvs1 = "PVS1"
        elif moi == "AR":
            pvs1 = "PVS1_strong"
        elif moi == "XL":
            pvs1 = "PVS1_Moderate"
        else:
            pvs1 = "PVS1_Supporting"
    row["2E"] = pvs1
    return row

# 2F分析，被良性区域完全覆盖
def section_2_loss_2f(row, region_hi_b_df):
    overlap_df_2f = region_hi_b_df[region_hi_b_df.apply(lambda x: common.complete_contain(x["# chrom"], x["start"], x["end"], row["#chrom"], row["start"], row["end"]), axis=1)]
    if len(overlap_df_2f) == 1:
        row["2F"] = 1
    return row

# 2G分析
def section_2_loss_2g(row, region_hi_b_df):
    # 2G分析，与良性区域重叠度大于75%
    overlap_df_2g = region_hi_b_df[region_hi_b_df.apply(lambda x: common.overlap_similarity(x["# chrom"], x["start"], x["end"], row["#chrom"], row["start"], row["end"]) > 0.75, axis=1)]
    if len(overlap_df_2g) > 1 and row["2F"] == 0:
        row["2G"] = 1
    return row

# 2H分析
## 引入Decipher，即当2B是1时，以Decipher进行2H评价
def section_2_loss_2h(row, decipher_hi_gene_list):
    if row["2B"] == 1:
        for gene in row["AllGenes"]:
            if gene in decipher_hi_gene_list:
                row["2H"] = 1
                break
    return row

def section_2_loss_cal(df, db_dict):
    df = df.copy()

    # 数据库
    gene_info_df = db_dict["GENE_INFO"]
    gene_info_df['hi_score'] = pd.to_numeric(gene_info_df['hi_score'], errors='coerce')
    region_info_df = db_dict["REGION_INFO"]
    region_info_df['hi_score'] = pd.to_numeric(region_info_df['hi_score'], errors='coerce')
    cds_info_df = db_dict["CDS_INFO"]
    cds_info_df['hi_score'] = pd.to_numeric(cds_info_df['hi_score'], errors='coerce')
    utr_info_df = db_dict["UTR_INFO"]
    utr_info_df['hi_score'] = pd.to_numeric(utr_info_df['hi_score'], errors='coerce')
    exon_info_df = db_dict["EXON_INFO"]
    exon_info_df['hi_score'] = pd.to_numeric(exon_info_df['hi_score'], errors='coerce')
    decipher_df = db_dict["DECIPHER"]
    clingen_df = db_dict["CLINGEN_GENE_DISEASE"]

    gene_hi_df = gene_info_df[gene_info_df["hi_score"] == float(3)]
    gene_hi_vus_df = gene_info_df[(gene_info_df['hi_score'] != float(1)) & (gene_info_df['hi_score'] != float(2)) & (gene_info_df['hi_score'] != float(3)) & (gene_info_df['hi_score'] != float(40))]
    region_hi_df = region_info_df[region_info_df["hi_score"] == float(3)]
    region_hi_vus_df = region_info_df[(region_info_df['hi_score'] != float(1)) & (region_info_df['hi_score'] != float(2)) & (region_info_df['hi_score'] != float(3)) & (region_info_df['hi_score'] != float(40))]
    region_hi_b_df = region_info_df[region_info_df["hi_score"] == float(40)]
    cds_hi_df = cds_info_df[cds_info_df["hi_score"] == float(3)]
    utr_hi_df = utr_info_df[utr_info_df["hi_score"] == float(3)]
    exon_hi_df = exon_info_df[exon_info_df["hi_score"] == float(3)]
    decipher_hi_gene_list = decipher_df[decipher_df["hi_gene"] == "Y"]["gene"].tolist()

    # 对clingen df获得基因和MOI的字典
    gene_moi_dict = clingen_df.set_index("GENE_SYMBOL")["MOI"].to_dict()

    # 初始化
    section_2_list = ["2A", "2B", "2C-1", "2C-2", "2D-1", "2D-2", "2D-3", "2D-4", "2E", "2F", "2G", "2H"]
    for section in section_2_list:
        df.loc[:, section] = 0

    # 处理df
    df.loc[:, "2A"] = df.apply(lambda x: section_2_2a(x, gene_hi_df, region_hi_df), axis=1)
    df.loc[:, "2B"] = df.apply(lambda x: section_2_2b(x, gene_hi_vus_df, region_hi_vus_df), axis=1)

    # 其他证据等级
    df = df.apply(lambda x: section_2_loss_other_evi(x, cds_hi_df, utr_hi_df, exon_hi_df), axis=1)
    df = df.apply(lambda x: section_2_loss_2e(x, gene_info_df, gene_moi_dict), axis=1)
    df = df.apply(lambda x: section_2_loss_2f(x, region_hi_b_df), axis=1)
    df = df.apply(lambda x: section_2_loss_2g(x, region_hi_b_df), axis=1)
    df = df.apply(lambda x: section_2_loss_2h(x, decipher_hi_gene_list), axis=1)
    return df

# GAIN模式分析
## 2C分析
def section_2_gain_2c(row, gene_ts_b_df):
    overlap_gain_2c = gene_ts_b_df[gene_ts_b_df.apply(lambda x: common.overlap_similarity(row["#chrom"], row["start"], row["end"], x["# chrom"], x["start"], x["end"]) > 0.9, axis=1)]
    if len(overlap_gain_2c) > 1:
        row["2C"] = 1
    return row

## 2D分析
def section_2_gain_2d(row, gene_ts_b_df):
    overlap_gain_2d = gene_ts_b_df[gene_ts_b_df.apply(lambda x: common.complete_contain(x["# chrom"], x["start"], x["end"], row["#chrom"], row["start"], row["end"]), axis=1)]
    if len(overlap_gain_2d) > 1:
        row["2D"] = 1
    return row

# 2E分析 || 可能会中断蛋白质编码，所以这里实际上应该对断点进行注释，看看是否会有Stop Gain的情况出现
## 但这里因为不想引入其他数据库和注释器，暂不考虑加入2E分析
def section_2_gain_2e(row):
    return row

# 2F分析
def section_2_gain_2f(row, gene_ts_b_df):
    overlap_gain_2f = gene_ts_b_df[gene_ts_b_df.apply(lambda x: common.complete_contain(row["#chrom"], row["start"], row["end"], x["# chrom"], x["start"], x["end"]), axis=1)]
    if len(overlap_gain_2f) == 1:
        row["2F"] = 1
    return row

# 2G分析
def section_2_gain_2g(row, gene_ts_b_df):
    overlap_gain_2g = gene_ts_b_df[gene_ts_b_df.apply(lambda x: common.overlap_intersect(row["#chrom"], row["start"], row["end"], x["# chrom"], x["start"], x["end"]), axis=1)]
    if len(overlap_gain_2g) > 1:
        row["2G"] = 1
    return row

# 2H分析
def section_2_gain_2h(row, gene_hi_df):
    overlap_gain_2h = gene_hi_df[gene_hi_df.apply(lambda x: common.complete_contain(row["#chrom"], row["start"], row["end"], x["# chrom"], x["start"], x["end"]), axis=1)]
    if len(overlap_gain_2h) == 1:
        row["2H"] = 1
    return row

# 2I分析
def section_2_gain_2i(row, gene_info_df, gene_moi_dict):
    # 获取当前的2E
    tmp_2e = row['2E']
    # 复用
    row = section_2_loss_2e(row, gene_info_df, gene_moi_dict)
    tmp_2i = row['2E']
    row['2E'] = tmp_2e
    row['2I'] = tmp_2i
    return row

# 2JKL分析 || 需要联系患者表型，不建立
def section_2_gain_2j_k_l(row):
    return row

def section_2_gain_cal(df, db_dict):
    df = df.copy()

    # 数据库
    gene_info_df = db_dict["GENE_INFO"]
    gene_info_df['ts_score'] = pd.to_numeric(gene_info_df['ts_score'], errors='coerce')
    region_info_df = db_dict["REGION_INFO"]
    region_info_df['ts_score'] = pd.to_numeric(region_info_df['ts_score'], errors='coerce')

    gene_ts_df = gene_info_df[gene_info_df["ts_score"] == float(3)]
    gene_hi_df = gene_info_df[gene_info_df["hi_score"] == float(3)]
    gene_ts_vus_df = gene_info_df[(gene_info_df['ts_score'] != float(1)) & (gene_info_df['ts_score'] != float(2)) & (gene_info_df['ts_score'] != float(3)) & (gene_info_df['ts_score'] != float(40))]
    gene_ts_b_df = gene_info_df[gene_info_df["ts_score"] == float(40)]
    region_ts_df = region_info_df[region_info_df["ts_score"] == float(3)]
    region_ts_vus_df = region_info_df[(region_info_df['ts_score'] != float(1)) & (region_info_df['ts_score'] != float(2)) & (region_info_df['ts_score'] != float(3)) & (region_info_df['ts_score'] != float(40))]

    # 对clingen df获得基因和MOI的字典
    clingen_df = db_dict["CLINGEN_GENE_DISEASE"]
    gene_moi_dict = clingen_df.set_index("GENE_SYMBOL")["MOI"].to_dict()

    # 初始化
    section_2_list = ["2A", "2B", "2C", "2D", "2E", "2F", "2G", "2H", "2I", "2J", "2K", "2L"]
    for section in section_2_list:
        df.loc[:, section] = 0

    # 处理df
    df.loc[:, "2A"] = df.apply(lambda x: section_2_2a(x, gene_ts_df, region_ts_df), axis=1)
    df.loc[:, "2B"] = df.apply(lambda x: section_2_2b(x, gene_ts_vus_df, region_ts_vus_df), axis=1)

    # 处理其他证据
    df = df.apply(lambda x: section_2_gain_2c(x, gene_ts_b_df), axis=1)
    df = df.apply(lambda x: section_2_gain_2d(x, gene_ts_b_df), axis=1)
    df = df.apply(lambda x: section_2_gain_2e(x), axis=1)   
    df = df.apply(lambda x: section_2_gain_2f(x, gene_ts_b_df), axis=1) 
    df = df.apply(lambda x: section_2_gain_2g(x, gene_ts_b_df), axis=1) 
    df = df.apply(lambda x: section_2_gain_2h(x, gene_hi_df), axis=1)
    df = df.apply(lambda x: section_2_gain_2i(x, gene_info_df, gene_moi_dict), axis=1)
    df = df.apply(lambda x: section_2_gain_2j_k_l(x), axis=1)
    return df

# 合并
def section_2_cal(df, db_dict, mode="LOSS"):
    print("[Section2] Start!")
    df = df.copy()
    if mode == "GAIN":
        df = section_2_gain_cal(df, db_dict)
    else:
        df = section_2_loss_cal(df, db_dict)
    print("[Section2] Done!")
    return df

# end
