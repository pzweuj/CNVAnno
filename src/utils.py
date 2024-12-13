# coding=utf-8
# pzw
# 20241210
# CNV ClinGen注释
# 参考
# https://cnvcalc.clinicalgenome.org/cnvcalc/cnv-loss
# https://cnvcalc.clinicalgenome.org/cnvcalc/cnv-gain

import sys
from src.config import hg19_db, hg38_db, DECIPHER, CLINGEN_GENE_DISEASE
import pandas as pd
from src.section_1 import section_1_cal
from src.section_2 import section_2_cal
from src.section_3 import section_3_cal

# 基础处理
# 格式如：chr1 1000 2000 DEL
def basic_process(bed_file):
    df = pd.read_csv(bed_file, sep='\t', header=0)

    # 如果列数大于4，报错退出
    if len(df.columns) > 4:
        print("Error: The input file has more than 4 columns.")
        sys.exit(1)
    
    # 检查标题，如果不是#chrom开头，则文件可能没有标题，此时将标题行设置为第一行，然后补充标题
    if not df.columns[0].startswith('#chrom'):
        df = pd.read_csv(bed_file, sep='\t', header=None)
        df.columns = ['#chrom', 'start', 'end', 'type']

    # 调整#chrom列，如果没有以chr开头，则加上chr，如果是MT，则改为chrM
    df['#chrom'] = df['#chrom'].apply(lambda x: 'chr' + str(x) if not str(x).startswith('chr') else str(x).replace('MT', 'chrM'))
    
    # type列，需要将Del，Dup, Deletion, Duplication等转换为DEL，DUP
    df['type'] = df['type'].apply(lambda x: 'DEL' if 'Del' in str(x) else ('DUP' if 'Dup' in str(x) else str(x)))
    return df

# 计算器
def risk_cal(row, mode="LOSS"):
    score = 0
    evidence_list = []
    point_dict = {
        "LOSS": {
            "section1": {"1A": 0, "1B": -0.6},
            "section2": {"2A": 1, "2B": 0,
                "2C-1": 0.9, "2C-2": 0.45,
                "2D-1": 0, "2D-2": 0.9, "2D-3": 0.45, "2D-4": 1,
                "2E": {"PVS1": 0.9, "PVS1_Strong": 0.45, "PVS1_Moderate": 0.3, "PVS1_Supporting": 0.15},
                "2F": -1, "2G": 0, "2H": 0.15},
            "section3": {"3A": 0, "3B": 0.45, "3C": 0.9}
        },
        "GAIN": {
            "section1": {"1A": 0, "1B": -0.6},
            "section2": {"2A": 1, "2B": 0,
                "2C": -1, "2D": -1, "2E": 0, "2F": -0.9,
                "2G": 0, "2H": 0,
                "2I": {"PVS1": 0.9, "PVS1_Strong": 0.45},
                "2J": 0, "2K": 0.45, "2L": 0},
            "section3": {"3A": 0, "3B": 0.45, "3C": 0.9}    
        }
    }

    if row["1A"] == 0:
        row["Evidences"] = "1B"
        row["Score"] = -0.6
        return row
    
    evidence_list.append("1A")

    point_dict_section2 = point_dict[mode]["section2"]
    if mode == "LOSS":
        for p in point_dict_section2:
            point_tmp = 0
            if p == '2E':
                point_tmp = point_dict_section2['2E'].get(row['2E'], 0)
            else:
                point_tmp = row[p]
            if point_tmp > 0:
                evidence_list.append(p)
                score += point_dict_section2[p]
                break
    else:
        for p in point_dict_section2:
            point_tmp = 0
            if p == '2I':
                point_tmp = point_dict_section2['2I'].get(row['2I'], 0)
            else:
                point_tmp = row[p]
            if point_tmp > 0:
                evidence_list.append(p)
                score += point_dict_section2[p]
                break

    point_dict_section3 = point_dict[mode]["section3"]
    for p in point_dict_section3:
        if row[p] > 0:
            evidence_list.append(p)
            score += point_dict_section3[p]
            break
    row["Score"] = float(score)
    row["Evidences"] = ";".join(evidence_list)
    return row

# 判断器
def judge(row):
    conclusion = "Pathogenic"
    score = row['Score']
    if score < -0.9:
        conclusion = "Benign"
    elif score < -0.85:
        conclusion = "Likely_benign"
    elif score < 0.9:
        conclusion = "VUS"
    elif score < 0.95:
        conclusion = "Likely_pathogenic"
    row["Conclusion"] = conclusion
    return row

# 测试处理
def cnv_anno_process(input_bed, output_txt, genome='hg19'):
    use_db = hg19_db
    if genome == 'hg38':
        use_db = hg38_db

    # 避免数据库的重复读取，在这里一次性读取
    db_dict = {}
    for key, value in use_db.items():
        db_dict[key] = pd.read_csv(value, sep='\t', header=0)
    db_dict["DECIPHER"] = pd.read_csv(DECIPHER, sep='\t', header=0)
    clingen_header =  [
        "GENE_SYMBOL", "GENE_ID", "DISEASE_LABEL", "DISEASE_ID", 
        "MOI", "SOP", "CLASSIFICATION", "ONLINE_REPORT", 
        "CLASSIFICATION_DATE", "GCEP"
    ]
    db_dict["CLINGEN_GENE_DISEASE"] = pd.read_csv(CLINGEN_GENE_DISEASE, sep=",", skiprows=6, names=clingen_header)

    # 开始分析
    df = basic_process(input_bed)

    # 初始化
    df.loc[:, "Evidences"] = "UnKnown"
    df.loc[:, "Score"] = 0
    df.loc[:, "Conclusion"] = "UnKnown"
    
    # section1 统一分析
    df_section1 = section_1_cal(df, db_dict["GENE_INFO"])

    # section2 开始要拆分DEL和DUP分别进行分析
    df_del = df_section1[df_section1["type"] == "DEL"]
    df_dup = df_section1[df_section1["type"] == "DUP"]
    df_del = section_2_cal(df_del, db_dict, "LOSS")
    df_dup = section_2_cal(df_dup, db_dict, "GAIN")

    # section3
    df_del = section_3_cal(df_del, db_dict["GENE_INFO"], "LOSS")
    df_dup = section_3_cal(df_dup, db_dict["GENE_INFO"], "GAIN")

    # 计算
    df_del = df_del.apply(lambda x: risk_cal(x, "LOSS"), axis=1)
    df_dup = df_dup.apply(lambda x: risk_cal(x, "GAIN"), axis=1)
    df_del = df_del.apply(judge, axis=1)
    df_dup = df_dup.apply(judge, axis=1)

    # 需要的列
    need_col = [
        "#chrom", "start", "end", "type",
        "Evidences", "Score", "Conclusion",
        "CodingGenes", "AllGenes", "CodingGenes_Num", "AllGenes_Num"
    ]
    df_del = df_del[need_col]
    df_dup = df_dup[need_col]
    df_combined = pd.concat([df_del, df_dup], ignore_index=True)
    df_sorted = df_combined.sort_values(by=['#chrom', 'start', 'end'])
    df_sorted.reset_index(drop=True, inplace=True)

    # 输出
    df_sorted.to_csv(output_txt, sep="\t", index=False, header=True)
