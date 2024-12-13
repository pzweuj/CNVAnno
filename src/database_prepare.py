# coding=utf-8
# pzw
# 20241204

import os
import gzip
import pandas as pd
from src import config

# MANE转录本对应字典
# https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4/MANE.GRCh38.v1.4.refseq_genomic.gff.gz
def mane_transcript_dict(mane_match_file):
    mane_dict = {}
    with gzip.open(mane_match_file, 'rt') as mane_match:
        for line in mane_match:
            if line.startswith('#'):
                continue
            lines = line.strip().split('\t')
            if lines[2] != 'mRNA':
                continue
            
            info = lines[8]
            attr_dict = {}
            for attr in info.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attr_dict[key] = value

            gene = attr_dict.get('gene', '.')

            mane_dict.setdefault(gene, {"refseq": ".", "ENST": "."})
            mane_dict[gene]['refseq'] = attr_dict.get('transcript_id', '.')

            # ensembl转录本
            enst = ""
            for enst in attr_dict.get('Dbxref', '.').split(','):
                if enst.startswith('Ensembl:'):
                    enst = enst.replace("Ensembl:", "")
                    break
            mane_dict[gene]['ENST'] = enst
    return mane_dict

# ClinGen基因剂量记录
# ftp://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh37.tsv
# ftp://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv
def clingene_dict(clingene_gene_curation):
    clingene_gene_dict = {}
    with open(clingene_gene_curation, 'r') as clingene_gene_curation:
        for line in clingene_gene_curation:
            if line.startswith('#'):
                continue
            lines = line.rstrip("\n").rstrip("\r").split('\t')
            gene = lines[0]
            # 使用三元运算符处理空值
            hi_score = lines[4] if lines[4] else '.'
            hi_des = lines[5] if lines[5] else '.'
            ts_score = lines[12] if lines[12] else '.'
            ts_des = lines[13] if lines[13] else '.'
            clingene_gene_dict.setdefault(gene, {
                'hi_score': hi_score,
                'hi_des': hi_des,
                'ts_score': ts_score,
                'ts_des': ts_des
            })
    return clingene_gene_dict

# 基因注释准备
# 只使用Mane，结果是后面计算1A时基本只会输出coding gene，不会输出all gene
# https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gff3.gz
# https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh37_mapping/gencode.v47lift37.annotation.gff3.gz
def gene_annotation_prepare(gff3_file, mane_match_file, clingene_gene_curation, bed_file_output):
    # 打开Mane匹配文件 || 注意这个文件只有GRCh38，但我们只用于转录本匹配
    mane_dict = mane_transcript_dict(mane_match_file)

    # 打开ClinGene基因注释文件，获取剂量记录
    clingene_gene_dict = clingene_dict(clingene_gene_curation)

    # 输出整理文件
    with gzip.open(gff3_file, 'rt') as gff3, open(bed_file_output, 'w') as bed:
        # 写入表头
        bed.write("# chrom\tstart\tend\tgene\tstrand\trefseq\tensembl\tgene_type\thi_score\thi_des\tts_score\tts_des\n")
        
        for line in gff3:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
                
            chrom, source, feature, start, end, score, strand, frame, attributes = fields
            
            # 解析属性字段
            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attr_dict[key] = value
                        
            if feature == 'transcript':
                gene = attr_dict.get('gene_name', '.')
                transcript_id = attr_dict.get('ID', '.')
                gene_type = attr_dict.get('gene_type', '.')
                trascript_dict = mane_dict.get(gene, {"refseq": ".", "ENST": "."})
                
                # 基因在Mane，这个时候如果转录本不是Mane，就跳过
                if gene in mane_dict:
                    # 查询当前的转录本是否与ENST匹配
                    if not transcript_id.split('.')[0] in trascript_dict['ENST']:
                        continue
                else:
                    # 不在Mane中的话，就需要是tag=Ensembl_canonical
                    if 'Ensembl_canonical' not in attr_dict.get('tag', '.'):
                        continue

                # 获得clingene基因注释
                gene_dict = clingene_gene_dict.get(gene, {'hi_score': '.', 'hi_des': '.', 'ts_score': '.', 'ts_des': '.'})
                hi_score = gene_dict.get('hi_score', '.')
                hi_des = gene_dict.get('hi_des', '.')
                ts_score = gene_dict.get('ts_score', '.')
                ts_des = gene_dict.get('ts_des', '.')

                # 获得mane基因注释
                refseq = trascript_dict['refseq']
                enst = transcript_id

                # 输出
                bed.write(f"{chrom}\t{start}\t{end}\t{gene}\t{strand}\t{refseq}\t{enst}\t{gene_type}\t{hi_score}\t{hi_des}\t{ts_score}\t{ts_des}\n")



# TS\HI基因和基因区域
# ftp://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh37.tsv
# ftp://ftp.clinicalgenome.org/ClinGen_region_curation_list_GRCh37.tsv
# 需要分开处理，基因水平需要计算到CDS区域、Exon区域，这样方便后续计算
def ts_hi_region_prepare(ts_hi_region_file, ts_hi_region_output_file):
    # 读取区域文件
    with open(ts_hi_region_output_file, 'w', encoding='utf-8') as ts_hi_region_output:
        # 标题
        ts_hi_region_output.write("# chrom\tstart\tend\thi_score\thi_des\tts_score\tts_des\tISCA ID\n")
        
        # 读取区域文件
        with open(ts_hi_region_file, 'r') as ts_hi_region:
            for line in ts_hi_region:
                if line.startswith("#"):
                    continue
                lines = line.rstrip("\n").rstrip("\r").split('\t')
                region = lines[3]
                
                if region == 'tbd':
                    continue

                chrom = region.split(':')[0]
                start = region.split(':')[1].split('-')[0]
                end = region.split(':')[1].split('-')[1]
                hi_score = lines[4]
                hi_des = lines[5]
                ts_score = lines[12]
                ts_des = lines[13]
                isca_id = lines[0]

                # 写入文件
                ts_hi_region_output.write(f"{chrom}\t{start}\t{end}\t{hi_score}\t{hi_des}\t{ts_score}\t{ts_des}\t{isca_id}\n")

# 外显子水平注释文件 || 只要Mane转录本，其他都不要了
def exon_annotation_prepare(gff3_file, mane_match_file, clingene_gene_curation, bed_file_output, cds_mode=False):
    mane_dict = mane_transcript_dict(mane_match_file)
    clingene_gene_dict = clingene_dict(clingene_gene_curation)

    # 输出整理文件
    with gzip.open(gff3_file, 'rt') as gff3, open(bed_file_output, 'w') as bed:
        # 写入表头
        bed.write("# chrom\tstart\tend\tgene\tstrand\trefseq\tensembl\tgene_type\tExon_Num\tFirst_Exon\tLast_Exon\thi_score\thi_des\tts_score\tts_des\n")
        
        # 必须遍历一遍才能获得最大的exon
        all_gene_exon_num = {}
        for line in gff3:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
                
            chrom, source, feature, start, end, score, strand, frame, attributes = fields
            
            # 解析属性字段
            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attr_dict[key] = value
            
            # 获得基因注释
            check_feature = "CDS" if cds_mode else "exon"
            if feature == check_feature:
                gene = attr_dict.get('gene_name', '.')
                transcript_id = attr_dict.get('transcript_id', '.')
                gene_type = attr_dict.get('gene_type', '.')
                trascript_dict = mane_dict.get(gene, {"refseq": ".", "ENST": "."})
                exon_num = attr_dict.get('exon_number', '.')

                # 如果没有exon_num, 则跳过
                if exon_num == '.':
                    continue

                # 基因在Mane，这个时候如果转录本不是Mane，就跳过
                if gene in mane_dict:
                    # 查询当前的转录本是否与ENST匹配
                    if not transcript_id.split('.')[0] in trascript_dict['ENST']:
                        continue
                else:
                    continue
            
                # 初始化字典
                gene_dict = clingene_gene_dict.get(gene, {'hi_score': '.', 'hi_des': '.', 'ts_score': '.', 'ts_des': '.'})
                hi_score = gene_dict.get('hi_score', '.')
                hi_des = gene_dict.get('hi_des', '.')
                ts_score = gene_dict.get('ts_score', '.')
                ts_des = gene_dict.get('ts_des', '.')
                all_gene_exon_num.setdefault(gene,
                    {
                        "enst": trascript_dict['ENST'], "refseq": trascript_dict['refseq'],
                        "first_exon": "999", "last_exon": "1",
                        "exon_list": {},
                        "hi_score": hi_score, "hi_des": hi_des, 
                        "ts_score": ts_score, "ts_des": ts_des,
                        "gene_type": gene_type, "strand": strand
                    }
                )

                # 判断是否第一个exon
                if int(exon_num) < int(all_gene_exon_num[gene]['first_exon']):
                    all_gene_exon_num[gene]['first_exon'] = exon_num

                # 判断是否最后一个exon
                if int(exon_num) > int(all_gene_exon_num[gene]['last_exon']):
                    all_gene_exon_num[gene]['last_exon'] = exon_num

                # 塞进所有内容
                all_gene_exon_num[gene]['exon_list'][exon_num] = {
                    "chrom": chrom, "start": start, "end": end
                }

        # 处理结果
        for gene, exon_dict in all_gene_exon_num.items():
            # 输出结果
            for exon in exon_dict['exon_list'].keys():
                location_dict = exon_dict['exon_list'][exon]
                chrom, start, end = location_dict['chrom'], location_dict['start'], location_dict['end']
                first_exon = "Y" if exon == exon_dict['first_exon'] else "N"
                last_exon = "Y" if exon == exon_dict['last_exon'] else "N"

                bed.write(f"{chrom}\t{start}\t{end}\t{gene}\t{exon_dict['strand']}\t{exon_dict['refseq']}\t{exon_dict['enst']}\t{exon_dict['gene_type']}\t{exon}\t{first_exon}\t{last_exon}\t{exon_dict['hi_score']}\t{exon_dict['hi_des']}\t{exon_dict['ts_score']}\t{exon_dict['ts_des']}\n")


# 从外显子水平和CDS水平注释文件中，获得UTR区域注释文件
# 使用外显子的last_exon，和CDS的last_exon，获得UTR3区域
# 使用外显子的first_exon，和CDS的first_exon，获得UTR5区域
def utr_annotation_prepare(exon_file, cds_file, output_file):
    gene_dict = {}
    
    exon_df = pd.read_csv(exon_file, sep='\t', header=0)
    first_exon_df = exon_df[exon_df['First_Exon'] == 'Y']
    last_exon_df = exon_df[exon_df['Last_Exon'] == 'Y']
    for index, row in first_exon_df.iterrows():
        gene = row['gene']
        chrom = row['# chrom']
        start = row['start']
        end = row['end']
        strand = row['strand']
        refseq = row['refseq']
        ensembl = row['ensembl']
        hi_score = row['hi_score']
        hi_des = row['hi_des']
        ts_score = row['ts_score']
        ts_des = row['ts_des']
        gene_dict.setdefault(gene, {"chrom": chrom, "first_exon_start": start, "first_exon_end": end, "strand": strand, "refseq": refseq, "ensembl": ensembl, "hi_score": hi_score, "hi_des": hi_des, "ts_score": ts_score, "ts_des": ts_des})
    for index, row in last_exon_df.iterrows():
        gene = row['gene']
        if gene in gene_dict:
            start = row['start']
            end = row['end']
            gene_dict[gene]['last_exon_start'] = start
            gene_dict[gene]['last_exon_end'] = end
    
    cds_df = pd.read_csv(cds_file, sep='\t', header=0)
    first_cds_df = cds_df[cds_df['First_Exon'] == 'Y']
    last_cds_df = cds_df[cds_df['Last_Exon'] == 'Y']
    for index, row in first_cds_df.iterrows():
        gene = row['gene']
        if gene in gene_dict:
            start = row['start']
            end = row['end']
            gene_dict[gene]['first_cds_start'] = start
            gene_dict[gene]['first_cds_end'] = end
    for index, row in last_cds_df.iterrows():
        gene = row['gene']
        if gene in gene_dict:
            start = row['start']
            end = row['end']
            gene_dict[gene]['last_cds_start'] = start
            gene_dict[gene]['last_cds_end'] = end
    
    # 计算UTR3和UTR5
    for gene in gene_dict:
        # 有些基因不是所有信息都有，也不要了
        if len(gene_dict[gene]) >= 16:
            gene_dict[gene]['utr3_start'] = gene_dict[gene]['last_cds_end']
            gene_dict[gene]['utr3_end'] = gene_dict[gene]['last_exon_end']
            gene_dict[gene]['utr5_start'] = gene_dict[gene]['first_exon_start']
            gene_dict[gene]['utr5_end'] = gene_dict[gene]['first_cds_start']
            if gene_dict[gene]['strand'] == "-":
                gene_dict[gene]['utr3_start'] = gene_dict[gene]['last_exon_start']
                gene_dict[gene]['utr3_end'] = gene_dict[gene]['last_cds_start']
                gene_dict[gene]['utr5_start'] = gene_dict[gene]['first_cds_end']
                gene_dict[gene]['utr5_end'] = gene_dict[gene]['first_exon_end']

    # 输出结果
    with open(output_file, 'w') as output:
        output.write("# chrom\tstart\tend\tlocation\tgene\trefseq\tensembl\tstrand\thi_score\thi_des\tts_score\tts_des\n")
        for gene, utr_dict in gene_dict.items():
            # 有些基因不是所有信息都有，也不要了
            if len(utr_dict) < 20:
                continue
            
            # 不存在UTR的基因，也过滤掉
            if (utr_dict['utr3_start'] == utr_dict['utr3_end']) or (utr_dict['utr5_start'] == utr_dict['utr5_end']):
                continue
            
            # 输出
            if utr_dict['strand'] == "+":
                output.write(f"{utr_dict['chrom']}\t{utr_dict['utr5_start']}\t{utr_dict['utr5_end']}\tUTR5\t{gene}\t{utr_dict['refseq']}\t{utr_dict['ensembl']}\t{utr_dict['strand']}\t{utr_dict['hi_score']}\t{utr_dict['hi_des']}\t{utr_dict['ts_score']}\t{utr_dict['ts_des']}\n")
                output.write(f"{utr_dict['chrom']}\t{utr_dict['utr3_start']}\t{utr_dict['utr3_end']}\tUTR3\t{gene}\t{utr_dict['refseq']}\t{utr_dict['ensembl']}\t{utr_dict['strand']}\t{utr_dict['hi_score']}\t{utr_dict['hi_des']}\t{utr_dict['ts_score']}\t{utr_dict['ts_des']}\n")
            else:
                output.write(f"{utr_dict['chrom']}\t{utr_dict['utr3_start']}\t{utr_dict['utr3_end']}\tUTR3\t{gene}\t{utr_dict['refseq']}\t{utr_dict['ensembl']}\t{utr_dict['strand']}\t{utr_dict['hi_score']}\t{utr_dict['hi_des']}\t{utr_dict['ts_score']}\t{utr_dict['ts_des']}\n")
                output.write(f"{utr_dict['chrom']}\t{utr_dict['utr5_start']}\t{utr_dict['utr5_end']}\tUTR5\t{gene}\t{utr_dict['refseq']}\t{utr_dict['ensembl']}\t{utr_dict['strand']}\t{utr_dict['hi_score']}\t{utr_dict['hi_des']}\t{utr_dict['ts_score']}\t{utr_dict['ts_des']}\n")

# 外显子水平与Clinvar交集，需要找到致病区域
# https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20241208.vcf.gz
# https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20241208.vcf.gz
def exon_clinvar_prepare(exon_file, clinvar_file):
    # 先对clinvar进行处理，仅保留Pathogenic和Likely_pathogenic
    BEDTOOLS = config.BEDTOOLS
    cmd = f"""
        zcat {clinvar_file} | grep -v '^#' | grep -E 'Pathogenic|Likely_pathogenic' | awk '{{print "chr"$1"\t"$2"\t"$2+length($4)-1}}' > {clinvar_file}.bed
        {BEDTOOLS} intersect -a {exon_file} -b {clinvar_file}.bed -c > {exon_file}.clinvar.bed
        head -n 1 {exon_file} | sed 's/$/\\tclinvar_count/' > {exon_file}.new.header
        cat {exon_file}.new.header {exon_file}.clinvar.bed > {exon_file}
        rm {exon_file}.new.header {exon_file}.clinvar.bed {clinvar_file}.bed
    """
    os.system(cmd)

# Decipher的剂量敏感文件
def decipher_sensitive_prepare(decipher_file, output_file):
    with open(output_file, 'w') as output:
        output.write("gene\thi_score\thi_index\thi_gene\n")
    
        # 读取文件
        with gzip.open(decipher_file, 'rt') as decipher:
            for line in decipher:
                if line.startswith('track name'):
                    continue
                lines = line.strip().split('\t')
                gene, hi_score, hi_index = lines[3].split("|")
                hi_gene = "Y" if float(hi_score) >= 0.9 else "N"
                output.write(f"{gene}\t{hi_score}\t{hi_index}\t{hi_gene}\n")

# ClinGen的基因疾病关联文件
# https://search.clinicalgenome.org/kb/gene-validity/download
# 注意，这个文件会以下载的日期进行命名，因此建议自行设置好名称

# 使用测试
# 准备基因注释文件
def prepare_all(genome_version):
    mane = config.ori_db["mane"]
    decipher = config.ori_db["decipher"]
    oridb_genome = config.ori_db[genome_version]
    predb_genome = config[genome_version + "_db"]

    gene_annotation_prepare(oridb_genome["gencode"],
                            mane,
                            oridb_genome["ClinGen_gene_curation"],
                            predb_genome["GENE_INFO"])

    # 准备区域注释文件
    ts_hi_region_prepare(oridb_genome["ClinGen_region_curation"], predb_genome["REGION_INFO"])

    # 准备外显子注释文件
    exon_annotation_prepare(oridb_genome["gencode"],
                            mane,
                            oridb_genome["ClinGen_gene_curation"],
                            predb_genome["EXON_INFO"])

    exon_annotation_prepare(oridb_genome["gencode"],
                            mane,
                            oridb_genome["ClinGen_gene_curation"],
                            predb_genome["CDS_INFO"], cds_mode=True)
    exon_clinvar_prepare(predb_genome["EXON_INFO"], oridb_genome["clinvar"])

    # 准备UTR注释文件
    utr_annotation_prepare(predb_genome["EXON_INFO"], predb_genome["CDS_INFO"], predb_genome["UTR_INFO"])

    # 准备Decipher的剂量敏感文件
    decipher_sensitive_prepare(decipher, config.DECIPHER)

# end
