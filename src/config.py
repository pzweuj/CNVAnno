# coding=utf-8
# pzw
# 20241210

# 工具
BEDTOOLS = "bedtools"

# 原始数据库
ori_db = {
    "hg19": {
        "ClinGen_gene_curation": "ClinGen_gene_curation_list_GRCh37.tsv", # ftp://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh37.tsv
        "ClinGen_region_curation": "ClinGen_region_curation_list_GRCh37.tsv", # ftp://ftp.clinicalgenome.org/ClinGen_region_curation_list_GRCh37.tsv
        "gencode": "gencode.v47lift37.annotation.gff3.gz", # https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gff3.gz
        "clinvar": "clinvar_20241208.vcf.gz", # https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20241208.vcf.gz

    },
    "hg38": {
        "ClinGen_gene_curation": "ClinGen_gene_curation_list_GRCh38.tsv", # ftp://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv
        "ClinGen_region_curation": "ClinGen_region_curation_list_GRCh38.tsv", # ftp://ftp.clinicalgenome.org/ClinGen_region_curation_list_GRCh38.tsv
        "gencode": "gencode.v47lift37.annotation.gff3.gz", # https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gff3.gz
        "clinvar": "clinvar_20241208.vcf.gz", # https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20241208.vcf.gz

    },
    "mane": "MANE.GRCh38.v1.4.refseq_genomic.gff.gz", # https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4/MANE.GRCh38.v1.4.refseq_genomic.gff.gz
    "decipher": "HI_Predictions_Version3.bed.gz"# https://www.deciphergenomics.com/downloads/HI_Predictions_Version3.bed.gz
}

# 数据库
# hg19
hg19_db = {
    "GENE_INFO": "data/GRCh37.gene.info.bed",
    "REGION_INFO": "data/GRCh37.region.info.bed",
    "CDS_INFO": "data/GRCh37.cds.info.bed",
    "UTR_INFO": "data/GRCh37.utr.info.bed",
    "EXON_INFO": "data/GRCh37.exon.info.bed"
}

# hg38
hg38_db = {
    "GENE_INFO": "data/GRCh38.gene.info.bed",
    "REGION_INFO": "data/GRCh38.region.info.bed",
    "CDS_INFO": "data/GRCh38.cds.info.bed",
    "UTR_INFO": "data/GRCh38.utr.info.bed",
    "EXON_INFO": "data/GRCh38.exon.info.bed",

}

DECIPHER = "data/Decipher.sensitive.txt"
CLINGEN_GENE_DISEASE = "data/Clingen-Gene-Disease-Summary-2024-12-13.csv" # https://search.clinicalgenome.org/kb/gene-validity/download

