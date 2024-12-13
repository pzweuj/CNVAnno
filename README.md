# CNVAnno
ClinGen CNV Anno (Section1~3)

## Requirements

```
bedtools
python3
pandas
```

## Database
Edit the config.py file in the src folder, download the corresponding databases, and configure the paths for each database. Databases with commented links need to be downloaded manually, while databases without commented links will be generated using commands.

```bash
python3 cnvanno.py database -g hg19
```

## Running
Organize the input file in the bed format, and then use the commands below to add comments.

```
chr1  1000  2000  DUP
chr2  3333  999999  DEL
```

command
```bash
python3 cnvanno.py anno -i input.bed -o output.txt -g hg19
```


