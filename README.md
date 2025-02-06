# Dante Labs VCF file converter for Genetic Lifehacks
This project provides Python tools to convert whole genome sequencing Dante Labs VCF files into txt format for use in the Genetic Lifehacks website. 
Python coding knowlege and Python coding environment is needed to use the scripts.

## Preparation
1. Get your 'filtered.snp.vcf.gz' and 'filtered.snp.vcf.gz.tbi' files from Dante Labs whole genome sequencing.
2. If your VCF file is not annotated (meaning it does not have rs###### numbers assigned to rows in your VCF file), download dbSNP files for VCF annotation:

```
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz.tbi
```

## Running annotation
annotate.py does the Dante Labs VCF annotation against dbSNP SNP information. Adjust the script to reflect your input and output file names and locations.
Note: this takes a long time (over an hour).

## Running conversion
convert.py takes your VCF file annotated with rs#### numbers and restructures into the format used by the Genetic Lifehacks website. The resulting TXT file output is what you use in the Genetic Lifehacks website. Adjust the script to reflect your input and output file names and locations.
