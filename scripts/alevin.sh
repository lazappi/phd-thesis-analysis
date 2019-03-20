module load salmon/0.12.0

salmon alevin \
    -l ISR \
    -1 fastqs/Orgs123/Org1/Org1_S1_L001_R1_001.fastq.gz \
    -2 fastqs/Orgs123/Org1/Org1_S1_L001_R2_001.fastq.gz \
    --chromium \
    -i transcriptomes/hg38/indexes/salmon-v0.12.0_v93.idx/ \
    -p 10 \
    -o Org1_alevin \
    --tgMap transcriptomes/hg38/hg38_tx2gene_versions_ensembl_v93.tsv \
    --mrna mt_genes.txt

