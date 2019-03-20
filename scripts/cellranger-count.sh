module load cellranger/3.0.1

cellranger count \
    --id=Org1 \
    --transcriptome=refdata-cellranger-GRCh38-3.0.0 \
    --fastqs=fastqs \
    --sample=Org1 \
    --nosecondary \
    --localcores=20 \
    --localmem=64

