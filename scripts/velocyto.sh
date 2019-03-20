module load samtools/1.8
source activate velocyto

velocyto run \
    --bcfile barcodes_Org1.tsv \
    --outputfolder velocyto \
    --sampleid Org1 \
    --mask hg38_rep_mask.gtf \
    --dtype uint16 \
    --verbose \
    Org1/outs/possorted_genome_bam.bam \
    refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf
