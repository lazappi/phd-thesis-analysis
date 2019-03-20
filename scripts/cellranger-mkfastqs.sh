module load bcl2fastq/2.17.1.14
module load cellranger/3.0.1

cellranger mkfastq \
    --run=illumina_run \
    --id=mkfastqs \
    --csv=samples.csv \
    --ignore-dual-index \
    --qc \
    --output-dir=fastqs \
    --project=Orgs123 \
    --localcores=20 \
    --localmem=64 \
    --use-bases-mask=Y100n,I8n,y100n

