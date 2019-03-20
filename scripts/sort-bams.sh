module load samtools/1.8

SAMPLES=( Org1 Org2 Org3 )
for SAMPLE in "${SAMPLES[@]}"
do
    SAMPLE_PATH="$SAMPLE/outs/possorted_genome_bam.bam"
    OUT_PATH="$SAMPLE/outs/cellsorted_possorted_genome_bam.bam"

    samtools sort \
        -t CB \
        -m 3GB \
        --threads 9 \
        --output-fmt BAM \
        -o $OUT_PATH \
        $SAMPLE_PATH

done

