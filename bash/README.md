#Pipeline

Scripts in order of execution:

1. `debarcode.sh` Use given barcodes to divide origin FASTQ files.
2. `map.sh` Align with BWA mem.
3. `dedup.sh` Deduplicate aligned bam files.