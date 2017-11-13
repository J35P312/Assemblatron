#accepts a bam file as first argument, prints an interleaved fastq to stdout
samtools view -h -F 2048 $1 | samtools view -F 1024 -Sh - | samtools view -Subh -F 256 - | samtools bam2fq -
