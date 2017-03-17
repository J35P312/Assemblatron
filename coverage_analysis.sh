#analyse the coverage across the genome
#first argument bam file of aligned contigs
#second argument - a gene list
#third argument output prefix
#requires bedtools and tiddit =P

TIDDIT_path=/home/jesper/Assemblatron/TIDDIT/bin/TIDDIT

$TIDDIT_path --cov -b $1 -o $3 -z 50
python compute_coverage.py $3.tab  > $3_genome_wide.tab

intersectBed -a $3.tab -b $2 > $3.panel.tab
python compute_coverage.py $3.tab  > $3_panel_wide.tab



