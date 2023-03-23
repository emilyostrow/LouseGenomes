module load spades

##filter for adapters
fastp -l 30 -i ssraw.fastq/ENO327FSo_S28_R1_001.fastq.gz -I ssraw.fastq/ENO327FSo_S28_R2_001.fastq.gz -o ssTrimmed/ENO327FSo_S28_R1.fastp.fastq -O ssTrimmed/ENO327FSo_S28_R2.fastp.fastq
##filter for length and genotype quality
/panfs/pfs.local/work/bi/bin/bbmap/bbduk.sh in1=../soTrimmed/ENO327FSo_S28_R1.fastp.fastq.gz in2=../soTrimmed/ENO327FSo_S28_R2.fastp.fastq.gz out1=ENO327FSo_S28_trimmedQual_R1.fastq.gz out2=ENO327FSo_S28_trimmedQual_R2.fastq.gz qtrim=r1 trimq=30 minlen=100 maq=30

##do an initial alignment using the adapter-filtered data as input
spades.py -1 ../ssTrimmed/ENO327FSo_S28_R1.fastp.fastq -2 ../ssTrimmed/ENO327FSo_S28_R1.fastp.fastq  -o ENO327FSo_S28_spades -t 12

##find the scaffold file in the spades output folder and make a copy of the scaffold file. Delete any contigs less than 10kb
##use the trimmed >10kb file as an untrusted contig and rerun spades with the quality and length filtered reads
spades.py -1 ../ssBbduk/ENO327FSo_S28_trimmedQual_R1.fastq.gz -2 ../ssBbduk/ENO327FSo_S28_trimmedQual_R2.fastq.gz --untrusted-contigs unfilteredSpades/oculatus_spades_scaffolds10kb.fasta -o ENO327FSo_S28_bbduk_10kbuntrusted -t 12
