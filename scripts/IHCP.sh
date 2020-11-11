#docker run -it -v $PWD:$PWD --name IHCP continuumio/anaconda3:latest bash
apt-get update
conda install -y -c bioconda bbmap fastqc bowtie2;
apt install -y lbzip2 samtools;
pip install multiqc;
WORK=/media/jochum00/Aagaard_Raid3/jochum_IHCP;

cd $WORK;


##############ADAPTER TRIMMING##############
mkdir -p $PWD/1_trim/adapters/
mkdir -p $PWD/1_trim/trimfail/
mkdir -p $PWD/1_trim/trimmed/paired/
mkdir -p $PWD/1_trim/trimmed/singletons/

for um1 in $(cat manifest)
do
	echo "decompressing $um1";
	lbzip2 -vdk $PWD/0_raw/$um1"_1.fastq.bz2";
	lbzip2 -vdk $PWD/0_raw/$um1"_2.fastq.bz2";
	echo "LOG:running bbmerge to identify adapters on $um1";
	bbmerge.sh \
	in=$PWD/0_raw/$um1"_1.fastq" \
	in2=$PWD/0_raw/$um1"_2.fastq" \
	outadapter=$PWD/1_trim/adapters/$um1".adapters.fa";
	echo "LOG:running bbduk adapter trimmimg on $um1";
    bbduk.sh \
	t=48 \
	in=$PWD/0_raw/$um1"_1.fastq" \
	in2=$PWD/0_raw/$um1"_2.fastq" \
	out1=$PWD/1_trim/trimmed/paired/$um1".trim_1.fastq" \
	out2=$PWD/1_trim/trimmed/paired/$um1".trim_2.fastq" \
	outm=$PWD/1_trim/trimfail/$um1".trimfail_1.fastq" \
	outm2=$PWD/1_trim/trimfail/$um1".trimfail_2.fastq" \
	outs=$PWD/1_trim/trimmed/singletons/$um1".trimmed.singletons.fastq" \
	ref=$PWD/1_trim/adapters/$um1.adapters.fa ktrim=r k=23 mink=11 hdist=1 tbo tpe qtrim=rl trimq=10 minlen=20;
done


##########QUALITY CONTROL#################
mkdir -p $PWD/1_trim/trimmed/paired/fastqc
mkdir -p $PWD/1_trim/trimmed/singletons/fastqc
fastqc --threads 48 --outdir $PWD/1_trim/trimmed/paired/fastqc mkdir $PWD/1_trim/trimmed/paired/*.fastq;
fastqc --threads 48 --outdir $PWD/1_trim/trimmed/singletons/fastqc $PWD/1_trim/trimmed/singletons/*.fastq;
multiqc 1_trim/trimmed/paired/ -o 1_trim/trimmed/paired/fastqc/
multiqc 1_trim/trimmed/singletons/ -o 1_trim/trimmed/singletons/fastqc/

jochum00/bowtie2_samtools_updated_versions 
###REMOVING HUMAN AND PHIX CONTAMINANTS##
mkdir -p $PWD/2_filter/aligned;
mkdir -p $PWD/2_filter/unaligned;
mkdir -p $PWD/2_filter/bam;
for um1 in $(cat manifest)
do
	echo "LOG: running bowtie2 on $um1";
	bowtie2 \
	--threads 48 \
	--seed 3141 \
	--very-sensitive \
	--un-conc $PWD/2_filter/unaligned/$um1.nonhuman_%.fastq \
	--al-conc $PWD/2_filter/aligned/$um1.human_%.fastq \
	--un  $PWD/2_filter/unaligned/$um1.singletons.nonhuman.fastq \
	--al $PWD/2_filter/aligned/$um1.singletons.human.fastq \
	-x /media/jochum00/Aagaard_Raid/reference_datasets/kneaddata/hg37/hg37dec_v0.1 \
	-1 $PWD/1_trim/trimmed/paired/$um1.trim_1.fastq \
	-2 $PWD/1_trim/trimmed/paired/$um1.trim_2.fastq \
	-U $PWD/1_trim/trimmed/singletons/$um1.trimmed.singletons.fastq |samtools view -bS -@48 > $PWD/2_filter/bam/$um1".bam";
	
	#-S $PWD/2_filter/sam/$um1".SE.sam";
	###Maybe do this later if you want the human reads###
	# echo converting $um1;
	# samtools view -@48 -bS $PWD/2_filter/sam/$um1".sam" > $PWD/2_filter/sam/$um1".bam";
	# echo sorting $um1;
	# samtools sort -@48 $PWD/2_filter/sam/$um1".bam" -o $PWD/2_filter/sam/$um1".sorted.bam";
	# echo indexing $um1;
	# samtools index -@48 $PWD/2_filter/sam/$um1".sorted.bam" $PWD/2_filter/sam/$um1".sorted.bai";
	echo "bowtie2 complete";
done
	
###TAXONOMIC CLASSIFICATION###
#docker run -it -v /media/jochum00:/media/jochum00 --name metaphlan3 biobakery/metaphlan:latest bash
###METAPHLAN3##
mkdir -p $PWD/3_metaphlan3/sams;
mkdir -p $PWD/3_metaphlan3/profiles;
mkdir -p $PWD/3_metaphlan3/bowtie2;
echo "Running metaphlan 3.0 on ${n}"

for fileBasename in $(cat manifest);
for fileBasename in $li;
do
	
	echo "Running metaphlan 3.0 on  $fileBasename";
	metaphlan \
	--force \
	$PWD/2_filter/unaligned/$fileBasename.nonhuman_1.fastq,\
	$PWD/2_filter/unaligned/$fileBasename.nonhuman_2.fastq,\
	$PWD/2_filter/unaligned/$fileBasename.singletons.nonhuman.fastq \
	--bowtie2db /media/jochum00/Aagaard_Raid/reference_datasets/metaphlan3/mpa_v30_CHOCOPhlAn_201901 \
	--input_type fastq \
	-s $PWD/3_metaphlan3/sams/$fileBasename.sam \
	--bowtie2out $PWD/3_metaphlan3/bowtie2/$fileBasename.bowtie2 \
	-o $PWD/3_metaphlan3/profiles/$fileBasename.profile.tsv \
	--nproc $(nproc)
	echo "metaphlan3 complete on $fileBasename";
	#echo "converting sam to bam on $fileBasename";
	#samtools view -bS -@48 $PWD/3_metaphlan3/sams/$fileBasename.sam> $PWD/3_metaphlan3/sams/$fileBasename.bam
done

###############KRAKEN2/BRACKEN##################################

#docker run -it -v /media/jochum00:/media/jochum00 --name k2 tbattaglia/kraken2 bash
WORK=/media/jochum00/Aagaard_Raid3/jochum_IHCP;
cd $WORK;

mkdir -p $PWD/4_kraken2/classified;
mkdir -p $PWD/4_kraken2/unclassified;
mkdir -p $PWD/4_kraken2/reports;
mkdir -p $PWD/4_kraken2/kraken2_out;

for n in  in $(cat manifest);
do
	echo "Running Kraken2 on $n";
	kraken2 \
	--db /media/jochum00/Jochum_Raid/reference_datasets/kraken2_std2 \
	--paired \
	--use-names \
	--threads $(nproc) \
	--classified-out  $PWD/4_kraken2/classified/$n.classified#.fastq \
	--unclassified-out $PWD/4_kraken2/unclassified/$n.unclassified#.fastq \
	--output 4_kraken2/kraken2_out/$n.kraken2 \
	--report 4_kraken2/reports/$n.kreport \
	$PWD/2_filter/unaligned/$n.nonhuman_1.fastq \
	$PWD/2_filter/unaligned/$n.nonhuman_2.fastq;
#	$PWD/2_filter/unaligned/$n.singletons.nonhuman.fastq; 
	echo "Kraken2 complete on $n";
done

###############BRACKEN################
for n in $(cat  manifest) 
do
bracken \
-d  /media/jochum00/Aagaard_Raid/reference_datasets/kraken2_db/ \
-i 4_kraken2/reports/$n.kreport \
-o 5_bracken/bracken_out/$n.bracken \
-w 5_bracken/reports/$n.breport \
-r 100 \
-t 1 \
-l S
done

###############FILTER_BRACKEN.PY################
for n in $(cat manifest); 
do
#NOTE: There are still human reads in this output that somehow passed through kneaddata
# We will now extract these reads and phix-174 using KrakenTools with taxid 2759 (Eukaryota),9606(human), and 10847(phix) 
python /media/jochum00/Aagaard_Raid/reference_datasets/KrakenTools/filter_bracken.out.py \
-i  5_bracken/bracken_out/$n.bracken \
-o  5_bracken/bracken_out/filt/$n.filt.bracken \
--exclude 9606 2759 10847
done


##############BIT-COMBINE-AND-ADD-LINEAGE###############
for n in $(cat manifest); 
do
# putting a few of the bracken output file names into a file 
ls 5_bracken/bracken_out/filt/$n.filt.bracken>>files;
cat files|cut -f1 -d ".">names;
paste -d "\t" files names>>5_bracken/filt.bracken-outputs.tsv;
rm files;
rm names;
# removing the info that is the same for all in this case
done

3.25 2.375 2.375 saving $85.00
###############EXTRACT_KRAKEN2_READS.PY################
ls |cut -f2 -d".">names
/media/jochum00/Aagaard_Raid/reference_datasets/KrakenTools/combine_kreports.py -r $(ls *.kreport) -o combined.kreport --no-headers --sample-names $(cat names)



for n 
kraken2 --db KRAKEN2DB --threads THREADNUM --report MYSAMPLE.KREPORT \
    --paired SAMPLE_1.FASTA SAMPLE_2.FASTA > MYSAMPLE.KRAKEN2
python kreport2mpa.py -r MYSAMPLE.KREPORT -o MYSAMPLE.MPA.TXT 


#NOTE: Now we convert the bracken report to mpa style
kreport2mpa.py \
-r 23862_6.paired.std.filt.breport \
-o 23862_6.paired.std.filt.mpa.percentages.breport \
--display-header \
--percentages \
--no-intermediate-ranks