    echo "LOG:running bbmerge to identify adapters on $fileBasename"
	bbmerge.sh \
	in=$fileBasename"_1.fastq" \
	in2=$fileBasename"_2.fastq" \
	outadapter=$PWD/1_trim/adapters/$fileBasename.adapters.fa;