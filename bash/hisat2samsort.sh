## pipeline mapping + samsort


bash hisat2_1by1.sh $1
bash samsort.sh $1
bash featurecounts.sh
