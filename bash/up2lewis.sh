#!/bin/bash
# upload fils to Lewis 
# command line: 
#    bash up2lewis.sh filetype filename destination

server="liuxu@lewis.rnet.missouri.edu"
server_home="/home/liuxu/"
script_dir="/home/liuxu/scripts/"
data_dir="/storage/hpc/data/liuxu/rootcells2007/"

src="/mnt/d/RNAseq/rootcells2007/"
#echo $from
echo "arguments: $1 $2 $3"
if [[ $1 == "data" ]];
then
	des=$data_dir
	
elif [[ $1 == "script" ]];
then
	des=$script_dir
	src="./"
	dos2unix $2
else
	echo "Error: undefined file type"
	exit 1
fi

if [[ $3 -eq 0 ]]
then
    echo "!!! use source name as destination file name;"
    desfile=$2
else
    desfile=$3
fi

echo "uploading |||| $src$2  ----> ${server}:$des$desfile"
rsync -ahvP $src$2 ${server}:$des$desfile
