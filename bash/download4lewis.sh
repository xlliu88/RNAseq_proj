#!/bin/bash
# download fils from Lewis 
# command line: 
#    bash download4Lewis.sh filetype filename destination

server="liuxu@lewis.rnet.missouri.edu"
server_home="/home/liuxu/"
script_lewis="/home/liuxu/scripts/"
data_lewis="/storage/hpc/data/liuxu/rootcells2007/"

data_local="/mnt/d/RNAseq/rootcells2007/"
script_local="/mnt/d/RNAseq/scripts"
#echo $from
echo "arguments: $1 $2 $3"
if [[ $1 == "data" ]];
then
	des=$data_local
	src=$data_lewis
	
elif [[ $1 == "script" ]];
then
	des=$script_local
	src=$script_lewis
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

echo "downloading |||| $server:$src$2  ----> $des$desfile"
rsync -ahvP ${server}:$src$2 $des$desfile
