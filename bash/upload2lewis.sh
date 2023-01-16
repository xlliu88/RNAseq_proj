#!/bin/bash
# move script to Lewis
# command line: 
#    bash upload2lewis.sh filetype filename destination

server="liuxu@lewis.rnet.missouri.edu"
server_home="/home/liuxu/"
script_dir="/home/liuxu/scripts/"
data_dir="/storage/hpc/data/liuxu/"

src="/mnt/d/RNAseq/"
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
	des=$server_home
fi

echo "uploading ## $src$2 to ## ${server}:$des"
rsync -ahvP $src$2 ${server}:${des}