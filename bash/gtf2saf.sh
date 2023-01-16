p="../genomes/gtf/"
gtf="Arabidopsis_thaliana.TAIR10.59.gtf"
saf="Arabidopsis_thaliana.TAIR10.59.2.saf"

egrep -h "\sexon\s" $p$gtf | cut -d" " -f1,4,5,7,9 | cut -d";" -f1 | cut -d'"' -f1-2 > $p$saf

