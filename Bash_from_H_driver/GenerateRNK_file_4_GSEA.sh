#!/bin/bash
DGE=$1
RNK=`echo $DGE | sed 's/.xls/.rnk/'`
sed 1d $DGE \
| sort -k7g \
| cut -d '_' -f2- \
| awk '!arr[$1]++' \
| awk '{OFS="\t"}
{ if ($6>0) printf "%s\t%4.3e\n", $1, 1/$7 ;
else printf "%s\t%4.3e\n", $1, -1/$7 }' \
| sort -k2gr > $RNK