#!/bin/bash
DGE=$1
RNK=`echo $DGE | sed 's/.xls/.rnk/'`
sed 1d $DGE \
| awk '{OFS="\t"}
{ if ($2>0) printf "%s\t%4.3e\n", $1, 1/$5 ;
else printf "%s\t%4.3e\n", $1, -1/$5 }' \
| sort -k2gr > $RNK