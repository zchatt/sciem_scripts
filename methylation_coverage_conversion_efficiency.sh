
## sciMET_N9 CG-context methylation coverage of human and mouse aligned reads for CGmap.gz format
## can be modified for other libraries and for CH context

zcat cluster{1..2304}.CGmap.gz |
awk '$1 ~ /^ch/ || $1 ~ /^NC/ {print $0}' |
awk 'BEGIN {total_meth=0;total=0} $4 == "CG" {total_meth+=$7;total+=$8} END {total_unmeth=total-total_meth;print total_meth, total_unmeth, total}'

## conversion efficiency of lambda for true positive controls only
## this is just an example

zcat cluster2302.CGmap.gz cluster4606.CGmap.gz |
awk '$1 == "Lambda_NEB" {print $0}' |
awk 'BEGIN {total_meth=0;total=0} {total_meth+=$7;total+=$8} END {meth_percent=total_meth/total;print total_meth, total, meth_percent}'
