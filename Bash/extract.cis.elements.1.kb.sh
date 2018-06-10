#Extract
downstream=1000
upstream=1000
cat RefSeq.annotation.noRandom.noChrM.bed | awk '{ if ($4 == "+") { print $1"\t"$2-'$upstream'"\t" $2+'$downstream'"\t"$4"\t"$5} else if ($4 == "-") { print $1"\t"$3+'$upstream'"\t"$3-'$downstream'"\t"$4"\t"$5}}' > RefSeq.10kb.extended.bed
