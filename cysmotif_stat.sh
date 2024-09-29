#Written by Andrey Shelenkov, CRIE (c) 2024
#Call is: cysmotif_stat.sh output_prefix input_fasta
#Output prefix will be used in all output files; 
#input file is a combined fasta output of cysmotif_searcher.pl
#1kP-Sample-List_2019-03.txt is available in 'data1k' folder on Github
#(https://github.com/fallandar/cysmotifsearcher)

nm=$1

sort -k1 1kP-Sample-List_2019-03.txt > 1kP-Sample-List_2019-03_sorted.txt

grep ">" $2  | awk -F'[-]' '{print $2}' | uniq -c | sed -e 's/^[ ]*//' | awk ' { t = $1; $1 = $2; $2 = t; print; } '| sed -e 's/[ ]/\t/' > ${nm}_motifs_tab.stat

join -j 1 ${nm}_motifs_tab.stat ../1kP-Sample-List_2019-03_sorted.txt -t $'\t' > ${nm}_motifs_tab1.stat
awk -F $'\t' '{a[$3] += $2} END{for (i in a) printf "%s\t%s\n", i, a[i]}' ${nm}_motifs_tab1.stat | sort -t$'\t' -k2 -nr > ${nm}_motifs_class1.stat
awk -F $'\t' '{a[$4] += $2} END{for (i in a) printf "%s\t%s\n", i, a[i]}' ${nm}_motifs_tab1.stat | sort -t$'\t' -k2 -nr > ${nm}_motifs_class2_family.stat
awk -F $'\t' '{a[$6] += $2} END{for (i in a) printf "%s\t%s\n", i, a[i]}' ${nm}_motifs_tab1.stat | sort -t$'\t' -k2 -nr > ${nm}_motifs_classx_parts.stat
grep "cformula" $2 | cut -f2,3,4 -d'*' | sed 's/cformula://' | sed  's/*/ /g' |  sort | uniq -c | sed 's/ \+/ /' |  tr ' ' '\t' | sort -k2,2 -k1,1nr | cut -c2- | awk -F'\t' 'BEGIN {OFS = FS} { print $2, $3, $4, $1 }' > ${nm}_cformula.stat
grep ">" $2  | cut -f1,4 -d '*' | sed 's/cformula://' | sed 's/\*/\t/' >${nm}_cformula.txt
sed -i 's/\"//' ${nm}_motifs_classx_parts.stat


#Alternative statistics collection for several output folders at once to collect Cys formula info
#cat *results/*final.fasta > ${nm}_all.fasta &&\
#get_fasta_from_list.pl -e -v -k 0 ${nm}_all.fasta > ${nm}_all_NOTCYSRICH.fasta &&\
#grep "cformula" ${nm}_all_NOTCYSRICH.fasta | cut -f2,3,4 -d'*' | sed 's/cformula://' | sed  's/*/ /g' |  sort | uniq -c | sed 's/ \+/ /' |  tr ' ' '\t' | sort -k2,2 -k1,1nr | cut -c2- | awk -F'\t' 'BEGIN {OFS = FS} { print $2, $3, $4, $1 }' > ${nm}_all_NOTCYSRICH_cysformula.stat

