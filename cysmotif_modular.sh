#Written by Andrey Shelenkov, CRIE (c) 2024

#Exemplary commands to reveal multidomain structure in translated transcriptome sequences
#Input file is supposed to contain translated transcriptomic sequences

#Step1
cysmotif_searcher.pl -i AAAA_translated.fasta.bz2 -m motifs_hairpinin.txt -t -b -n 55 -l 600 –c

#Take the motif file from previous run and use it as input; -u options allows to skip output sequence conversion to upper case
#and search for second motif; first motif is skipped by the program since it is in lower case

#Step2
cysmotif_searcher.pl -i AAAA_translated_motifs_orfonly_withM.fasta -m motifs_hairpinin.txt -t -b -n 55 -l 600 -c –u

#Select the output file of step2 for artificial CYSRICH group (the sequences containing additional cysteines after motif go to this group)
#and repeat the procedure
#Step3
cysmotif_searcher.pl -i AAAA_translated_motifs_orfonly_withM_CYSRICH.fasta -m motifs_hairpinin.txt -t -b -n 55 -l 600 -c -u

#Take CYSRICH output file from previous step and continue iterations
#Step4
cysmotif_searcher.pl -i AAAA_translated_motifs_orfonly_withM_CYSRICH_CYSRICH.fasta -m motifs_hairpinin.txt -t -b -n 55 -l 600 -c -u

#Continue until the step on which no CYSRICH sequences will be found

