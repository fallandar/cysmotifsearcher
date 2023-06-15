# cysmotifsearcher
Cysmotif searcher

Written by Andrey Shelenkov, Central Reseach Institute of Epidemiology (www.crie.ru).<br>
Last update to the program package was made on April 12, 2023. Current version is 3.3.3.<br>Made avaialble in github on November 12, 2017.
<p>
 data1K folder contains the results of pipeline running on transcriptomes from 1kP project (https://sites.google.com/a/ualberta.ca/onekp/ , see https://doi.org/10.3390/antibiotics9020060 for the paper describing this analysis)
<p>
Cysmotif searcher is a set of Perl scripts that performs profile search to reveal peptide sequences possessing cysteine motifs common to various families of AMPs and other cysteine-rich peptides. It can be run on any Linux machine with Perl installed, and can also be executed on Windows machines, but with some limitations (SPADA and SignalP cannot be integrated into the pipeline in this case). Motifs to search for were derived from literature, and then were supplemented and further refined during several iterations of searching and refining steps.
<p>
Scientific papers describing the Cysmotif searcher pipeline and its applications to real biological data analysis can be found using the links below:<br>
 <ol><li>A. Shelenkov, A. Slavokhotova, T. Odintsova (2020) Predicting Antimicrobial and Other Cysteine-Rich Peptides in 1267 Plant Transcriptomes, Antibiotics (Basel), 9(2):60, https://doi.org/10.3390/antibiotics9020060, https://www.ncbi.nlm.nih.gov/pubmed/32032999 (please cite this paper if you use the results provided in data1k folder)</li>
 <li>A.A. Shelenkov, A.A. Slavokhotova, and T. I. Odintsova (2018) Cysmotif Searcher Pipeline for Antimicrobial Peptide Identification in Plant Transcriptomes, Biochemistry (Moscow), 83(11), 1424-1432, https://www.ncbi.nlm.nih.gov/pubmed/30482154 (please cite this paper if you publish the results obtained using Cymotif searcher for your own data)</li>
 <li>A.A. Slavokhotova, A. A. Shelenkov et al. (2017) Defense peptide repertoire of Stellaria media predicted by high throughput next generation sequencing, Biochimie, 135:15-27, https://www.ncbi.nlm.nih.gov/pubmed/28038935 </li>
 <li>A. A. Slavokhotova, A. A. Shelenkov, T. I. Odintsova (2015) Prediction of Leymus arenarius (L.) antimicrobial peptides based on de novo transcriptome assembly, Plant Mol. Biol., 89(3):203-14, https://www.ncbi.nlm.nih.gov/pubmed/26369913</li></ol><p><p>

General usage: cysmotif_searcher.pl [OPTIONS] -m motifs.txt -i input.fasta<br>
You can always get help by typing cysmotif_searcher.pl --help<p>

List of options is as follows:<br>
-i FILE1                set input fasta file to FILE1 (REQUIRED)<br>
-m FILE2                set input motif file to FILE2 (REQUIRED)<p>

-b                      delete everything after first whitespace or star (*) in input fasta headers (useful for trinity)    (default=not active)<br>
-c                      print motifs in CYSRICH group in lowercase (default=off)<br>
-k NUM                  remove output sequences that are a subset of other sequences (between spada and cysmotif, NUM is max lg diff for seqs)(default=keep all)<br>
-f                      print results for motifs only to one file (default=each motif to separate file)<br>
-g                      start with using signalP (input file supposed to be like orfonly_with_M)<br>
-l LG                   set max length for mature peptide to LG (default=150)<br>
-k NUM                  remove output sequences that are a subset of other sequences (between spada and cysmotif, NUM is max lg diff for seqs)(default=keep all)<br>
-n NUM                  set number of threads to NUM (default=1)<br>
-p,--prefix PREF        set output file prefix to PREF (default=input filename without extension)<br>
-s                      run spada to get additional cys-rich seqs (should be installed separately, also needs signalP) (default=not active)<br>
-t                      skip translation of input sequences (input file must contain amino acid seqs) (default=not active)<br>
-y DIR                  get spada results from DIR (already calculated), -s option is set automatically<br>
-z FILE                 set nucleotide fasta file to FILE for spada input (needed when translation is skipped in main program), -s option is set automatically<br>
-h, -help, --help        displays help and nothing else<br>
-v, --version            shows current version number of the program<br>
<p>
Requirements:

In order to run the package, yoy should have Perl 5.8 or later installed on your machine (https://www.perl.org/get.html).

In order to check for presence of signal peptides in motifs (which is very important step of filtration) you should have SignalP program (http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp) installed and should specify the path to its executable in the script under $signalp variable at the beginning of the script.

In order to include SPADA in the computational pipeline, you should download spada from https://github.com/orionzhou/SPADA and follow the installation procedure described there. Then you should specify the path to SPADA installation in cysmotif_searcher.pl under $spada_dir variable and to SPADA executable (spada.pl) under $spada variable.
<p>
 To test your installation of cysmotifsearcher you can use the demo input file demo.fasta. Run the following command:
 <br>
 <b>cysmotif_searcher.pl -i demo.fasta -m motifs.txt</b>
 <br>
 If you do not have any results (motifs) in the output files, then something is wrong with the installation.<br> If you do not have any results for your own input files, this could simply mean that cysmotif searcher has not found any motifs in your data. You can try additional searching with SPADA option turned on.
 <p>
 Licensed under the GNU GENERAL PUBLIC LICENSE Version 3 (the "License");
 you may not use this file except in compliance with the License. You may obtain a copy of the License at https://www.gnu.org/licenses/gpl.html
