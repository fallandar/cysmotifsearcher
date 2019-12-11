#!/usr/bin/perl
use File::Copy 'move';
use threads;
use Config;
use FindBin; 
use lib "$FindBin::Bin";
use cysmotif;

$program="cysmotif_searcher.pl";  
$version="3.3.1";
$last_update="December 11, 2019";
$comment="Written by Andrew Shelenkov, VIGG of RAS";

#$signalp="/home/fallandar/spada_soft/signalp-4.1/signalp";
#$signalp="/mss2/export/Projects/Chrysanthemum-miRNA/raw_data/signalp-4.1/signalp";
#$signalp="/data2/andrew/bin/signalp-4.1/signalp";
$signalp="/export/home/shelenkov/bin/signalp-4.1/signalp";
$signalpx="/export/home/shelenkov/bin/signalp-5.0/bin/signalp";
$spada="/export/home/shelenkov/soft/spada/spada.pl";
$spada_dir="/export/home/shelenkov/soft/spada";

@ids=();
@seqs=();

sub format_pattern
{
	my $out = shift;
	$out=~s/\[ABD-Z\]/X/g;
	$out=~s/\?\^\:\(//;
	$out=~s/[\(\)]//g;
	$out=~s/\?\-xism\://;
	return $out;
}#format_pattern

sub translate6
{
	my ($start,$end)=@_;
	my ($name,$seq,$protein,$codon,$i,$j,$n,$p);
	my $toprint="";
	
	for ($i=$start; $i<$end;$i++)
	{
 	    $name=$ids[$i];
   		$seq=$seqs[$i];
       $p=1;
       while ($p>=-1) 
       {
       	#translate seq in 6 frames and print results to output file
         $n=1;
         for ($n=1;$n<=3;$n++)
         {
         	$protein='';
           for($j=$n-1;$j<(length($seq)-2);$j+=3)
           {
           $codon=substr($seq,$j,3);
           $protein.=&codon2aa($codon);
           }
           $toprint.= sprintf "%s %d\n%s\n",$name,$n*$p,$protein;
         }
         $p-=2;
         $seq=revcomp($seq);
       }
  }
 return $toprint;
}#translate6

sub print_log
{
 my $datestring=localtime();
 printf LOG "[$datestring] %s",$_[0];
 unless ($_[1])
 {printf LOG "\n";}
 return 0;
}#print_log

sub cformula
{
	my $t2x=lc $_[0];
  my $numc = $t2x =~ tr/c/x/;
  my $offset = 0; my @posc=();
  while (1)
  {
  	my $position = index($t2x, "x", $offset);
  	last if ($position < 0);
  	push @posc,$position;
  	$offset=$position+1;
  }
  my $formc="";
  for (my $i=0;$i<=$numc-2;$i++)
  {
  	my $vl=-1-$posc[$i]+$posc[$i+1];
  	if ($vl>0)
  	{$formc.=sprintf "CX\{%d\}",$vl;}
  	else
  	{$formc.="C";}
  }
  $formc.="C";
  return $formc;
}#cformula

#default values and init
$delete_after_blank=0;
#$noendc=1;
$no_spada=1;
$num_threads=1;
$prefix="";
$print_files=1;
$skip_translation=0;
$spada_nucl="";
$spada_precalc_dir="";
$start_with_sp=0;
$keep_subseqs=1; $keep_subseqs_lg=1;
%origseq=();
%tabseq=();

@arr=localtime;
$dat=sprintf '%04d%02d%02d', $arr[5]+1900, $arr[4]+1, $arr[3]; 
$print_header = sprintf "#Made with %s, version %s, started on %s#\n# %s #\n\n",$program,$version,$dat,$comment; 
$start = time;
$n01=0; $posn01=-1;
$limit_length=150;

if (($ARGV[0] eq "--version") || ($ARGV[0] eq "-v") )
{ #print the number of version and date of release
  print "You are using $program version $version\n";
  print "Last update: $last_update\n";
  print "$comment\n";
  exit;
} 
if ($ARGV[0] eq "-h" || $ARGV[0] eq "-help" || $ARGV[0] eq "--help")
{
  print "You are using $program version $version \n";
  print "$comment\n\n";
  printf "Usage1: $program [OPTIONS] -m motif.txt -i input.fasta\n";
  printf "Usage2: $program -g -i input_m.fasta\n";
  printf "-i FILE1                set input fasta file to FILE1 (REQUIRED)\n";  
  printf "-m FILE2                set input motif file to FILE2 (REQUIRED)\n\n";  
  printf "-b                      delete everything after first whitespace or star (*) in input fasta headers (useful for trinity) (default=not active)\n";
  printf "-k NUM                  remove output sequences that are a subset of other sequences (between spada and cysmotif, NUM is max lg diff for seqs)(default=keep all)\n";
  printf "-f                      print results for motifs only to one file (default=each motif to separate file)\n";
  printf "-g                      start with using signalP (input file supposed to be like orfonly_with_M)\n";  
  printf "-l LG                   set max length for mature peptide to LG (default=150)\n";  
  printf "-n NUM                  set number of threads to NUM (default=1)\n";  
  printf "-p,--prefix PREF        set output file prefix to PREF (default=input filename without extension)\n";
  printf "-s                      run spada to get additional cys-rich seqs (should be installed separately, also needs signalP) (default=not active)\n";  
  printf "-t                      skip translation of input sequences (input file must contain amino acid seqs) (default=not active)\n";
  printf "-y DIR                  get spada results from DIR (already calculated), -s option is set automatically\n";
  printf "-z FILE                 set nucleotide fasta file to FILE for spada input (needed when translation is skipped in main program), -s option is set automatically\n";
  print "-h, -help, --help        displays help and nothing else\n";
  print "-v, --version            shows current version number of the program\n";
  exit;
}
 $arg=join " ",@ARGV;
 $nfile=0; 
while ($ARGV[0]=~ /^-/)
 #reading options and keys
{
 if ($ARGV[0]=~ /-b/)
 #delete everything after first whitespace in input fasta headers
 {$delete_after_blank=1;}
 #if ($ARGV[0]=~ /-c/)
 #allow presence of cysteines in mature peptide after motif (default=not allowed)
 #{$noendc=0;}
 elsif ($ARGV[0]=~ /-f/)
 #print results for each motif to separate file
 {$print_files=0;}
 elsif ($ARGV[0]=~ /-g/)
 #start with signalp; skip all previous steps
  {$start_with_sp=1; $skip_translation=1;}
 elsif ($ARGV[0]=~ /-i/)
 #set input file
 {
 	shift @ARGV; 
 	$filename_source=$ARGV[0];
 	die "ERROR: Input file (-i arg) $filename_source not found: $!\n" unless (-e $filename_source);
 	$nfile++;
 }
 elsif ($ARGV[0]=~ /-k/)
 #remove spada seqs that are subseq of cysmotif (or vice versa)
 {
 	$keep_subseqs=0;
 	shift @ARGV; 
 	$keep_subseqs_lg=$ARGV[0];
 	die "ERROR: Value for -k option should be >=1! ($keep_subseqs_lg)\n" unless ($keep_subseqs_lg>=1);
 }
 elsif ($ARGV[0]=~ /-l/)
 #set max length for mature peptide
 {
 	shift @ARGV; 
 	$limit_length=$ARGV[0];
 	die "ERROR: Wrong length is set with -l arg ($limit_length). Must be >=50.\n" unless ($limit_length >=50);
 }
 elsif ($ARGV[0]=~ /-m/)
 #set motif file
 {
 	shift @ARGV; 
 	$motif_file=$ARGV[0];
 	die "ERROR: Motif file (-m arg) $motif_file not found: $!\n" unless (-e $motif_file);
 	$nfile++;
 }
  elsif ($ARGV[0]=~ /-n/)
 #set number of threads
 {
 	shift @ARGV; 
 	$num_threads=$ARGV[0];
 	die "ERROR: Wrong number of threads. Must be >=1\n" unless ($num_threads>=1);
 }
 elsif ($ARGV[0]=~ /-p/ || index($ARGV[0],"--prefix")==0)
 #set output filename prefix
 {
 	shift @ARGV; 
 	$prefix=$ARGV[0];
 }
 elsif ($ARGV[0]=~ /-s/)
 #run spada
 {$no_spada=0;}
 elsif ($ARGV[0]=~ /-t/)
 #skip translation
 {
 	$skip_translation=1;
 	$no_spada=1; #spada need nucleotide seqs
 }
 elsif ($ARGV[0]=~ /-y/)
 #set dir with spada results (if already present)
 {
 	shift @ARGV;
 	$spada_precalc_dir=$ARGV[0];
 	die "ERROR: Spada results dir (-y arg) $spada_precalc_dir not found: $!\nCheck dir or run spada once more.\n" unless (-d $spada_precalc_dir);
 	$no_spada=0;
 }
 elsif ($ARGV[0]=~ /-z/)
 #set nucl file for spada
 {
 	shift @ARGV;
 	$spada_nucl=$ARGV[0];
 	die "ERROR: Input nucleotide spada file (-z arg) $spada_nucl not found: $!\n" unless (-e $spada_nucl);
 	$no_spada=0;
 }
 else 
 {
 	die "ERROR: Unrecognized input option: $ARGV[0]. Processing stopped. Check $program --help for available options.\n";
 }
	 shift @ARGV;
}
die "ERROR: Not enough input files specified!\nUsage: $program [OPTIONS] -m motif.txt -i input.fasta\nRun '$program -h' for additional info.\n" if (($nfile<2)&&(!$start_with_sp));

$statn="";
unless ($prefix)
{$filename_noext=$filename_source; $filename_noext=~s/\.bz2$//;$filename_noext=~s/\.gz$//; $filename_noext=~ s/\.[^.]+$//; $prefix=$filename_noext; $statn="${filename_noext}_stat.txt"}
else
{$filename_noext=$prefix; $statn="${filename_noext}_stat.txt";}

#create names
%names=(mpre=>"${filename_noext}_motifs_pre.fasta",morf=>"${filename_noext}_motifs_orfonly.fasta", morfM=>"${filename_noext}_motifs_orfonly_withM.fasta",
     mhtml=>"${filename_noext}_motifs.html",gff=>"${filename_noext}_motifs_orfonly_withM_signalp.gff",tabbed=>"${filename_noext}_motifs_final_tabbed.txt",
     final=>"${filename_noext}_motifs_final.fasta",logf=>"${filename_noext}.log",stats=>$statn);
     
unless ($no_spada) 
 {
 	$names{withspada}="${filename_noext}_motifs_final_with_spada_uniq.fasta";
 	$names{spada}="${filename_noext}_spada.fasta";
 	$names{spadalog}="${filename_noext}_spada.log";
 	$names{spadatabbed}="${filename_noext}_motifs_all_tabbed.txt";
 }
open (LOG,">$names{logf}") || die "Can not open log file for writing: $!\n";
#flush LOG
my $previous_default = select(STDOUT);
select(LOG);
$|++;
select($previous_default); 
$|=1;

printf LOG "###Starting %s, version %s, last update %s###\n### %s ###\n",$program,$version,$last_update,$comment; 
print_log "Command line: $arg\n";

if (-e $signalp)
{
	#get version of signalp
	$signalp_version=4;
	open (INN, "$signalp -version 2>/dev/null|");
	$line=<INN>;
	close IN;
	if (index($line,"version 5.0")!=-1)
	{$signalp_version=5;}
}

unless ($start_with_sp)
{
open (IN,$motif_file) || die "Can not open $motif_file: $!\n";
@patterns=();
@patterns_nm=();
%patterns_nm_hash=();
$num_motifs=0;

  print_log ("Reading motif file $motif_file... ",1);
  #read motifs to be searched from file
  while (<IN>)
  {
   unless (/^\#/)
   #skip commented motifs
   {
   	next if (index($_,"C")==-1);
   	chomp;
   	s/ //g;
   	s/X([0-9]{1,2})\-([0-9]{1,2})/\[ABD\-Z\]\{$1,$2\}/g;
    s/X([0-9]{1,2})/\[ABD\-Z\]\{$1\}/g;
    s/X/\[ABD\-Z\]/g;
    
    if (/\t/)
    {
     #new format
     @arr=split /\t/,$_;
     $patterns[$num_motifs]=qr/($arr[1])/;
     $patterns_nm[$num_motifs]=$arr[0];
     if ($arr[0] eq "N01") {$n01=1; $posn01=$num_motifs;}
     $patterns_nm_hash{$arr[0]}=$num_motifs;
    }
    else
    {
    	#old format
     $patterns[$num_motifs]=qr/($_)/;
     $tmp=sprintf "motif%02d",$num_motifs+1;
     $patterns_nm[$num_motifs]=$tmp;
     $patterns_nm_hash{$tmp}=$num_motifs;
     #print "$num_motifs $patterns_nm[$num_motifs]\n";
    }
   	$num_motifs++;
   }	
  }
  close IN;
  print LOG "done. $num_motifs motifs read.\n";
}
$bignum=$num_motifs+10;
@tarr=(0..$num_motifs-1,$bignum-1);
$patterns_nm[$bignum-1]="CYSRICH"; $patterns_nm_hash{"CYSRICH"}=$bignum-1;
#read contigs and translate by 6 reading frames
unless ($skip_translation)
{
 print_log "Starting translation of input seqs in 6 reading frames...";
 $names{translated}=$filename_noext."_translated.fasta";
 if ($filename_source=~/\.gz$/)
 {open ($fh,"zcat $filename_source |") || die "Can not open $filename_source: $!\n";}
 elsif ($filename_source=~/\.bz2$/)
 {open ($fh,"pbzip2 -c -d $filename_source |") || die "Can not open $filename_source: $!\n";}
 else
 { open ($fh,$filename_source) || die "Can not open $filename_source: $!\n";}
 
while (read_fasta_sequence($fh, \%sequence_data_tmp)) 
{
	$id=$sequence_data_tmp{header};
	$seq=$sequence_data_tmp{seq}; 
	
	if ($delete_after_blank)
	{
	 $posb=index($id," ");
	 if ($posb!=-1)
	 {$id=substr($id,0,$posb);}
  }
	$seq=~s/\s+//g;
  $seq=~s/\*//g;
  $seq=~s/[^A-Za-z]//g;
  $seq2=uc ($seq);
  
  push @ids,$id;
  push @seqs,$seq2;
}
close $fh;
 unless (@ids == @seqs)
 {
 	print_log "Corrupted input fasta file!";
  die "Corrupted input fasta file!\n" ;
 }
unless ($Config{useithreads})
#threads are not supported
{$num_threads=1; print_log "Threads are not supported. Working in 1 thread.";}

open (OUT,">$names{translated}") || die "Can not open $names{translated} for writing:$!\n";

$n=@ids;
if ($num_threads > 1)
{
	#run translation in threads
	@arr_comb=();
  $j=int($n/($num_threads)+0.5);
	for ($i=0; $i<$n; $i+=$j)
	{push @arr_comb, $i;}
  push @arr_comb, $n;
  $sz=@arr_comb;
	
  @threads=();
  for ($i=0;$i<=$sz-2;$i++)
  {$threads[$i]=threads->create(\&translate6,$arr_comb[$i],$arr_comb[$i+1]);}
  for ($i=0;$i<=$sz-2;$i++)
  {print OUT $threads[$i]->join;}
}
else
{
	#no threads
	print OUT translate6(0,$n);
}
close OUT;
print_log "Translation finished.";
}#tr
else
{	
	$names{translated}=$filename_source;
	print_log "Translation skipped.";
}

unless ($start_with_sp)
{
	print_log "Searching for motifs...";
  #open files for motif printing
  if ($print_files)
  {
  	%filehandle=();
    for $i (@tarr)
    {
    	$name = sprintf "%s_%s.fasta",$filename_noext,$patterns_nm[$i];
  	  open $filehandle{$i}, ">", $name or die "Can't open $name for writing: $!";
    }
  }
  if ($names{translated}=~/\.gz$/)
   {open ($fh,"zcat $names{translated} |") || die "Can not open $names{translated}: $!\n";}
  elsif ($names{translated}=~/\.bz2$/)
   {open ($fh,"pbzip2 -c -d $names{translated} |") || die "Can not open $names{translated}: $!\n";}
  else
  { open ($fh,$names{translated}) || die "Can not open $names{translated}: $!\n";}
  
  open (OUT,">$names{mpre}") || die "Can not open $names{mpre} for writing: $!\n";
  open (OUTX,">$names{morf}") || die "Can not open $names{morf} for writing: $!\n";
  open (OUTXX,">$names{morfM}") || die "Can not open $names{morfM} for writing: $!\n";
  open (OUTH,">$names{mhtml}") || die "Can not open $names{mhtml} for writing: $!\n";
  
  printf OUTH "<HTML>\n<BODY>\n";
  print OUTH "$print_header<p>";
  @max=(0, 0, 0);
  @min=(10000,10000,10000);
  %uniq=();
  $num_seqs=0;
  
  while (read_fasta_sequence($fh, \%sequence_data_tmp)) 
  {
  	    $name=$sequence_data_tmp{header};
  	    $seq=uc($sequence_data_tmp{seq}); 
  	    if ($skip_translation && !$num_seqs && $seq=~/^[ACTG]+$/)
  	    {
  	    	#nucleic acid sequence found in translated file - exiting
  	    	  foreach $file (keys (%names))
            {
            	if ($file ne "logf" && $file ne "translated")
            	{
            		unlink $names{$file};
            	}
            }
            if ($print_files)
            {
              for $i (@tarr)
              {
            	  close $filehandle{$i};
              }
              for $i (@tarr)
              {
              		$name = sprintf "%s_%s.fasta",$filename_noext,$patterns_nm[$i];
              		if (-z $name)
                  #delete empty motif files
              		{unlink "$name";}
              		else
              		{$names{$name}=$name;}
              }
            }
          print LOG "Sequences in $names{translated} seem not to be amino acid sequences! Check your input file!\n";
          close LOG;  
  	    	die "Sequences in $names{translated} seem not to be amino acid sequences! Check your input file!\n";
  	    }
  	    	if ($delete_after_blank)
        	{
	         $posb=index($name," ");
	         if ($posb!=-1)
	         {$name=substr($name,0,$posb);}
	         $posb=index($name,"*");
	         if ($posb!=-1)
	         {$name=substr($name,0,$posb);}
          }
	      $name=~s/\*/_/g;
    	  $k=0;
    	  $match=""; $matchx=""; $match_orf=""; $numm=-1; $final_endc=0;
    		for ($i=0;$i<=$num_motifs-1;$i++) 
    		{
    			#search for motifs in sequence
          if( $seq=~/$patterns[$i]/ ) 
          {
   		      $tmp=$1;
   		      	if ($i == $posn01)
   		       	  {$num1=length $2; $num2=length $3; $num3=length $4;}
   		       	  
   		      if ( ($seq=~/_([A-Z]*${tmp}([A-Z]*))_/) || ($seq=~/_([A-Z]*${tmp}([A-Z]*)$)/)|| ($seq=~/(^[A-Z]*${tmp}([A-Z]*))_/) || ($seq=~/(^[A-Z]*${tmp}([A-Z]*)$)/)) 
     		     {
     		     	$endc=0;
     		     	#collect uniq seqs in ORF
     		     	$tmporf=$1;
     		     	$seqafter=$2;
     		     	#cysteine found after motif
     		     	if (index($seqafter,"C")!=-1)
     		     	{$endc=1;}
     		     	unless ($k)
     		     	#first motif found
     		     	{	$match=$tmp;$match_orf=$tmporf;$numm=$i;$final_endc=$endc;}	
     		    	else
     		    	{
     		    		$matchx=$tmp; $x=0; $y=0;
     		    		$x++ while ($match =~ m/C/g);
     		    		$y++ while ($matchx =~ m/C/g);
     		    		#new motif from this orf has more cysteins in it
     		    		if ($x<$y && !$endc)
     		    		{$match=$matchx;$match_orf=$tmporf;$numm=$i;$num1=0;$num2=0;$num3=0;$final_endc=0;}
     		      }
     		     }
     		    $k++ unless ($endc); 
          }
        }#motifs
       if ($numm+1)
       {
       	#at least one motif found 
        if ($n01 && $numm == $posn01)
        {
           #get consensus for 1st motif	
           if ($max[0]<$num1) {$max[0]=$num1;} if ($min[0]>$num1) {$min[0]=$num1;}
   		     if ($max[1]<$num2) {$max[1]=$num2;} if ($min[1]>$num2) {$min[1]=$num2;}
     		   if ($max[2]<$num3) {$max[2]=$num3;} if ($min[2]>$num3) {$min[2]=$num3;}
     		   @tofile=(sprintf ("%s*%s*%s,pattern:CX{%d}CXXXCX{%d}CX{%d}CXC\n", $name,$patterns_nm[$numm],format_pattern($patterns[$numm]),$num1,$num2,$num3),$match,$seq,$numm+1);
     		}
   		   else
   		  {
   		  	unless ($final_endc)
   		  	 {@tofile=(sprintf ("%s*%s*%s\n", $name,$patterns_nm[$numm],format_pattern($patterns[$numm])),$match,$seq,$numm+1);}
   		  	else
   		  	#cys after motif - CYSRICH
   		  	{@tofile=(sprintf ("%s*%s*%s\n", $name,"CYSRICH",""),$match,$seq,$bignum);}
   		  }
   		   $uniq{$match_orf}=[@tofile];
   		}
      	if (($num_seqs+1) % 50000 == 0)
      	{printf LOG "%d sequences processed for motifs\n",$num_seqs+1;} 
      	$num_seqs++;
      	
  }#IN
  close $fh;
  print_log "Initital search finished. Printing results...";
  
       %motif_stat=(); %motif_stat_m=(); %motif_stat_final=();
       $motif_all=0; $motif_all_m=0; $motif_cysrich=0;
       foreach $key (keys (%uniq))
       {
	
       	$tmp=$uniq{$key}[2];
       	$tmp2=$uniq{$key}[1];
       	$tmp3=$key;
       	unless ($uniq{$key}[3]==$bignum)
       	{
       	 $tmp2x=lc($tmp2);
       	 $tmp=~s/$tmp2/$tmp2x/;
       	 $tmp3=~s/$tmp2/$tmp2x/;
       	}
       	else
       	#CYSRICH - no motif in lower case 
       	{$tmp2x=$tmp2;}
       	#print full seqs with motifs in lower case
       	printf OUT "%s%s\n",$uniq{$key}[0],$tmp;
      	#print orfs only
       	printf OUTX "%s%s\n",$uniq{$key}[0],$tmp3;
        if ($tmp3=~/(M.*${tmp2x}.*$)/)
       	#there exists methionine before motif
       	{
       		$motif_stat_m{$uniq{$key}[3]}++; $motif_all_m++;
       		printf OUTXX "%s%s\n",$uniq{$key}[0],$1;
       	}
       	
       	$motif_stat{$uniq{$key}[3]}++;$motif_all++;
       		if ($print_files)
  	      {
  	      	printf {$filehandle{$uniq{$key}[3]-1}} "%s%s\n",$uniq{$key}[0],$tmp;
  	      }
  	      
        $tmp=~s/[c]/<font color=red>c<\/font>/g;
   		  $tmp=~s/[C]/<font color=red>C<\/font>/g;
   		  printf OUTH "<BR>%s\n<BR>%s\n<BR>",$uniq{$key}[0],$tmp;
       }
  close OUT;
  close OUTX;
  close OUTXX;
  printf OUTH "<BR>;full range pattern: CX{%d,%d}CXXXCX{%d,%d}CX{%d,%d}CXC\n\n",$min[0],$max[0],$min[1],$max[1],$min[2],$max[2] if ($n01 && $max[0]+$max[1]+$max[2]>0);
  printf OUTH "<\/BODY>\n<\/HTML>\n";
  close OUTH;
  
  if ($print_files)
  {
    for $i (@tarr)
    {
  	  close $filehandle{$i};
    }
    for $i (@tarr)
    {
    		$name = sprintf "%s_%s.fasta",$filename_noext,$patterns_nm[$i];
    		if (-z $name)
        #delete empty motif files
    		{unlink "$name";}
    		else
    		{$names{$name}=$name;}
    }
  }
}
else
{
	print_log "Starting processing from signalP analysis...";
	$names{morfM}=$filename_source;
}   
if (-e $signalp)
#run signalp
{
	if ($signalp_version!=5)
	{
		$names{gff}=$filename_noext."_signalp.gff";
   	#make compatible file for stupid signalp
   	open (XX,$names{morfM});
   	open (XXO,">${filename_noext}_signalp.tmp");
   	while (<XX>)
   	{
   		if (/>/)
   		{
   			chomp;
   			@arr=split /\*/;
   			@arr2=split /\-/,$arr[0];
   			printf XXO "%s\n",$arr2[0]."-".$arr2[1]."-".$arr2[2];
   			#printf XXO "%s\n",$arr[0];
   		}
   		else
   		{printf XXO $_;}
   	}
   	close XX;
   	close XXO;
   	#run
    $cmd=sprintf "%s -c 70 -n %s %s >/dev/null 2>&1",$signalp,$names{gff},"${filename_noext}_signalp.tmp";
  }
  else
  {
  	($names{gff}=$names{morfM})=~s/.fasta$/.gff3/;
  	($pmt=$names{morfM})=~s/.fasta$/_summary.signalp5/;
  	$cmd=sprintf "%s -gff3 -org euk -fasta %s >/dev/null 2>&1",$signalp,$names{morfM};
  	
  }
 print_log "Running signalP as ".$cmd;
 system $cmd || print STDERR "signalp failed: $!\n";
 if ($signalp_version!=5)
  {unlink "${filename_noext}_signalp.tmp" ;}
 else 
  {unlink $pmt;}
}
else 
{	print_log "WARNING: signalp was not found in $signalp. Skipping.";}

if (-e $names{gff})
#processing signalp results
{
	print_log "Processing signalP results - reading ".$names{gff};
	open IN,$names{gff};
   %sig=();
   while (<IN>)
   {
   	unless (/^\#\#/)
   	{
   		@arr=split /\t/;
   		$sig{$arr[0]}=$arr[4];
   	}
   }
   close IN;
   sort keys(%sig);
   open ($fh, $names{morfM});
   
   unless ($start_with_sp)
   {$names{end}=sprintf "%s_motifs_orfonly_withM_signalp_uniq_lt%d.fasta",$filename_noext,$limit_length; }
   else
   {$names{end}=sprintf "%s_signalp_uniq_lt%d.fasta",$filename_noext,$limit_length;}
   open OUT3, ">$names{end}" || die "Can not open $names{end} for writing: $!\n";
   
   open OUT4, ">$names{tabbed}" || die "Can not open $names{tabbed} for writing: $!\n";
   printf OUT4 "\t\tseq\tmotif_name\tseq_name\tCformula\tmotif\n";
      
   unless ($start_with_sp)
   {$names{signalp}=$filename_noext."_motifs_orfonly_withM_signalp.fasta";}
   else
   {$names{signalp}=$filename_noext."_signalp.fasta";}
   
   open OUT2, ">$names{signalp}" || die "Can not open $names{signalp} for writing: $!\n";
   
   $motif_sigp=0; %uniq_signalp=(); 
   $motif_sigp_cysrich=0;
   if ($start_with_sp) {$num_sp=0;}
   while (read_fasta_sequence($fh, \%sequence_data_tmp)) 
  {
  	$name=$sequence_data_tmp{header};
	  $seq=$sequence_data_tmp{seq};
   	
   	if ($start_with_sp) {$num_sp++;}
   	$flag=0;
   	foreach $key (sort keys (%sig))
   	{
   		if (index($name,"$key")!=-1)
   		{
   			#seq belongs to set of seqs found in gff file earlier
   			$flag=1;
   			$val=$sig{$key};
   			last;
   		}
   	}
   	if ($flag==1)	
   	 {
   	 	 #get splitting point
   			$t1=substr($seq,0,$val);
   			$t2=substr($seq,$val);
   			$cysrich=0;
   			if ($t1!~/[a-z]/ && $t2!~/[a-z]/)
   			#CYSRICH
   			{$cysrich=1;}
   			unless ($t1=~/[a-z]/)
   			#not inside motif
   			{
   				printf OUT2 "%s\n%s %s\n",$name,$t1,$t2; 
   				if ($cysrich)
   				{$motif_sigp_cysrich++;}
   				else
   				{$motif_sigp++;}
   				
   				unless ((length($t2)>$limit_length) || (exists $uniq_signalp{$seq}))
   				{ 
   				 $uniq_signalp{$seq}=$name; 
   				 #get cysteine formula for mature peptide
   				 $formc=cformula($t2);
   				 #print to combined final fasta file
   				 $nmt=sprintf "%s*cformula:%s",$name,$formc; 
   				 $seqt=sprintf "%s %s",$t1,$t2;
   				 print OUT3 "$nmt\n$seqt\n";
   				 $origseq{$nmt}=$seqt;
   				 $arr[0]=~s/>//;
   				 #print to tab-delimited output file
   				 $topr=sprintf "\t\t%s %s\t%s\t%s\t%s\t%s\n",$t1,$t2,$arr[1],$arr[0],$formc,$arr[2];
   				 print OUT4 $topr;
   				 $tabseq{$nmt}=$topr;
   				}
   			}
   	 }
   }
   close $fh;
   close OUT2;
   close OUT3;
   close OUT4;
}
else 
{
	print_log "WARNING: signalp results were not found!"; 
	$names{end}=$filename_noext."_motifs_orfonly_withM.fasta";
	delete $names{gff};
}

#remove seqs those are subseqs of other seqs (for cysmotif only now)
  $cmd=sprintf "cysmotif_fasta_cleaner.pl -l %d %s > %s", $keep_subseqs_lg,$names{end},"${filename_noext}_tmp1.fasta";
 	system $cmd;
 	open $fh,"${filename_noext}_tmp1.fasta" || print LOG "WARNING: Can not open filtered motif file - skipping it: $!\n";
 	open (OUT2,">$names{final}") || die "Can not open $names{final} for writing: $!\n"; 
 	%motif_stat_final=(); $motif_sigp_un_lt_cysrich=0; $motif_sigp_un_lt=0;
  while (read_fasta_sequence($fh, \%sequence_data_tmp)) 
  {
  	#read sequences from filtered file and restore splitting points and motifs from %origseq
  	$name=$sequence_data_tmp{header};
  	#get final statistics for motifs (some of the previously found could be already filtered out)
  	@arr3=split /\*/,$name;
  	
  	$motif_stat_final{$patterns_nm_hash{$arr3[1]}+1}++;
  	if ($arr3[1] eq "CYSRICH")
     {$motif_sigp_un_lt_cysrich++;}
    else
    {$motif_sigp_un_lt++;}
  	
	  $seq=$sequence_data_tmp{seq};
	  $seq=~s/\s+//g;
    $seq=~s/\*//g;
    $seq=~s/[^A-Za-z]//g;
    $seq2=uc($seq);
    printf OUT2 "%s\n%s\n",$name,$origseq{$name};
	}
	close $fh;
	close OUT2;
	unlink "${filename_noext}_tmp1.fasta";

#print statistics
print_log "Calculating statistics.";
open OUT,">$names{stats}" || die "Can not open $names{stats} for writing: $!\n";
unless ($start_with_sp)
{
  print OUT $print_header;
  printf OUT "Results for %s:\nTotal sequences read:%d\nNumber of motifs that were searched for:%d\n\nMotifs found:\n",$filename_source,$num_seqs,$num_motifs;
  printf OUT "motif_id\tnumber_of_seq_motifs\tnumber_of_seq_motifs_with_M\tfinal_number_of_motifs\tmotif_pattern\n";
  foreach $key (sort {$a<=>$b} keys(%motif_stat))
  {
  	printf OUT "%s\t%d\t%d\t%d\t%s\n",$patterns_nm[$key-1],$motif_stat{$key},$motif_stat_m{$key},$motif_stat_final{$key},format_pattern($patterns[$key-1]);
  }
  printf OUT "\nTotal sequences with motifs: %d\n",$motif_all-$motif_stat{$bignum};
  printf OUT "Total sequences with motifs and methionine before motif: %d\n",$motif_all_m-$motif_stat_m{$bignum};
}
else
{	 printf OUT "Results for %s\nTotal sequences read:%d\n\n",$filename_source,$num_sp; }

if (-e $names{signalp})
{ 
	printf OUT "Total sequences with motifs and methionine before motif that passed signalP: %d\n",$motif_sigp;
	printf OUT "Total sequences with motifs and methionine before motif that passed signalP, unique and having length <=%d: %d\n",$limit_length,$motif_sigp_un_lt;
	printf OUT "\nTotal cystein-rich sequences without correct motifs, but that have passed all other filters: %d\n",$motif_sigp_un_lt_cysrich;
}
	if ($n01 && $max[0]+$max[1]+$max[2]>0)
	{printf OUT "\n\nfull range pattern for 1st motif: CX{%d,%d}CXXXCX{%d,%d}CX{%d,%d}CXC\n\n",$min[0],$max[0],$min[1],$max[1],$min[2],$max[2];}

$dirname=$prefix."_results";
unless ($no_spada)
{
 #run spada and integrate results
 if ($spada_precalc_dir)
 #use pre-calculated spada results
 {
 	print_log "Getting pre-calculated spada results from '$spada_precalc_dir'.";
 	$cmd=sprintf "ln -s %s %s",$spada_precalc_dir,"${dirname}_spada";
 	system $cmd;
 }
 else
 {
   if ($spada_nucl)
   #use different (nucleotide) file for spada as input; useful when translation is skipped for main programs
   {$cmd=sprintf "%s -e 0.001 -d %s -s 0 -t 10 -f %s -c %s/cfg.txt >%s 2>&1",$spada,"${dirname}_spada",$spada_nucl,$spada_dir,$names{"spadalog"};}
   else
   {$cmd=sprintf "%s -e 0.001 -d %s -s 0 -t 10 -f %s -c %s/cfg.txt >%s 2>&1",$spada,"${dirname}_spada",$filename_source,$spada_dir,$names{"spadalog"};} 
   print_log "Starting spada as ".$cmd;
   system $cmd;
   
   if (get_num_str_by_grep($names{spadalog},"Pipeline successfully completed")==1)
    {print_log "Spada finished";}
   else
    #problems with spada
    {print_log "ERROR: Spada not finished or finished with errors. Spada results (if any) are unreliable." unless ($spada_precalc_dir);}
 }
 #get signalp splitting point
 %spsplit=();
 $cmd= sprintf "sort -k 2 %s > %s && sort %s > %s","${dirname}_spada/31_model_evaluation/61_final.gtb","${filename_noext}_1","${dirname}_spada/31_model_evaluation/35_sigp.tbl","${filename_noext}_2";
 system $cmd;
 open SP, "join -t \$'\t' -1 2 -2 1 ${filename_noext}_1 ${filename_noext}_2 \| cut -f2,22 |";
 while (<SP>)
 {
 	chomp;
 	@arr= split /\t/;
 	$spsplit{$arr[0]}=$arr[1];
 }
 close SP;
 unlink ("${filename_noext}_1","${filename_noext}_2");
  
 printf OUT "\nTotal sequences found by spada: %d\n", get_num_str_by_grep("${dirname}_spada/61_final.fasta",">");
 
 #combine results with spada and get unique sequences; all other factors held equal, spada results are removed when duplication found 
 system "cat $names{final} ${dirname}_spada/61_final.fasta > ${filename_noext}_tmp1.fasta";
  
 unless ($keep_subseqs)
 {
 	$cmd=sprintf "cysmotif_fasta_cleaner.pl -l %d %s > %s", $keep_subseqs_lg,"${filename_noext}_tmp1.fasta","${filename_noext}_tmp.fasta";
 	system $cmd;
 	unlink "${filename_noext}_tmp1.fasta";
 }
 
 open $fh,"${filename_noext}_tmp.fasta" || print LOG "WARNING: Can not open joined spada fasta file - skipping it: $!\n";
 %uniqseq=();
 
 while (read_fasta_sequence($fh, \%sequence_data_tmp)) 
  {
  	
  	$name=$sequence_data_tmp{header};
	  $seq=$sequence_data_tmp{seq};
	  $seq=~s/\s+//g;
    $seq=~s/\*//g;
    $seq=~s/[^A-Za-z]//g;
    $seq2=uc($seq);
    $uniqseq{$seq2}=$name;
    unless (exists $origseq{$name} || exists $tabseq{$name})
    #spada
     {
     	#add space at signal peptide end as in non-spada seqs
     	$n1=$name; $n1=~s/^>//;
     	$t1=substr($seq,0,$spsplit{$n1}); 
     	$t2=substr($seq,$spsplit{$n1});
     	$origseq{$name}=$t1." ".$t2;
     	$tabseq{$name}=sprintf "\t\t%s\t%s\t%s\t%s\t%s\n",$t1." ".$t2,"SPADA",rename_spada($n1,1),cformula($t2),""; 
     }
	}
	close $fh;
	open (OUT2,">$names{withspada}") || die "Can not open $names{withspada} for writing: $!\n";
	open (OUTT, ">$names{spadatabbed}") || die "Can not open $names{spadatabbed} for writing: $!\n";
   printf OUTT "\t\tseq\tmotif_name\tseq_name\tCformula\tmotif\n";
	$nu=0;
 
 foreach $key (sort {$uniqseq{$a} cmp $uniqseq{$b}} keys (%uniqseq))
  {
  	printf OUT2 "%s\n%s\n",rename_spada($uniqseq{$key},1),$origseq{$uniqseq{$key}}; 
  	$nu++;
  	print OUTT $tabseq{$uniqseq{$key}}
  }
  close OUT2;
  close OUTT;
 unlink "${filename_noext}_tmp.fasta";
 printf OUT "Total unique sequences including spada results: %d\n",$nu;
}
close OUT;

print_log "Moving files to output dir $dirname...";
unless ($start_with_sp)
{
  if (!$skip_translation && $names{translated}!~/\.gz$/ && $names{translated}!~/\.bz2$/)
  {
  	#archive file with newly translated seqs
  	$sz=-s $names{translated};
  	if ($sz > 10000000)
  	{system "pbzip2 -f $names{translated}";}
  }
 	#do not move translated file
  delete $names{translated};

	unless (-d $dirname)
	{
   $cmd=sprintf "mkdir %s",$dirname;
   system $cmd;
  }
  
   @flist=();
  foreach $file (keys (%names))
  {
  	if ($file ne "final" && $file ne "tabbed" && $file ne "stats" && $file ne "logf" && $file ne "withspada" && $file ne "spadatabbed" && $file ne "spada")
  	{
  		push @flist,$names{$file};
  		delete $names{$file};
  	}
  }
  unless ($no_spada)
  	{
  		push @flist,$names{final}; delete $names{final};
  	  push @flist,$names{tabbed}; delete $names{tabbed};
  	  $cmd=sprintf "cp %s %s", "${dirname}_spada/61_final.fasta",$names{spada};
      system $cmd;
      push @flist,$names{spada}; delete $names{spada};
      if ($spada_precalc_dir) 
      #delete copy of spada results - do not need them anymore
      {system "rm -rf ${dirname}_spada";}
      #delete spada files with translated seqs
      else 
      { system "rm -f ${dirname}_spada/01_preprocessing/11_refseq_trans6.fas ${dirname}_spada/01_preprocessing/12_orf_genome.fas";}
    }

  #zip auxiliary files
  system "zip -9 -q $filename_noext @flist";
  $names{zip}=$filename_noext.".zip";
  system "rm -f @flist";
  
  #move files to output dir
  foreach $file (keys (%names))
  {
    move $names{$file},$dirname || print STDERR "Error with moving $names{$file} : $!\n";
  }
  $cmd=sprintf "cp %s %s",$motif_file,$dirname;
  system $cmd;

}
 print LOG "\n";
 print_log(program_time($start,time,"Pipeline",1,$delay)); 
 print LOG "Results can be found in '${dirname}' dir.\n";
 unless ($no_spada) 
 { 
 	if ($spada_precalc_dir)
 	{$outdir=$spada_precalc_dir;}
 	else
 	{$outdir="${dirname}_spada";}
 	
 	 print LOG "Additional spada data can be found in '$outdir' dir.\n";
 }
 close LOG;

