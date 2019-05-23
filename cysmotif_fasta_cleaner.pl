#!/usr/bin/perl

use cysmotif qw (read_fasta_sequence rename_spada);
$program="cysmotif_fasta_cleaner.pl";  
$version="1.0";
$last_update="May 23, 2019";
$comment="Written by Andrew Shelenkov, VIGG of RAS";
$usage="$program [OPTIONS] cysmotif_fasta_file(s)...";

$batch_mode=0;
$difflg=5; #default length difference
$spada_rename=0; #rename spada seqnames

if ($ARGV[0] eq "-h" || $ARGV[0] eq "-help" || $ARGV[0] eq "--help")
{
  print "You are using $program version $version \n";
  printf "Usage: $usage\n";
}

if (($ARGV[0] eq "--version") || ($ARGV[0] eq "-v") )
{ #print the number of version and date of release
  print "You are using $program version $version\n";
  print "Last update: $last_update\n";
  print "$comment\n";
  exit;
} 

while ($ARGV[0]=~ /^-/)
 #reading options and keys
{
 if ($ARGV[0]=~ /-b/)	
 #set batch mode of input file processing; output will go to *uniq_filt.fasta
 {$batch_mode=1;}
 if ($ARGV[0]=~ /-l/)
 #change default length to check sequences
 {$difflg=$ARGV[1]; shift @ARGV; die "Wrong lg set with -l option. Should be >=1" if ($difflg<1);}
 if ($ARGV[0]=~ /-s/)	
 #set batch mode of input file processing; output will go to *uniq_filt.fasta
 {$spada_rename=1;}
	
	shift @ARGV;
} 

while ($ARGV[0])
{
	$filename_source=shift;
	
	if ($batch_mode)
	{
	 ($filename_out=$filename_source)=~s/.gz$//;
	  $filename_out=~s/.bz2$//;
	  $filename_out=~s/uniq.fasta$/uniq_filt.fasta/;
	}
	
if ($filename_source=~/\.gz$/)
 {open ($fh,"zcat $filename_source |") || die "Can not open $filename_source: $!\n";}
 elsif ($filename_source=~/\.bz2$/)
 {open ($fh,"pbzip2 -c -d $filename_source |") || die "Can not open $filename_source: $!\n";}
 else
 { open ($fh,$filename_source) || die "Can not open $filename_source: $!\n";}
 
 %uniqseq=();
while (read_fasta_sequence($fh, \%sequence_data_tmp)) 
{
	#get uniq sequences
	$id=$sequence_data_tmp{header};
	$seq=$sequence_data_tmp{seq}; 
	
	  $seq=~s/\s+//g;
    $seq=~s/\*//g;
    $seq=~s/[^A-Za-z]//g;
    $seq2=uc($seq);
    $uniqseq{$seq2}=$id;
}
close $fh;

#remove subseqs
 	foreach $key (keys (%uniqseq))
 	{
 		foreach $key2 (keys (%uniqseq))
 		{
 			next if ($key2 eq $key);
 			next unless ($uniqseq{$key2} && $uniqseq{$key});
 			if ((index ($key,$key2)==0) && (length ($key) - length ($key2) <=$difflg ) )
 			{
 				$flag=1;
 				if ($spada_rename)
 				{
 					($tmp1=$uniqseq{$key2})=~s/_0M_1$//;
 					
 					if ($tmp1=~/^>crp[0-9]+_(.*)$/)
 					{$tmp2=$1;}
 					#print STDERR "EEE:tmp1=$tmp1 tmp2=$tmp2 $uniqseq{$key}\n";
 					$flag=0 if (index ($uniqseq{$key},$tmp2)==-1);
 				}
 				if ($flag)
 				{
 			  	#print STDERR "Found subseq of $uniqseq{$key2} in $uniqseq{$key}. Deleting $uniqseq{$key2}\n$uniqseq{$key}\n$key\n$uniqseq{$key2}\n$key2\n\n-------------------\n\n";
 				 $uniqseq{$key2}=undef;
 				}
 			}
 		}
 	}
 	open OUT,">$filename_out" || die "Can not open output $filename_out: $!\n" if ($batch_mode);;
 	foreach $key (keys (%uniqseq))
 	{
 		if ($uniqseq{$key})
 		{
 			$topr=$uniqseq{$key};
 			if ($spada_rename)
 			{$topr=rename_spada($uniqseq{$key});}
 			
 		 if ($batch_mode)
 		 {print OUT "$topr\n$key\n";}
 		 else
 		 {print "$topr\n$key\n";}
 		}
 	}
 	close OUT if ($batch_mode);;
}
 