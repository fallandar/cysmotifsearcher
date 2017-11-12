package cysmotif;
#Last update: May 04, 2017
BEGIN{}
use Exporter;
our @ISA    = qw(Exporter);
our @EXPORT = qw(codon2aa get_num_str_by_grep program_time read_fasta_sequence revcomp);

our @ACGT=("A","C","G","T"); 

sub codon2aa
{
my($codon)=@_;
$codon=uc $codon;
my(%g)=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y',
'TAA'=>'_','TAG'=>'_','TGC'=>'C','TGT'=>'C','TGA'=>'_','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P',
'CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M',
'ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V',
'GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G',
'GGG'=>'G','GGT'=>'G');
 if(exists $g{$codon})
 {
  return $g{$codon};
 }
 elsif ($codon=~/N/)
 {
 	return '_';
 }
 else
 {
  print STDERR "Bad codon \"$codon\"!!\n";
 }
}#codon2aa

sub get_num_str_by_grep 
{
	#get and return number of occurences of given string (2nd param) in given file (plain, gz, bz2) (1st param)
 my ($filename,$str) = @_;
 if ($filename=~/\.gz$/)
  {open (GREP,"zcat $filename \| grep -c \"$str\" |") || print LOG "WARNING: Can not parse $filename: $!\n";}
 elsif	($filename=~/\.bz2$/)
  {open (GREP,"pbzip2 -c -d $filename \| grep -c \"$str\" |") || print LOG "WARNING: Can not parse $filename: $!\n";}
 else
  {open (GREP,"grep -c \"$str\" $filename |") || print LOG "WARNING: Can not parse $filename: $!\n";} 	
 my $ln=<GREP>; 
 chomp $ln;
 close GREP;
 return $ln;
} #get_num_seqs_fasta

sub read_fasta_sequence
{
   my ($fh, $seq_info) = @_;

   $seq_info->{seq} = undef; # clear out previous sequence

   # put the header into place
   $seq_info->{header} = $seq_info->{next_header} if $seq_info->{next_header};

   my $file_not_empty = 0; 
   while (<$fh>) {
      $file_not_empty = 1;
      next if /^\s*$/;  # skip blank lines
      chomp;    

      if (/^>/) { # fasta header line
         my $h = $_;    
         #$h =~ s/^>//;  
         if ($seq_info->{header}) {
            $seq_info->{next_header} = $h;
            return $seq_info;   
         }              
         else { # first time through only
            $seq_info->{header} = $h;
         }              
      }         
      else {    
         s/\s+//;  # remove any white space
         $seq_info->{seq} .= $_;
      }         
   }    

   if ($file_not_empty) {
      return $seq_info;
   }    
   else {
      # clean everything up
      $seq_info->{header} = $seq_info->{seq} = $seq_info->{next_header} = undef;

      return;   
   }    
}#read_fasta_sequence

sub revcomp
{
	#get reversed complementary sequence for input seq
	my $tmp=reverse $_[0];
  $tmp=~tr/ATCGatcg/TAGCtagc/;
  return $tmp;
}#revcomp 

sub program_time
{
	#calculate time interval of program (or any part of program) execution
	my ($start,$finish,$pr,$mode,$del)=@_;
	
	my ($duration,$min,$sec,$hour); my $info="";
  $duration = $finish - $start - $del*$mode;
  $min=int($duration/60); $sec=$duration % 60;
  if ($min >= 60)
  {$hour=int($min/60); $min=$min % 60; $info=sprintf "%s finished succesfully in\t%d\thours\t%d\tminutes\t%d\tseconds.\n",$pr,$hour,$min,$sec;}
  else
  {$info=sprintf "%s finished succesfully in\t%d\tminutes\t%d\tseconds.\n",$pr,$min,$sec;}
  
  if ($del && $mode)
  {
   $min=int($del/60); $sec=$del % 60;
   if ($min >= 60)
    {$hour=int($min/60); $min=$min % 60; $info.=sprintf "(Starting delay was set to\t%d hours\t%d minutes\t%d\tseconds.)\n",$hour,$min,$sec;}
   else
    {$info.=sprintf "(Starting delay was set to\t%d\tminutes\t%d\tseconds.)\n",$min,$sec;}
  }
  return $info;
}#program_time
 
1;
END{}