#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# FastaSplit.pl                                             |
#                                                           |
#-----------------------------------------------------------+
#   AUTHOR: James C. Estill                                 |
#  CONTACT: JamesEstill at gmail.com                        |
#  STARTED: 04/19/2006                                      |
# MODIFIED: 01/18/2007                                      |
#                                                           |
# DESCRIPTION:                                              |
#  Takes an input file in fasta format and splits it into   |
#  a number of fasta files. The number of fasta files that  |
#  are created is a varaible pased by the user at the       |
#  command line. Useful for creating the FASTA files that   |
#  will be used in BLAST queries sent to a cluster.         |
#                                                           |
# EXAMPLE USAGE:                                            |
#                                                           |
#  FastaSplit.pl InputFastaFile NumFile                     |
#  FastaSplit.pl ArabidopsisGenes.fasta 6 Arabidopsis       |
#                                                           |
#-----------------------------------------------------------+

$pad = 0;

# TO DO
# Modify this to use PERL switches

use Bio::SeqIO;      # Use Bio::SeqIO for the reading and 
                     # writing of seq files in different formats
use Cwd;             # Use the Cwd module to get the current working directory

# Get command-line arguments, or die with a usage statement
my $usage = "Fasta2Dir.pl infile NumFiles OutRoot \n";
my $infile = shift or die $usage;     
                     # The full path of the infile that will be transformed
my $NumFiles = shift or die $usage;                    
                     # The format that the input file is in. 
                     # Any valid bioperl format (ie. abi, ace, gcg, genbank, 
                     # fasta, swiss, tigr etc.) 
my $outfile = shift or die $usage;                         
                     # The name that will be used as the prefix for 
                     # the output files. (ie. Output0000001, Output0000002, 
                     #..etc)
my $CurrentDir = cwd();
my $OutputDir = $CurrentDir."/".$outfile."/";              
                     # Set the output dirctory that the files 
                     # will be created in
mkdir $OutputDir;    # Make the output directory, if the output directory 
                     # exists the existing dir will not be overwritten
my @SeqsOut;         # Array to hold references to seq out objects

$SeqNum = 1;         # Counter to keep track of the number of seqs
$OutCount = 1;       # The OutCount will be used to keep track of
                     # which output file the fasta file will go to.

# create one SeqIO object to read in,and another to write out
$infileformat = "fasta";
my $seq_in = Bio::SeqIO->new('-file' => "<$infile",
                             '-format' => $infileformat);

#-----------------------+
# OPEN A SEPARATE SEQ   |
# OUT OBJECT FOR EACH   |
# FILE.                 |
#-----------------------+
for ( $i = 1; $i<=$NumFiles; $i++ )
{

    my $num;
    print "Creating File $i \n";
    if ($i == 100) {exit;}

    if ($pad) {
	$num = sprintf("%3d", $i);  
	# Pad number with zeros so that the total length
	# of the string is 3 characaters (ie. 012)
	$num=~ tr/ /0/;
    }
    else {
	$num = $i;
    }

    $OutFilePath = $OutputDir.$outfile.$num.".fasta";
    $SeqsOut[$i] =  Bio::SeqIO->new('-file' => ">$OutFilePath",
                                    '-format' => "fasta");
    
}
#exit;

while (my $inseq = $seq_in->next_seq)
{

      print "SeqIn:".$SeqNum."\tSeqOut:".$OutCount." \n";    
                     # Print the number of the intput and output file nums

      # Write the seq data out to the appropriate fasta file
      $SeqsOut[$OutCount]->write_seq($inseq);
      #$seq_out->write_seq($inseq);  # Write the individual fasta file

      $SeqNum++;                    # Increment SeqNumber
      $OutCount++;
      if ( $OutCount > $NumFiles) {$OutCount = 1;}

}

print "The files have all been placed in the directory: \n";
print $OutputDir."\n";

exit;

#
# 07/14/2007
# - Existing program moved to the jperl site
