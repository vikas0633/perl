#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_genemark2gff.pl - Convert genemark.hmm output to gff  |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_at_sourceforge.net                   |
# STARTED: 10/30/2007                                       |
# UPDATED: 01/12/2010                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Converts output from the genemark.hmm gene prediction    |
#  program to the GFF format.                               |
#                                                           |
# USAGE:                                                    |
#  cnv_genemark2gff.pl -i infile.out -o outfile.gff         |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

package DAWGPAWS;

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;
use File::Copy;
use Getopt::Long;
use Bio::Tools::Genemark;      # Genemark program from Bio::Tools
# The following needed for printing help
use Pod::Select;               # Print subsections of POD documentation
use Pod::Text;                 # Print POD doc as formatted text file
use IO::Scalar;                # For print_help subfunction
use IO::Pipe;                  # Pipe for STDIN, STDOUT for POD docs
use File::Spec;                # Convert a relative path to an abosolute path

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev: 948 $ =~ /(\d+)/;
# Get GFF version from environment, GFF2 is DEFAULT
my $gff_ver = uc($ENV{DP_GFF}) || "GFF2";

#-----------------------------------------------------------+
# VARIABLE SCOPE                                            |
#-----------------------------------------------------------+

# GENERAL PROGRAM VARS
my @inline;                    # Parse of an input line
my $msg;                       # Message printed to the log file

# DIR PATHS
my $infile;                    # Input path of *.out file to convert
my $outfile;                   # Output gff file path
my $prefix;                    # Prefix name for gff output file

# FILE PATHS
my $rm_path;                   # Full path to the repeatmasker binary
my $ap_path;                   # Full path to the apollo program

my $src_prog = "GeneMarkHMM";        # Source program and matrix
my $src_seq = "unknown_src";         # Set source seq to unknown by default
my $param;                     # Parameter file used (ie wheat)

# BOOLEANS
my $show_help = 0;             # Show program help
my $show_version = 0;          # Show program version
my $show_man = 0;              # Show program manual page using peldoc
my $show_usage = 0;            # Show program usage command             
my $quiet = 0;                 # Boolean for reduced output to STOUT
my $apollo = 0;                # Path to apollo and apollo variables
my $test = 0;
my $verbose = 0;
my $debug = 0;                 # Run the program in debug mode 

# PROGRAM COMMAND STRINGS
#my $cmd_repmask;               # Command to run RepeatMasker
#my $cmd_make_gff_db;           # Make the gff file for an individual database
#my $cmd_make_gff_el;           # Appears to not be used

# COUNTERS AND INDEX VARS
my $i;                         # Used to index the config file

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(
		    # REQUIRED ARGUMENTS
		    "i|infile=s"   => \$infile,
                    "o|outfile=s"  => \$outfile,
		    # ADDITIONAL OPTIONS
		    "gff-ver=s"    => \$gff_ver,
		    "program=s"    => \$src_prog, 
		    "seqname=s"    => \$src_seq,
                    "parameter=s"  => \$param,
		    # BOOLEANS
		    "verbose"      => \$verbose,
		    "debug"        => \$debug,
		    "test"         => \$test,
		    "usage"        => \$show_usage,
		    "version"      => \$show_version,
		    "man"          => \$show_man,
		    "h|help"       => \$show_help,
		    "q|quiet"      => \$quiet,);

#-----------------------------+
# STANDARDIZE GFF VERSION     |
#-----------------------------+
unless ($gff_ver =~ "GFF3" || 
	$gff_ver =~ "GFF2") {
    # Attempt to standardize GFF format names
    if ($gff_ver =~ "3") {
	$gff_ver = "GFF3";
    }
    elsif ($gff_ver =~ "2") {
	$gff_ver = "GFF2";
    }
    else {
	print "\a";
	die "The gff-version \'$gff_ver\' is not recognized\n".
	    "The options GFF2 or GFF3 are supported\n";
    }
}

#-----------------------------+
# SHOW REQUESTED HELP         |
#-----------------------------+
if ( ($show_usage) ) {
#    print_help ("usage", File::Spec->rel2abs($0) );
    print_help ("usage", $0 );
}

if ( ($show_help) || (!$ok) ) {
#    print_help ("help",  File::Spec->rel2abs($0) );
    print_help ("help",  $0 );
}

if ($show_man) {
    # User perldoc to generate the man documentation.
    system ("perldoc $0");
    exit($ok ? 0 : 2);
}

if ($show_version) {
    print "\ncnv_genemark2gff.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}

genemark_to_gff ($src_seq, $src_prog, $infile, $outfile, $param);

exit;

#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub genemark_to_gff {
    
    my ($gm_src_seq, $gm_src_prog, $gm_in_path, $gff_out_path, 
	$parameter) = @_;

    my $attribute;
    my $source = $gm_src_prog;

    #-----------------------------+
    # OPEN THE GFF OUTFILE        |
    #-----------------------------+
    # Default to STDOUT if no argument given
    if ($gff_out_path) {
	open (GFFOUT, ">$gff_out_path") ||
	    die "ERROR: Can not output gff output file:\n$gff_out_path\n"
    }
    else {
	open (GFFOUT, ">&STDOUT") ||
	    die "Can not print to STDOUT\n";
    }
    if ($gff_ver =~ "GFF3") {
	print GFFOUT "##gff-version 3\n";
    }

    #-----------------------------+
    # OPEN THE GENEMARK INFILE    |
    #-----------------------------+
    my $gm_obj;
    if ($gm_in_path) {
	$gm_obj = Bio::Tools::Genemark->new(-file => $gm_in_path);
    }
    else {
	$gm_obj = Bio::Tools::Genemark->new(-fh => \*STDIN);
    }
    
    if ($parameter) {
	$gm_src_prog = $gm_src_prog.":".$parameter;
    }

    my $rna_count = 0;
    my $gene_num = 0;

    while(my $gene = $gm_obj->next_prediction()) {
       
	$gene_num++;     
	my $gene_id = sprintf("%05d", $gene_num); # Pad the gene number
	my $gene_name = $gm_src_prog."_gene_".$gene_id."\n";
	$gene_id = "gene".$gene_id;

	# The following is a duplicate of gene_num
	$rna_count++;
	my $rna_id = sprintf("%05d", $rna_count);
	# End duplicate

	my @exon_ary = $gene->exons();

	my $num_exon = @exon_ary;

	#-----------------------------+
	# Gene span for GFF3          |
	#-----------------------------+
	if ($gff_ver =~ "GFF3") {

	    #-----------------------------+
	    # GENE                        |
            #-----------------------------+
	    $attribute = "ID=".$source."_".$gene_id;

	    # Get gene start and end
	    my $gene_start = $gene->start();
	    my $gene_end = $gene->end();
	    if ($gene_start > $gene_end) {
		$gene_end =  $gene->start();
		$gene_start = $gene->end();
	    }

	    # Get gene score
	    my $gene_score;
	    if ($gene->score()) {
		$gene_score = $gene->score();
	    }
	    else {
		$gene_score = ".";
	    }

	    # Get gene strand
	    my $gene_strand = $gene->strand()."\t";
	    if ($gene_strand =~ "-1") {
		$gene_strand = "-";
	    }
	    else {
		$gene_strand = "+";
	    }

	    print GFFOUT  $gm_src_seq."\t".     # Seqname
		$source."\t".                   # Source
		"gene\t".                       #feature
		"$gene_start\t".                # start
		"$gene_end\t".                  # end
		"$gene_score\t".                # score
		"$gene_strand\t".               # strand
		".\t".                          # frame
		$attribute.                     # attribute
		"\n";

	}

	#-----------------------------+
	# EXONS                       |
	#-----------------------------+
	my $exon_num = 0;
	for my $ind_exon (@exon_ary) {
	    
	    $exon_num++;
	    my $exon_id = sprintf("%05d", $exon_num); # Pad the exon number
	    $exon_id = "exon".$exon_id;

	    # Set Feature
	    my $feature = "exon";

	    #-----------------------------+
	    # GET START AND END           |
	    #-----------------------------+
	    my $start = $ind_exon->start;
	    my $end = $ind_exon->end;
	    my $strand = $ind_exon->strand;
	    if ($start > $end) {
		$end =  $ind_exon->start();
		$start = $ind_exon->end();
	    }

	    #-----------------------------+
	    # FORMAT STRAND               |
	    #-----------------------------+
	    if ($strand =~ '-1') {
		$strand = "-"; 
	    }
	    elsif ($strand =~ '1') {
		$strand = "+";
	    }
	    else {
		$strand = ".";
	    }

	    #-----------------------------+
	    # GET SCORE                   |
	    #-----------------------------+
	    # It does not look like a score is assigned
	    # but will put this here in case one shows up
	    # in the future
	    my $score = $ind_exon->score() || ".";

	    #-----------------------------+
	    # SET ATTRIBUTE               |
	    #-----------------------------+
	    if ($gff_ver =~ "GFF3") {
		$attribute = "ID=".$source."_".$gene_id."_".$exon_id.
		    "\;Parent=".$source."_".$gene_id;
	    }
	    else {
		$attribute = $gene_id;
	    }

	    #-----------------------------+
	    # GET FRAME
	    #-----------------------------+
	    #my $frame = $ind_exon->frame();
	    my $frame = ".";

	    #-----------------------------+
	    # PRINT GFF OUTPUT            |
	    #-----------------------------+
	    # The following produces output for Apollo
	    print GFFOUT $gm_src_seq."\t".   # seq name
		$gm_src_prog."\t".           # source
		$feature."\t".               # feature
		$start."\t".                 # start
		$end."\t".                   # end
		$score."\t".                 # score
		$strand."\t".                # strand
		$frame."\t".                 # frame
		$attribute.                  # attribute		
		"\n";      

	} # End of for each exon

    } # End of for each gene model

    $gm_obj->close();

    if ($gff_out_path) {
	close (GFFOUT);
    }

}


sub print_help {
    my ($help_msg, $podfile) =  @_;
    # help_msg is the type of help msg to use (ie. help vs. usage)
    
    print "\n";
    
    #-----------------------------+
    # PIPE WITHIN PERL            |
    #-----------------------------+
    # This code made possible by:
    # http://www.perlmonks.org/index.pl?node_id=76409
    # Tie info developed on:
    # http://www.perlmonks.org/index.pl?node=perltie 
    #
    #my $podfile = $0;
    my $scalar = '';
    tie *STDOUT, 'IO::Scalar', \$scalar;
    
    if ($help_msg =~ "usage") {
	podselect({-sections => ["SYNOPSIS|MORE"]}, $0);
    }
    else {
	podselect({-sections => ["SYNOPSIS|ARGUMENTS|OPTIONS|MORE"]}, $0);
    }

    untie *STDOUT;
    # now $scalar contains the pod from $podfile you can see this below
    #print $scalar;

    my $pipe = IO::Pipe->new()
	or die "failed to create pipe: $!";
    
    my ($pid,$fd);

    if ( $pid = fork() ) { #parent
	open(TMPSTDIN, "<&STDIN")
	    or die "failed to dup stdin to tmp: $!";
	$pipe->reader();
	$fd = $pipe->fileno;
	open(STDIN, "<&=$fd")
	    or die "failed to dup \$fd to STDIN: $!";
	my $pod_txt = Pod::Text->new (sentence => 0, width => 78);
	$pod_txt->parse_from_filehandle;
	# END AT WORK HERE
	open(STDIN, "<&TMPSTDIN")
	    or die "failed to restore dup'ed stdin: $!";
    }
    else { #child
	$pipe->writer();
	$pipe->print($scalar);
	$pipe->close();	
	exit 0;
    }
    
    $pipe->close();
    close TMPSTDIN;

    print "\n";

    exit 0;
   
}

__END__

=head1 NAME

cnv_genemark2gff.pl - Convert genemark output to gff format

=head1 VERSION

This documentation refers to program version $Rev: 948 $

=head1 SYNOPSIS

=head2 Usage

    cnv_genemark2gff.pl -i infile.genemark -o outfile.gff

=head2 Required Arguments

    --infile     # Path to the input file to translate
                 # If not provided, assumes input from STDIN
    --outfile    # Path to the output gff file
                 # If not provided, writes output to STDOUT
    --seqname    # The id of the sequence analyzed

=head1 DESCRIPTION

Converts the output from the genemark.hmm program to the gff format. This
has been tested to work with gmhmme2 and gmhmme3. All exons are currently
tagged as 'exon' in the gff output file. This is for compatibility with
the Apollo genome annotation curation program.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path to the genemark file to translate to gff. If an infile is not
specified, then the program will expect input from standard input.

=item -o,--outfile

Path to the gff output file. If an outfile is not specified, the progrm
will write the gff file to standard output.

=back

=head1 OPTIONS

=over 2

=item --gff-ver

The GFF version for the output. This will accept either gff2 or gff3 as the
options. By default the GFF version will be GFF2 unless specified otherwise.
The default GFF version for output can also be set in the user environment
with the DP_GFF option. The command line option will always override the option
defined in the user environment. 

=item --seqname

This is the value listed as the source sequence in the gff output file. While
not a specifically required variable, the default value for this in unknown.
This will generally be set to the BAC ID or contig ID.

=item --program

This is the source program name used in the gff output file. By default this
is set to be GeneMarkHMM. This option allows you to set the source program
to any value that you would want.

=item --usage

Short overview of how to use program from command line.

=item --help

Show program usage with summary of options.

=item --version

Show program version.

=item --man

Show the full program manual. This uses the perldoc command to print the 
POD documentation for the program.

=item --verbose

Run the program with maximum output.

=item --test

Run the program without doing the system commands.

=back

=head1 EXAMPLES

=head2 Typical Use

The typical use of this program will be to parse a file produce from the
genemark.hmm program.

    cnv_genemark2gff.pl -i HEX2493A05_genemark_hv.out --seqname HEX2493A05
                        -o HEX2493A05_genemark_hv.gff

This will produce a gff output file similar to the following:

    HEX2493A05 GeneMarkHMM	exon	683	1393	.  +  .	RNA0001
    HEX2493A05 GeneMarkHMM	exon	1736	2084	.  +  .	RNA0001
    HEX2493A05 GeneMarkHMM	exon	2195	2515	.  +  .	RNA0001
    HEX2493A05 GeneMarkHMM	exon	2696	2803	.  +  .	RNA0001
    HEX2493A05 GeneMarkHMM	exon	2918	3035	.  +  .	RNA0001
    HEX2493A05 GeneMarkHMM	exon	3058	3131	.  +  .	RNA0001
    HEX2493A05 GeneMarkHMM	exon	3219	3502	.  +  .	RNA0001
    HEX2493A05 GeneMarkHMM	exon	3552	3559	.  +  .	RNA0002
    HEX2493A05 GeneMarkHMM	exon	3711	3801	.  +  .	RNA0002
    HEX2493A05 GeneMarkHMM	exon	3947	4711	.  +  .	RNA0002
    ...

The --seqname option used above allows you to specify the value written in the 
first column of the gff file. If the --seqname was not specified
like the following:

    cnv_genemark2gff.pl -i HEX2493A05_genemark_hv.out
                        -o HEX2493A05_genemark_hv.gff

The gff output would be similar to the following:

    unknown_src GeneMarkHMM	exon	683	1393	.  +  .	RNA0001
    unknown_src GeneMarkHMM	exon	1736	2084	.  +  .	RNA0001
    unknown_src GeneMarkHMM	exon	2195	2515	.  +  .	RNA0001
    unknown_src GeneMarkHMM	exon	2696	2803	.  +  .	RNA0001
    unknown_src GeneMarkHMM	exon	2918	3035	.  +  .	RNA0001
    unknown_src GeneMarkHMM	exon	3058	3131	.  +  .	RNA0001
    unknown_src GeneMarkHMM	exon	3219	3502	.  +  .	RNA0001
    unknown_src GeneMarkHMM	exon	3552	3559	.  +  .	RNA0002
    unknown_src GeneMarkHMM	exon	3711	3801	.  +  .	RNA0002
    unknown_src GeneMarkHMM	exon	3947	4711	.  +  .	RNA0002
    ...

=head2 Specify the Training Matrix Used

It is also possible to designate the second column of the gff output file
using the --program option.
This can be used to specify the training data use for gene predictions.
This will allow you to later separate gene models for different
training data sets.
For example if I used the wheat training matrix, I may do the following:

    cnv_genemark2gff.pl -i HEX2493A05_genemark_hv.out --seqname HEX2493A05
                        -o HEX2493A05_genemark_hv.gff 
                        --program GeneMark:wheat

This will produce output similar to the following:

    HEX2493A05 GeneMark:wheat	exon	683	1393	.  +  .	RNA0001
    HEX2493A05 GeneMark:wheat	exon	1736	2084	.  +  .	RNA0001
    HEX2493A05 GeneMark:wheat	exon	2195	2515	.  +  .	RNA0001
    HEX2493A05 GeneMark:wheat	exon	2696	2803	.  +  .	RNA0001
    HEX2493A05 GeneMark:wheat	exon	2918	3035	.  +  .	RNA0001
    HEX2493A05 GeneMark:wheat	exon	3058	3131	.  +  .	RNA0001
    HEX2493A05 GeneMark:wheat	exon	3219	3502	.  +  .	RNA0001
    HEX2493A05 GeneMark:wheat	exon	3552	3559	.  +  .	RNA0002
    HEX2493A05 GeneMark:wheat	exon	3711	3801	.  +  .	RNA0002
    HEX2493A05 GeneMark:wheat	exon	3947	4711	.  +  .	RNA0002
    ...

=head1 DIAGNOSTICS

The error messages that can be generated will be listed here.

=head1 CONFIGURATION AND ENVIRONMENT

=over 2

=item DP_GFF

The DP_GFF variable can be defined in the user environment to set
the default GFF version output. Valid settings are 'gff2' or
'gff3'.

=back

=head1 DEPENDENCIES

=head2 Required Software

=over

=item * genemark

This program parses output from the genemark.hmm program. This program is
available from http://opal.biology.gatech.edu/GeneMark/

=back

=head2 Required Perl Modules

=over

=item * Bio::Tools::Genemark

This module is part of bioperl. Information on installing biperl is
available from:
http://bioperl.open-bio.org/wiki/Installing_BioPerl

=back

=head1 BUGS AND LIMITATIONS

=over 2

=item * No bugs currently known 

If you find a bug with this software, file a bug report on the DAWG-PAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=204962

=back

=head2 Limitations

=over

=item * Limited testing with different versions of genemark.hmm

This program is known to work with output produced from gmhmme2 and gmhmme3.

=back

=head1 SEE ALSO

Thhis program is part of the DAWGPAWS package of genome
annotation programs. See the DAWGPAWS web page 
( http://dawgpaws.sourceforge.net/ )
or the Sourceforge project page 
( http://sourceforge.net/projects/dawgpaws ) 
for additional information about this package.

=head1 REFERENCE

Your use of the DAWGPAWS programs should reference the following manuscript:

JC Estill and JL Bennetzen. 2009. 
"The DAWGPAWS Pipeline for the Annotation of Genes and Transposable 
Elements in Plant Genomes." Plant Methods. 5:8.

The use of the cnv_genemark2gff.pl  program should also reference the 
appropriate GeneMark manuscript listed at:

http://exon.biology.gatech.edu/references.html

=head1 LICENSE

GNU General Public License, Version 3

L<http://www.gnu.org/licenses/gpl.html>

THIS SOFTWARE COMES AS IS, WITHOUT ANY EXPRESS OR IMPLIED
WARRANTY. USE AT YOUR OWN RISK.

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 10/30/2007

UPDATED: 04/15/2010

VERSION: $Rev: 948 $

=cut

#-------------------------------------------------------+
# HISTORY                                               |
#-------------------------------------------------------+
# 10/30/2007
# - Main body of program established
# 01/20/2009
# - Large cleanup of code, previous version was a copy
#   and paste job from a repeatmasker program and needed
#   massive cleanup. 
# - Changed to print_help subfunction that extracts
#   help from POD documentation.
# - Dropped the apollo convert subfunction
# 01/21/2009
# - Added support for input from STDIN if -i option
#   not given.
# 03/27/2009
# - Renamed --src-seq to --seqname
# - Fixed bug where every gene model was being placed in the positive strand
# 01/12/2010
# - Added option for gff3 output
# - Using a program name for the source that can have a suffix for the
# - Added default printing to STDOUT if gffout not defined
# - Modified gene id for the gff2 version
# - Updated POD to take these changes into account
# 04/15/2010
# -Updated reference and license
