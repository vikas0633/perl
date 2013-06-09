#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_geneid.pl - Run geneid gene prediction in batch mode|
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 11/06/2008                                       |
# UPDATED: 03/11/2010                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Short Program Description                                |
#                                                           |
# USAGE:                                                    |
#  batch_geneid.pl -i in_dir/ -o base_out_dir/              |
#                                                           |
# VERSION: $Rev: 948 $                                      |
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
use Getopt::Long;
# The following needed for printing help
use Pod::Select;               # Print subsections of POD documentation
use Pod::Text;                 # Print POD doc as formatted text file
use IO::Scalar;                # For print_help subfunction
use IO::Pipe;                  # Pipe for STDIN, STDOUT for POD docs
use File::Spec;                # Convert a relative path to an abosolute path

#use Bio::Tools::Geneid;

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev: 948 $ =~ /(\d+)/;
# Get GFF version from environment, GFF2 is DEFAULT
my $gff_ver = uc($ENV{DP_GFF}) || "GFF2";

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $indir;
my $outdir;
my $geneid_param_file;         # File containing the parameters for geneid

# The path to the binary for running the geneid program
# by defalt this will assume it is in the user's Path variable
# under the default name
my $geneid_path = $ENV{GENEID_BIN} || "geneid";  

# BOOLEANS
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_test = 0;                  # Run the program in test mode


#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED
		    "i|indir=s"     => \$indir,
                    "o|outdir=s"    => \$outdir,
		    "p|param=s"     => \$geneid_param_file,
		    # OPTIONS
		    "gff-ver=s"    => \$gff_ver,
		    "geneid-path=s" => \$geneid_path,
		    "q|quiet"       => \$quiet,
		    "verbose"       => \$verbose,
		    # ADDITIONAL INFORMATION
		    "usage"         => \$show_usage,
		    "test"          => \$do_test,
		    "version"       => \$show_version,
		    "man"           => \$show_man,
		    "h|help"        => \$show_help,);

#-----------------------------+
# SHOW REQUESTED HELP         |
#-----------------------------+
if ( $show_usage ) {
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
    print "\nbatch_geneid.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}

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
# CHECK REQUIRED ARGS         |
#-----------------------------+
if ( (!$indir) || (!$outdir) || (!$geneid_param_file)  ) {
    print "\a";
    print STDERR "\n";
    print STDERR "ERROR: An input directory was not specified at the".
	" command line\n" if (!$indir);
    print STDERR "ERROR: An output directory was not specified at the".
	" command line\n" if (!$outdir);
    print STDERR "ERROR: A geneid paremter file was not specified at".
	" command line\n" if (!$geneid_param_file);
    print_help ("usage", $0 );
}

#-----------------------------+
# CHECK FOR SLASH IN DIR      |
# VARIABLES                   |
#-----------------------------+
# If the indir does not end in a slash then append one
# TO DO: Allow for backslash
unless ($indir =~ /\/$/ ) {
    $indir = $indir."/";
}

unless ($outdir =~ /\/$/ ) {
    $outdir = $outdir."/";
}

#-----------------------------+
# Get the FASTA files from the|
# directory provided by the   |
# var $indir                  |
#-----------------------------+
opendir( DIR, $indir ) || 
    die "Can't open directory:\n$indir"; 
my @fasta_files = grep /\.fasta$|\.fa$/, readdir DIR ;
closedir( DIR );

my $count_files = @fasta_files;

#-----------------------------+
# SHOW ERROR IF NO FILES      |
# WERE FOUND IN THE INPUT DIR |
#-----------------------------+
if ($count_files == 0) {
    print STDERR "\a";
    print STDERR "\nERROR: No fasta files were found in the input directory\n".
	"$indir\n".
	"Fasta files must have the fasta or fa extension.\n\n";
    exit;
}

print STDERR "NUMBER OF FILES TO PROCESS: $count_files\n" if $verbose;

#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+
for my $ind_file (@fasta_files) {


    my $name_root;
    
    #-----------------------------+
    # Get the root name of the    |
    # file to mask                |
    #-----------------------------+
    if ($ind_file =~ m/(.*)\.masked\.fasta$/) {
	# file ends in .masked.fasta
	$name_root = "$1";
    }
    elsif ($ind_file =~ m/(.*)\.hard\.fasta$/) {
	# file ends in .masked.fasta
	$name_root = "$1";
    }
    elsif ($ind_file =~ m/(.*)\.fasta$/ ) {	    
	# file ends in .fasta
	$name_root = "$1";
    }  
    elsif ($ind_file =~ m/(.*)\.fa$/ ) {	    
	# file ends in .fa
	$name_root = "$1";
    } 
    else {
	$name_root = $ind_file;
    }
    
    #-----------------------------+
    # CREATE ROOT NAME DIR        |
    #-----------------------------+
    my $name_root_dir = $outdir.$name_root."/";
    unless (-e $name_root_dir) {
	mkdir $name_root_dir ||
	    die "Could not create dir:\n$name_root_dir\n";
    }
    
    #-----------------------------+
    # CREATE THE GENEID DIR       |
    #-----------------------------+
    my $geneid_dir = $name_root_dir."geneid/";
    unless (-e $geneid_dir) {
	mkdir $geneid_dir ||
	    die "Could not create genscan out dir:\n$geneid_dir\n";
    }

    #-----------------------------+
    # CREATE GFF OUTDIR           |
    #-----------------------------+
    # This will hold the gff file modified from the
    # original gff fine
    my $gff_dir = $name_root_dir."gff/";
    unless (-e $gff_dir) {
	mkdir $gff_dir ||
	    die "Could not create genscan out dir:\n$gff_dir\n";
    }


    #-----------------------------+
    # RUN THE GENID PROGRAM       |
    #-----------------------------+
    my $geneid_out = $geneid_dir.$name_root.".geneid.gff";

    my $geneid_cmd = "$geneid_path -G -P $geneid_param_file ".
	$indir.$ind_file.
	" > ".$geneid_out;
    
    print STDERR "$geneid_cmd\n" if $verbose;
    
    system( $geneid_cmd ) unless $do_test;


    #-----------------------------+
    # CONVERT GENEID GFF TO       |
    # APOLLO TOLERATED GFF        |
    #-----------------------------+
    if (-e $geneid_out) {
	my $gff_outfile = $gff_dir.$name_root.".geneid.gff";
	geneid2gff ("geneid", $geneid_out, $gff_outfile);

    }


}


exit 0;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub geneid2gff {

    my ($source, $geneid_output, $gff_out, $seq_id, $src_suffix, $do_append) 
	=  @_;
    
    # Array to hold all of the geneid results
    my @geneid_results;

    my $attribute;

    #-----------------------------+
    # OPEN THE GFF OUTFILE        |
    #-----------------------------+
     # Default to STDOUT if no argument given
    if ($gff_out) {
	if ($do_append) {
	    open (GFFOUT, ">>$gff_out") ||
		die "ERROR: Can not open gff outfile:\n $gff_out\n";
	}
	else {
	    open (GFFOUT,">$gff_out") ||
		die "ERROR: Can not open gff outfile:\n $gff_out\n";
	    if ($gff_ver =~ "GFF3") {
		print GFFOUT "##gff-version 3\n";
	    }
	    
	}
    } 
    else {
	open (GFFOUT, ">&STDOUT") ||
	    die "Can not print to STDOUT\n";
	if ($gff_ver =~ "GFF3") {
	    print GFFOUT "##gff-version 3\n";
	}
	
    }

    #-----------------------------+
    # SET PROGRAM SOURCE          |
    #-----------------------------+
    unless ($source) {
	$source = "geneid";
    }
    if ($src_suffix) {
	$source = $source.":".$src_suffix;
    }

    open (GENEID, "<$geneid_output") ||
	die "Can not open genid gff file";

    #-----------------------------------------------------------+
    # GET DATA FROM THE GENEID OUTPUT FILE                      |
    #-----------------------------------------------------------+
    my $gene_strand;
    my $i = -1;         # indexing gene count
    my $j;              # indexing exon count for individual gene
    while (<GENEID>) {
	
	#next if m/^\#/;
	# May also want to choose to ignore comment lines ...
	
	# Parsing line in form of 
	# # Gene 1 (Forward). 20 exons. 1122 aa. Score = 63.52 
	if ( m/^\#/) {
	    if (m/\# Gene (\d*) \((.*)\)\. (\d*) exons\. (\d*) aa\. Score \= (\d*)\.(\d*)/) {
		$i++;
		$j=-1;
		my $geneid_num = $1;
		my $geneid_dir = $2;
		my $geneid_num_exons = $3;
		my $geneid_num_aa = $4;
		my $geneid_score = $5.".".$6;
		
		$geneid_results[$i]{num} = $1;
		$geneid_results[$i]{dir} = $2;
		$geneid_results[$i]{num_exons} = $3;
		$geneid_results[$i]{num_aa} = $4;
		$geneid_results[$i]{gene_score} = $5.".".$6;
		
		
		#-----------------------------+
		# GET GENE STRAND             |
		#-----------------------------+
		if ($geneid_dir =~ "Forward") {
		    $geneid_results[$i]{gene_strand} = "+";
		}
		elsif ($geneid_dir =~ "Reverse") {
		    $geneid_results[$i]{gene_strand} = "-";
		}
		else {
		    $geneid_results[$i]{gene_strand} = "?";
		}
		
#		# Test variables passed to strings
#		print STDERR "///////////////////////////////////////\n";
#		print STDERR "Gene ".$geneid_num."\t";
#		print STDERR "Dir: ".$geneid_dir."\t";
#		print STDERR "Exons: ".$geneid_num_exons."\t";
#		print STDERR "AA: ".$geneid_num_aa."\t";
#		print STDERR "Score: ".$geneid_score."\t";
#		print STDERR "Strand:".$gene_strand."\t";
#		print STDERR "\n";
#		print STDERR "///////////////////////////////////////\n";
#		
#		# test variables passed to array of hashes
#		print STDERR $geneid_results[$i]{num}."\t";
#		print STDERR $geneid_results[$i]{dir}."\t";
#		print STDERR $geneid_results[$i]{num_exons}."\t";
#		print STDERR $geneid_results[$i]{num_aa}."\t";
#		print STDERR $geneid_results[$i]{gene_score}."\t";
#		print STDERR "\n";
#		print STDERR "///////////////////////////////////////\n";
		
	    }
	}
	# parsing of exon information
	else {
	    # Split GFF results by tab
	    if (m/(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)/) {
		$j++;
		if ($j ==0) {
		    $geneid_results[$i]{gene_start} = $4;
		}

		$geneid_results[$i]{exon}[$j]{seq_id} = $1;
		$geneid_results[$i]{exon}[$j]{prog_source} = $2;
		$geneid_results[$i]{exon}[$j]{feat_type} = $3;
		$geneid_results[$i]{exon}[$j]{start} = $4;
		$geneid_results[$i]{exon}[$j]{end} = $5;
		$geneid_results[$i]{exon}[$j]{score} = $6;
		$geneid_results[$i]{exon}[$j]{strand} = $7;
		$geneid_results[$i]{exon}[$j]{frame} = $8;
		$geneid_results[$i]{exon}[$j]{attribute} = $9;		
		# keep overwriting until correct
		$geneid_results[$i]{gene_end} = $5;
		$geneid_results[$i]{gene_attribute} = $9;
		$geneid_results[$i]{gene_seq_id} = $1;
	    }
	    else {
		print STDERR "Warning:Line not in expected format:\n".$_."\n";
	    }
	}

    } # End of while GENEID

    #-----------------------------------------------------------+
    # PRINT OUTPUT FROM THE ARRAY OF HASHES
    #-----------------------------------------------------------+
    my $parent_id;
    for my $href ( @geneid_results ) {

	# PRINT GENE SPAN INFORMATION
	
	#-----------------------------+
	# PRINT GENE VALS FOR GFF3    |
	#-----------------------------+
	if ($gff_ver =~ "GFF3") {


	    $parent_id = $href->{gene_attribute};
	    unless ($seq_id) {
		$seq_id = $href->{gene_seq_id};
	    }

	    print GFFOUT $seq_id."\t".                # seq id
		"geneid\t".
		"gene\t".
		$href->{gene_start}."\t".    # start
		$href->{gene_end}."\t".      # end
		$href->{gene_score}."\t".    # score
		$href->{gene_strand}."\t".        # strand
		".\t".                       # Frame
		"ID=".$parent_id."\t".      # attribute
		"\n";

	}

	#-----------------------------+
	# EXONS
	#-----------------------------+
	my $exon_count = 0;
	for my $ex ( @{ $href->{exon} } ) {
	    $exon_count++;

	    

	    if ($gff_ver =~ "GFF3") {
	    $attribute = "ID=".$ex->{attribute}."_exon_".$exon_count.
		";Parent=".$parent_id;
	    }
	    else {
		$attribute = $ex->{attribute};
	    }

	    print GFFOUT $ex->{seq_id}."\t".
		"geneid\t".
		"exon\t".
		$ex->{start}."\t".
		$ex->{end}."\t".
		$ex->{score}."\t".
		$ex->{strand}."\t".
		$ex->{frame}."\t".
		$attribute."\t".
		"\n";


	}

    }

    # PRINT TO GFF OUT
    
    close GENEID;
    close GFFOUT;
    
    
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


sub seqid_encode {
    # Following conventions for GFF3 v given at http://gmod.org/wiki/GFF3
    # Modified from code for urlencode in the perl cookbook
    # Ids must not contain unescaped white space, so spaces are not allowed
    my ($value) = @_;
    $value =~ s/([^[a-zA-Z0-9.:^*$@!+_?-|])/"%" . uc(sprintf "%lx" , unpack("C", $1))/eg;
    return ($value);
}

sub gff3_encode {
    # spaces are allowed in attribute, but tabs must be escaped
    my ($value) = @_;
    $value =~ s/([^[a-zA-Z0-9.:^*$@!+_?-| ])/"%" . uc(sprintf "%lx" , unpack("C", $1))/eg;
    return ($value);
}


1;
__END__

=head1 NAME

batch_geneid.pl - Run the genid program in batch mode.

=head1 VERSION

This documentation refers to program version $Rev: 948 $

=head1 SYNOPSIS

=head2 Usage

    batch_geneid.pl -i InDir -o OutDir

=head2 Required Arguments

    --infile        # Path to the input directory
    --outfie        # Path to the output directory
    --param         # Path to the GeneID parameter file

=head1 DESCRIPTION

This program will run the geneid gene annotation program on a 
directory of fasta files.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path of the directory containing the sequences to process.

=item -o,--outdir

Path of the directory to place the program output.

=item -p, --param

Path to the GeneID Parameter file.

=back

=head1 OPTIONS

=over 2

=item --gff-ver

The GFF version for the output. This will accept either gff2 or gff3 as the
options. By default the GFF version will be GFF2 unless specified otherwise.
The default GFF version for output can also be set in the user environment
with the DP_GFF option. The command line option will always override the option
defined in the user environment. 

=item --geneid-path

The path to the geneid program. If a path is not supplied at the command
line the program will assume that the directory containing the geneid
program is included in the user PATH parameter.

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

=item -q,--quiet

Run the program with minimal output.

=item --test

Run the program without doing the system commands. This will
test for the existence of input files.

=back

=head1 DIAGNOSTICS

Error messages generated by this program and possible solutions are listed
below.

=over 2

=item ERROR: No fasta files were found in the input directory

The input directory does not contain fasta files in the expected format.
This could happen because you gave an incorrect path or because your sequence 
files do not have the expected *.fasta extension in the file name.

=item ERROR: Could not create the output directory

The output directory could not be created at the path you specified. 
This could be do to the fact that the directory that you are trying
to place your base directory in does not exist, or because you do not
have write permission to the directory you want to place your file in.

=back

=head1 CONFIGURATION AND ENVIRONMENT

=head2 Configuration File

This program currently does not make use of a configuration file.

=head2 Environment

This program does not make use of variables in the user environment.

=head1 DEPENDENCIES

=head2 Required Software

=over

=item * geneid

This program requies the geneid gene annotation program to work. This program
is available from:
http://genome.imim.es/software/geneid/index.html

=back

=head2 Required Perl Modules

=over

=item * Getopt::Long

This module is required to accept options at the command line.

=back

=head1 BUGS AND LIMITATIONS

Any known bugs and limitations will be listed here.

=head2 Bugs

=over 2

=item * No bugs currently known 

If you find a bug with this software, file a bug report on the DAWG-PAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=204962

=back

=head2 Limitations

=over

=item * Limited testing on versions of geneid

This program is known to work on genid version March_1_2005. This program
has not been tested with other versions of geneid.

=back

=head1 SEE ALSO

The program is part of the DAWG-PAWS package of genome
annotation programs. See the DAWG-PAWS web page 
( http://dawgpaws.sourceforge.net/ )
or the Sourceforge project page 
( http://sourceforge.net/projects/dawgpaws ) 
for additional information about this package.

=head1 REFERENCE

A manuscript is being submitted describing the DAWGPAWS program. 
Until this manuscript is published, please refer to the DAWGPAWS 
SourceForge website when describing your use of this program:

JC Estill and JL Bennetzen. 2009. 
The DAWGPAWS Pipeline for the Annotation of Genes and Transposable 
Elements in Plant Genomes.
http://dawgpaws.sourceforge.net/

=head1 LICENSE

GNU General Public License, Version 3

L<http://www.gnu.org/licenses/gpl.html>

THIS SOFTWARE COMES AS IS, WITHOUT ANY EXPRESS OR IMPLIED
WARRANTY. USE AT YOUR OWN RISK.

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 11/06/2008

UPDATED: 03/11/2010

VERSION: $Rev: 948 $

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
# 11/06/2008
# - Main body of program written
# 01/19/2009
# - Added Rev to SVN propset
# 03/12/2010
# - Testing swtich to Bio::Tools::Geneid for parsing
# - This did not work, erros like
#Argument "Internal" isn't numeric in sort at /Users/jestill/code/bioperl-live//Bio/RangeI.pm line 412, <GEN2> line 24.
#Argument "geneid_v1.3" isn't numeric in sort at /Users/jestill/code/bioperl-live//Bio/RangeI.pm line 410, <GEN2> line 25.
#Argument "Internal" isn't numeric in sort at /Users/jestill/code/bioperl-live//Bio/RangeI.pm line 412, <GEN2> line 25.
#Argument "geneid_v1.3" isn't numeric in sort at /Users/jestill/code/bioperl-live//Bio/RangeI.pm line 410, <GEN2> line 26.
#Argument "Internal" isn't numeric in sort at /Users/jestill/code/bioperl-live//Bio/RangeI.pm line 412, <GEN2> line 26.
