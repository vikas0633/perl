#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_genscan.pl - Run genscan gene prediction program    |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill at gmail.com                         |
# STARTED: 07/31/2007                                       |
# UPDATED: 04/28/2009                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Run the genscan gene prediction program in batch mode.   |
#  Runs genscan as well as converts output to gff format.   |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+
# minor modification

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;
use File::Copy;
use Getopt::Long;
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

#-----------------------------+
# VARS USING ENV              |
#-----------------------------+
# If not provided here, use the default assumtion that
# the file location are in the user's path and will
# use Maize as the default library.
my $genscan_path = $ENV{DP_GENSCAN_BIN} ||
    "genscan";
my $lib_path = $ENV{DP_GENSCAN_LIB} ||
    "Maize.smat";

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $logfile;                   # Path to a logfile to log error info
my $indir;                     # Directory containing the seq files to process
my $outdir;                    # Directory to hold the output
my $msg;                       # Message printed to the log file

my $search_name;               # Name searched for in grep command
my $bac_out_dir;               # Dir for each sequnce being masked
my $name_root;                 # Root name to be used for output etc

# BOOLEANS
my $show_help = 0;             # Show program help
my $show_version = 0;          # Show program version
my $show_man = 0;              # Show program manual page using peldoc
my $show_usage = 0;            # Show program usage command             
my $quiet = 0;                 # Boolean for reduced output to STOUT
my $apollo = 0;                # Path to apollo and apollo variables
my $test = 0;
my $verbose = 0;

# COUNTERS
my $num_proc = 1;              # Number of processors to use
my $param;
my $program;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(
		    # REQUIRED
		    "i|indir=s"    => \$indir,
                    "o|outdir=s"   => \$outdir,
		    # OPTIONS
		    "param=s"      => \$param,
		    "program=s"    => \$program,
		    "lib-path"     => \$lib_path,
		    "gff-ver=s"    => \$gff_ver,
		    "genscan-path" => \$genscan_path,
		    "logfile=s"    => \$logfile,
		    "apollo"       => \$apollo,
		    # ADDITIONAL INFORMATION
		    "verbose"      => \$verbose,
		    "test"         => \$test,
		    "usage"        => \$show_usage,
		    "version"      => \$show_version,
		    "man"          => \$show_man,
		    "h|help"       => \$show_help,
		    "q|quiet"      => \$quiet,);

my $proc_num = 0;

#//////////////////////
my $file_num_max = 2;
my $file_num = 0;
#\\\\\\\\\\\\\\\\\\\\\\

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
    print "\nbatch_genscan.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}

#-----------------------------+
# CHECK REQUIRED ARGS         |
#-----------------------------+
if ( (!$indir) || (!$outdir) ) {
    print "\a";
    print STDERR "\n";
    print STDERR "ERROR: An input directory was not specified at the".
	" command line\n" if (!$indir);
    print STDERR "ERROR: An output directory was specified at the".
	" command line\n" if (!$outdir);
    print_help ("usage", $0 );
}

#-----------------------------+
# OPEN THE LOG FILE           |
#-----------------------------+
#if ($logfile) {
#    # Open file for appending
#    open ( LOG, ">>$logfile" ) ||
#	die "Can not open logfile:\n$logfile\n";
#    my $time_now = time;
#    print LOG "==================================\n";
#    print LOG "  batch_genscan.pl\n";
#    print LOG "  JOB: $time_now\n";
#    print LOG "==================================\n";
#}

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
my $num_files = @fasta_files;

#-----------------------------+
# SHOW ERROR IF NO FILES      |
# WERE FOUND IN THE INPUT DIR |
#-----------------------------+
if ($num_files == 0) {
    print "\a";
    print "\nERROR: No fasta files were found in the input directory\n".
	"$indir\n".
	"Fasta files must have the fasta or fa extension.\n\n";
    exit;
}

#-----------------------------+
# CREATE THE OUT DIR          |
# IF IT DOES NOT EXIST        |
#-----------------------------+
unless (-e $outdir) {
    print "Creating output dir ...\n" unless $quiet;
    mkdir $outdir ||
	die "Could not create the output directory:\n$outdir";
}

#-----------------------------+
# RUN GENSCAN AND PARSE       |
# RESULTS FOR EACH SEQ IN THE |
# fasta_files ARRAY           |
#-----------------------------+

for my $ind_file (@fasta_files) {
    
    $proc_num++;
    $file_num++;
    #if ($file_num == $file_num_max){exit;}

    #-----------------------------+
    # GET THE ROOT NAME OF THE    |
    # FASTA FILE                  |
    #-----------------------------+
    if ($ind_file =~ m/(.*)\.hard\.fasta$/) {
	# file ends in .hard.fasta
	# This is hard masked fasta files
	$name_root = "$1";
    }
    elsif ($ind_file =~ m/(.*)\.masked\.fasta$/) {
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

    my $infile_path = $indir.$ind_file;
    

    #-----------------------------+
    # CREATE ROOT NAME DIR        |
    #-----------------------------+
    my $name_root_dir = $outdir.$name_root;
    unless (-e $name_root_dir) {
	mkdir $name_root_dir ||
	    die "Could not create dir:\n$name_root_dir\n"
    }

    #-----------------------------+
    # CREATE GENSCAN OUTDIR       |
    #-----------------------------+
    # Dir to hold gene prediction output from local software
    my $genscan_dir = $outdir.$name_root."/gene/";
    unless (-e $genscan_dir) {
	mkdir $genscan_dir ||
	    die "Could not create genscan out dir:\n$genscan_dir\n";
    }
    
    my $gff_dir = $outdir.$name_root."/gff/";
    unless (-e $gff_dir) {
	mkdir $gff_dir ||
	    die "Could not create gff out dir:\n$gff_dir\n";
    }

    my $out_path = $genscan_dir.$name_root.".genscan.out";
    my $gff_path = $gff_dir.$name_root.".genscan.gff";

    my $genscan_cmd = "$genscan_path $lib_path $infile_path -v > $out_path";

    print STDERR "=======================================\n" if $verbose;
    print STDERR "Running Genscan for $name_root\n" if $verbose;
    print STDERR " File $file_num of $num_files\n" if $verbose;
    print STDERR "=======================================\n" if $verbose;
    print STDERR "$genscan_cmd\n" if $verbose;

    system($genscan_cmd) unless $test;

    print STDERR "Converting Genscan output\n";
    if (-e $out_path) {
	genscan2gff($out_path, $gff_path, $name_root, $program, $param);
    }
    else {
	print STDERR "ERROR: Could not find genscan output at:\n$out_path\n"
    }
    
} # End of for each file in the input folder

#close LOG if $logfile;

exit;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+
sub genscan2gff {

    my ($infile, $outfile, $name_root, $dp_source, $dp_source_suffix) = @_;


    # Encode sequence name to be legal
    $name_root = seqid_encode ($name_root);

    my $gff_source;
    unless ($dp_source) {
	$gff_source = "genscan";
    }
    else {
	$gff_source = $dp_source;
    }
    
    if ($dp_source_suffix) {
	$gff_source = $dp_source.":".$dp_source_suffix;
    }

    my @genscan_results;

    my %exon_type = ('Sngl', 'Single Exon',
		     'Init', 'Initial Exon',
		     'Intr', 'Internal Exon',
		     'Term', 'Terminal Exon');
    
    if ($infile) {
	open (IN, "<".$infile) || 
	    die "Can not open genscan file for input:\n$infile\n";
    }
    else {
	open (IN, "<&STDIN") || 
	    die "Can not STDIN for input\n";
    }


    if ($outfile) {
	open (GFFOUT, ">".$outfile) ||
	    die "Can not open GFF putout file for output:\n$outfile";
    } 
    else {
	open (GFFOUT, ">&STDOUT") ||
	    die "Can not print to STDOUT\n";
    }

    # PRINT GFF VERSION HEADER
    if ($gff_ver =~ "GFF3") {
	print GFFOUT "##gff-version 3\n";
    }

    while (<IN>) {

	# The following appears to just fetch the exons
	
	# Last line before predictions contains nothing but spaces and dashes
	if (/^\s*-[-\s]+$/)  {

	    my $cur_gene_name;
	    my $prev_gene_name = "null";
	    my $exon_count = 0;
	    my $i = -1;                     # i indexes gene count
	    my $j = -1;                     # j indexes exon count

	    while (<IN>)  {
		#my %feature;

		# TO DO: Add use of polyy and promoter
		#        at the moment and promoter at the moment
		if (/init|term|sngl|intr/i) {
		    
		    my @f  = split;
		    
		    my ($gene, $exon) = split (/\./, $f[0]); 
		    my $cur_gene_name = $gene;

		    #-----------------------------+
		    # PUT START < END             |
		    #-----------------------------+
		    my $start;
		    my $end;
		    if ($f[2] eq '+') {
			$start  = $f[3];
			$end = $f[4];
		    } elsif ($f[2] eq '-') {
			$start  = $f[4];
			$end = $f[3];
		    }


		    if ($cur_gene_name =~ $prev_gene_name) {
			# IN SAME GENE MODEL
			$j++;  # increment exon count
		    } else {
			# IN NEW GENE MODEL
			$j=-1;
			$i++;   # increment gene count
			$j++;   # increment exon count

			$genscan_results[$i]{gene_strand} = $f[2];
			$genscan_results[$i]{gene_name} = $gff_source."_".
			    "gene_".
			    $cur_gene_name;
			$genscan_results[$i]{gene_start} = $start;

		    }
		    # set previous name to current name after comparison
		    $prev_gene_name = $cur_gene_name;

		    # Continue to overwrite end
		    $genscan_results[$i]{gene_end} = $end;

		    # LOAD EXON INFORMATION TO ARRAY
		    $genscan_results[$i]{exon}[$j]{exon_id} = $exon;
		    $genscan_results[$i]{exon}[$j]{start} = $start;
		    $genscan_results[$i]{exon}[$j]{end} = $end;
		    $genscan_results[$i]{exon}[$j]{score} = $f[12];
		    $genscan_results[$i]{exon}[$j]{strand} = $f[2];
		    # Probability of exon
		    $genscan_results[$i]{exon}[$j]{p} = $f[11];
		    # Coding region score
		    $genscan_results[$i]{exon}[$j]{cod_rg} = $f[10];
		    $genscan_results[$i]{exon}[$j]{exon_type} = 
			$exon_type{$f[1]};
		    
		} elsif (/predicted peptide/i) {
		    last;   
		}
	    } # End of second while statement
	} # End of if seach command
    } # End of while INPUT


   #-----------------------------+
   # PRINT GFF OUTPUT            |
   #-----------------------------+
    my $parent_id;
    for my $href ( @genscan_results ) {
	
	# If GFF3 need to print the parent gene span
	if ($gff_ver =~ "GFF3") {
	    $parent_id = $href->{gene_name};

	    print GFFOUT $name_root."\t".                # seq id
		$gff_source."\t".
		"gene\t".
		$href->{gene_start}."\t".    # start
		$href->{gene_end}."\t".      # end
		".\t".    # score
		$href->{gene_strand}."\t".        # strand
		".\t".                       # Frame
		"ID=".$parent_id."\t".      # attribute
		"\n";

	}

	#-----------------------------+
	# EXONS
	#-----------------------------+
	my $exon_count = 0;
	my $attribute;
	for my $ex ( @{ $href->{exon} } ) {

	    $exon_count++;
	    
	    if ($gff_ver =~ "GFF3") {
		$attribute = "ID=".$href->{gene_name}.
		    "_exon_".
		    $ex->{exon_id}.
		    ";Parent=".$parent_id;
	    }
	    else {
		$attribute = $href->{gene_name};
	    }
	    
	    print GFFOUT $name_root."\t".
		$gff_source."\t".
		"exon\t".
		$ex->{start}."\t".
		$ex->{end}."\t".
		$ex->{score}."\t".
		$ex->{strand}."\t".
		".\t".
		$attribute."\t".
		"\n";

	    
	}
	

	
    }

} #End of genscan_2_gff subfunction


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

batch_genscan.pl - Run genscan and parse results to a gff format file. 

=head1 VERSION

This documentation refers to batch_genscan.pl version $Rev: 948 $

=head1 SYNOPSIS

=head2 Usage

    batch_genscan.pl -i DirToProcess -o OutDir

=head2 Required Variables

    -i, --indir    # Directory of fasta files to process
    -o, --outdir   # Path to the base output directory

=head1 DESCRIPTION

Run the GENSCAN gene prediction program in batch mode.
This will run genscan as well as convert the output to gff format.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path of the directory containing the sequences to process.

=item -o,--outdir

Path of the directory to place the program output.

=back

=head1 OPTIONS

=over 2

=item --genscan-path

The full path to the genscan binary.

=item --lib-path

The full path to the library file.

=item --logfile

Path to a file that will be used to log program status.
If the file already exists, additional information will be concatenated
to the existing file.

=item --usage

Short overview of how to use program from command line.

=item --help

Show program usage with summary of options.

=item --version

Show program version.

=item --man

Show the full program manual. This uses the perldoc command to print the 
POD documentation for the program.

=item -q,--quiet

Run the program with minimal output.

=item --test

Run the program without doing the system commands.

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

The batch_genscan.pl program currently does not make use of any 
external configuration files.

=head2 User Environment

The following environment variables may be used to designate the full
path to the GENSCAN binary and the default library path to use.

=over 

=item DP_GENSCAN_BIN

The path to the genscan binary file.

=item DP_GENSCAN_LIB

The path to the default library to use for gene predictions using genscan.

=back

The following example shows the lines that would need to be added to your
.bashrc file in the bash shell to specify the location of the genscan
binary and the Maize library file in your home directory.

    export DP_GENSCAN_BIN='/home/yourname/apps/genscan/genscan'
    export DP_GENSCAN_LIB='/home/yourname/apps/genscan/Maize.smat'

=head1 DEPENDENCIES

=head2 Required Software

=over

=item * GENSCAN

This program obviously requires the GENSCAN gene model prediction program.
Information for downloading the GENSCAN program is available from:
http://genes.mit.edu/GENSCANinfo.html

=back

=head2 Required Perl Modules

=over

=item * File::Copy

This module is required to copy the results.

=item * Getopt::Long

This module is required to accept options at the command line.

=back

=head1 BUGS AND LIMITATIONS

=head2 Bugs

=over 2

=item * No bugs currently known 

If you find a bug with this software, file a bug report on the DAWG-PAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=204962

=back

=head2 Limitations

=over

=item * Config file must use UNIX format line endings

The config file must have UNIX formatted line endings. Because of
this any config files that have been edited in programs such as
MS Word must be converted to a UNIX compatible text format before
being used with batch_blast.

=back

=head1 SEE ALSO

The batch_genscan.pl program is part of the DAWG-PAWS package of genome
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

GNU GENERAL PUBLIC LICENSE, VERSION 3

http://www.gnu.org/licenses/gpl.html  

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 07/31/2007

UPDATED: 03/15/2010

VERSION: $Rev: 948 $

=cut

#-------------------------------------------------------+
# HISTORY                                               |
#-------------------------------------------------------+
#
# 07/31/2007
# - Program started and bulk of functions added and
#   working copy of the program finished
# 
# 12/11/2007
# - Added SVN tracking of Rev
# - Moved POD documentation to the end of the file
# - Added use strict
# - Added print_help subfunction that extracts help
#   and usage statements from the POD documentation
# - Update POD documentation
# - Changed the path to the genscan binary and the default
#   libary to variables set in the user environment
#
# - Changed Version reference to $0 instead of hard coded
#   program name
# - Additional updates to POD documentation
#
# 04/28/2009
# - Putting output in gff dir
#
# 03/15/2010
# - Adding support for GFF3 output
# - Loading data to array before printing to GFF
# - Dropping LOG file, using STDOUT instead
# - Rewrote subfunction for GFF conversion
#
# TO DO: 
# Switch to a config file of
# LIBNAME tab   LIBPATH
