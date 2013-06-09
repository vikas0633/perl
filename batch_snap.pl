#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_snap.pl - Run SNAP gene prediction in batch mode    |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 03/29/2010                                       |
# UPDATED: 04/08/2010                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Run the SNAP gene predition program in batch mode.       |
#                                                           |
# USAGE:                                                    |
#  batch_snap.pl -i indir -o outdir -c snap_config.cfg      |
#                                                           |
# VERSION: $Rev$                                            |
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

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/;
# Get GFF version from environment, GFF2 is DEFAULT
my $gff_ver = uc($ENV{DP_GFF}) || "GFF2";

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $indir;                     # Files to process
my $outdir;                    # Base output dir
my $config_file;               # Configuration file
my @snap_params = ();          # SNAP Parameters

# BOOLEANS
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_test = 0;                  # Run the program in test mode


# Path to the snap program
my $sn_path = $ENV{SNAP_PATH} ||
    "snap";
my $program = "SNAP";
my $param;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|indir=s"   => \$indir,
                    "o|outdir=s"  => \$outdir,
                    "c|config=s"  => \$config_file,
		    # ADDITIONAL OPTIONS
		    "snap-bin"    => \$sn_path,
		    "program=s"   => \$program,
		    "gff-ver=s"   => \$gff_ver,
		    "q|quiet"     => \$quiet,
		    "verbose"     => \$verbose,
		    # ADDITIONAL INFORMATION
		    "usage"       => \$show_usage,
		    "test"        => \$do_test,
		    "version"     => \$show_version,
		    "man"         => \$show_man,
		    "h|help"      => \$show_help,);

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
    print "\nbatch_mask.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}

#-----------------------------+
# CHECK REQUIRED ARGS         |
#-----------------------------+
if ( (!$indir) || (!$outdir) || (!$config_file) ) {
    print "\a";
    print STDERR "\n";
    print STDERR "ERROR: An input directory was not specified at the".
	" command line\n" if (!$indir);
    print STDERR "ERROR: An output directory was not specified at the".
	" command line\n" if (!$outdir);
    print STDERR "ERROR: A configuration file was not specified at the".
	" command line\n" if (!$config_file);
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

#-----------------------------+
# PARSE THE CONFIG FILE       |
#-----------------------------+
my $i=0;
my $config_line_num=0;

open (CONFIG, "<$config_file") ||
    die "ERROR Can not open the config file:\n $config_file";

while (<CONFIG>) {
    chomp;
    $config_line_num++;
    unless (m/^\#/) {
       	my @in_line = split (/\t/);           # Implicit split of $_ by tab
	my $num_in_line = @in_line; 
	
	# Can have just a name to run LTR_Finder with default settings
	# or can have two columns with additional parameter options
	# parameter options in the second columns.
	# I will currently stick with the two column config file
	# since there are so many options availabe with LTR_FINDER
	# that a multiple column config file would get messy.
	if ($num_in_line == 2 ||
	    $num_in_line == 3) { 
	    $snap_params[$i][0] = $in_line[0];            # Name
	    $snap_params[$i][1] = $in_line[1];            # HMM File Path
	    $snap_params[$i][2] = $in_line[2] || "NULL";  # Options
	    $i++;
	} # End of if $num_in_line is 10
	else {
	    print "\a";
	    print STDERR "WARNING: Unexpected number of line in config".
		" file line $config_line_num\n$config_file\n";
	}

   } # End of unless comment line
} # End of while CONFIG file
close CONFIG;

# Number of parameter sets specified in the config file
my $num_par_sets = $i;

if ($num_par_sets == 0) {
    print "\a";
    print STDERR "ERROR: No parameter sets were found in the config file:\n".
	"$config_file\n";
    exit;
}

my $num_proc_total = $count_files * $num_par_sets;

print STDERR "$num_proc_total find_ltr runs to process\n\n" if $verbose;

#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+
my $file_num = 0;
my $proc_num = 0;
for my $ind_file (@fasta_files) {

    $file_num++;
    my $name_root;
    # Get root file name
    if ($ind_file =~ m/(.*)\.masked\.fasta$/) {
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
	    die "Could not create dir:\n$name_root_dir\n"
    }

    #-----------------------------+
    # CREATE SNAP OUTDIR           |
    #-----------------------------+
    # Dir to hold gene prediction output from local software
    my $snap_dir = $name_root_dir."snap/";
    unless (-e $snap_dir) {
	mkdir $snap_dir ||
	    die "Could not create ltr_finder out dir:\n$snap_dir\n";
    }
    
    #-----------------------------+
    # CREATE GFF OUTDIR           |
    #-----------------------------+
    my $gff_dir = $name_root_dir."gff/";
    unless (-e $gff_dir) {
	mkdir $gff_dir ||
	    die "Could not create gff out dir:\n$gff_dir\n";
    }




    #-----------------------------+
    # FOR EACH PARAM SET  IN THE  |
    # CONFIG FILE                 |
    #-----------------------------+
    my $gff_count=0;               # Count of the number of gff files

    for ($i=0; $i<$num_par_sets; $i++) {	
	
	$proc_num++;
	my $sn_param_name = $snap_params[$i][0];
	my $sn_matrix = $snap_params[$i][1];
	my $sn_opt = chomp($snap_params[$i][2]) 
	    || "NULL";
	
	my $snap_in = $indir.$ind_file;
	# The SNAP result file
	my $sn_res = $snap_dir.$name_root."_snap_".
	    $sn_param_name.".txt";
	# SNAP result in GFF format
	my $sn_gff = $gff_dir.$name_root."_snap_".
	    $sn_param_name.".gff";
	
	# RUN THE SNAP PROGRAM
	my $sn_cmd = $sn_path;
	$sn_cmd = $sn_path." -gff";
	unless ($sn_opt =~ "NULL") {
	    $sn_cmd = $sn_cmd." ".$sn_opt
	}
	$sn_cmd = $sn_cmd." ".$sn_matrix." ".
	    $snap_in." > ".$sn_res;
	
	print STDERR "CMD: ".$sn_cmd."\n" if $verbose;

	# CONVERT THE SNAP OUTPUT TO GFF
	system ($sn_cmd) unless $do_test;

	if (-e $sn_res) {
	    snap2gff ($program, $sn_res, $sn_gff, $name_root, 
		      $sn_param_name, 0);
	}

    }
    # Running the SNAP program
    # name         - the base name to use for output
    # O.sativa.hmm - the parametr file to use
    # HEX*fasta    - the sequence file to annotate
    # gff          - Produce GFF formatted output
    #                This is GFF3
    #
    #snap -name snap_def -gff O.sativa.hmm HEX3045G05.masked.fasta 



}





exit 0;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub snap2gff {
    
    my ($source, $snap_in, $gffout, $seq_id, $src_suffix, $do_append ) = @_;
    # Array to hold all snap results for a single contig
    my @snap_results;
    my $prev_gene_name;

    # Rows are the individual genes
    #   exons array nested in snap results will contain results
    #   for each individual gene model
    # Starting i and j at -1 so that increments start at 0
    my $i = -1;
    my $j = -1;
    my $model_num = 0;

    if ($src_suffix) {
	$source = $source.":".$src_suffix;
    }
    
    my $attribute;

    #-----------------------------+
    # OPEN INPUT FILE HANDLE      |
    #-----------------------------+
    if ($snap_in) {
	open (INFILE, "<$snap_in") ||
	    die "ERROR: Can not open LTR_FINDER result file\n $snap_in\n";
	
    }
    else {
	print STDERR "Expecting input from STDIN\n";
	open (INFILE, "<&STDIN") ||
	    die "Can not accept input from standard input.\n";
    }
    
    #-----------------------------+
    # OPEN GFF OUTFILE            |
    #-----------------------------+
    # Default to STDOUT if no argument given
    if ($gffout) {
	if ($do_append) {
	    open (GFFOUT, ">>$gffout") ||
		die "ERROR: Can not open gff outfile:\n $gffout\n";
	}
	else {
	    open (GFFOUT,">$gffout") ||
		die "ERROR: Can not open gff outfile:\n $gffout\n";
	}
    } 
    else {
	open (GFFOUT, ">&STDOUT") ||
	    die "Can not print to STDOUT\n";
    }

    if ($gff_ver =~ "GFF3") {
	print GFFOUT "##gff-version 3\n";
    }

    $prev_gene_name = "NULL";
    # PROCESS FILE
    while (<INFILE>) {
	
#	print STDERR $_;

	my @gff_parts = split;
	my $num_gff_parts = @gff_parts;
	my $exon_count = 0;
	
	my $cur_gene_name = $gff_parts[8];

	#-----------------------------+
	# MAKE START < END            |
	#-----------------------------+
	my $start;
	my $end;
	if ($gff_parts[3] < $gff_parts[4]) {
	    $start = $gff_parts[3];
	    $end = $gff_parts[4];
	} else {
	    $end = $gff_parts[4];
	    $start = $gff_parts[3];
	}


	if ($cur_gene_name =~ $prev_gene_name) {
	    # IN SAME GENE MODEL
	    $j++;  # increment exon count
	} else {
	    # IN NEW GENE MODEL
	    $j=-1;
	    $i++;   # increment gene count
	    $j++;   # increment exon count
	    
	    $snap_results[$i]{gene_strand} = $gff_parts[6];
	    $snap_results[$i]{gene_name} = $source."_".
		"gene_".
		$cur_gene_name;
	    $snap_results[$i]{gene_start} = $start;
	    $snap_results[$i]{gene_end} = $end;
	}
	
	$prev_gene_name = $cur_gene_name;
	
	# UPDATE GENE START AND END
	if ($start < $snap_results[$i]{gene_start} ) {
	    $snap_results[$i]{gene_start} = $start;
	}
	if ($end > $snap_results[$i]{gene_start} ) {
	    $snap_results[$i]{gene_end} = $end;
	}
	
	# LOAD EXON INFORMATION TO ARRAY
	if  ($seq_id) {
	    $snap_results[$i]{seq_id} = $seq_id;
	}
	else {
	    $snap_results[$i]{seq_id} = $gff_parts[0];
	}
	$snap_results[$i]{exon}[$j]{exon_id} = $j + 1;
	$snap_results[$i]{exon}[$j]{start} = $start;
	$snap_results[$i]{exon}[$j]{end} = $end;
	$snap_results[$i]{exon}[$j]{score} = $gff_parts[5];
	$snap_results[$i]{exon}[$j]{strand} = $gff_parts[6];
	$snap_results[$i]{exon}[$j]{frame} = $gff_parts[7];
	$snap_results[$i]{exon}[$j]{exon_type} = $gff_parts[2];

	print STDERR "MODEL NUMBER: $model_num\n" if $verbose;

    }

    # DONE WITH INFILE
    close INFILE;
    
    #-----------------------------+
    # PRINT GFFOUT FROM ARRAY     |
    #-----------------------------+
    my $parent_id;
    for my $href ( @snap_results ) {
		
	# If GFF3 need to print the parent gene span
	# and a fake transcript gene
	if ($gff_ver =~ "GFF3") {
	    $parent_id = $href->{gene_name};
	    
	    print GFFOUT $href->{seq_id}."\t".                # seq id
		$source."\t".
		"gene\t".
		$href->{gene_start}."\t".    # start
		$href->{gene_end}."\t".      # end
		".\t".    # score
		$href->{gene_strand}."\t".        # strand
		".\t".                       # Frame
		"ID=".$parent_id.      # attribute
		"; Name=".$parent_id.
		"\n";
	    
	    # Added fake transcript for Apollo
	    # 04/08/2010
	    # trs added fro "transcript
	    $parent_id = $href->{gene_name}."_trs";
	    print GFFOUT $href->{seq_id}."\t".                # seq id
		$source."\t".
		"transcript\t".
		$href->{gene_start}."\t".    # start
		$href->{gene_end}."\t".      # end
		".\t".    # score
		$href->{gene_strand}."\t".        # strand
		".\t".                       # Frame
		"ID=".$parent_id.      # attribute
		"; Name=".$parent_id.
		"\n";
	    
	    
	}

	my $exon_count = 0;

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
	    
	    # Currently not reporting UTRs
	    # May want to exclude UTRs from gene span reported above
		
	    print GFFOUT $href->{seq_id}."\t".                # seq id
		$source."\t".
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

    # DONE
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

1;
__END__

=head1 NAME

Name.pl - Short program description. 

=head1 VERSION

This documentation refers to program version 0.1

=head1 SYNOPSIS

  USAGE:
    Name.pl -i InDir -o OutDir

    --infile        # Path to the input file
    --outfie        # Path to the output file

=head1 DESCRIPTION

This is what the program does

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path of the directory containing the sequences to process.

=item -o,--outdir

Path of the directory to place the program output.

=item -c, --config

Path to a config file. This is a tab delimited text file
indicating the required information for each of the databases to blast
against. Lines beginning with # are ignored.

=back

=head1 OPTIONS

=over 2

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

=head1 Additional Options

=over 2

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

The location of the configuration file is indicated by the --config option
at the command line.
This is a tab delimited text file
indicating required information for each of the databases to blast
against. Lines beginning with # are ignored, and data are in six 
columns as shown below:

=over 2

=item Col 1. Blast program to use [ tblastx | blastn | blastx ]

The blastall program to use. DAWG-PAWS will support blastn,
tblastx, and blastx format.

=item Col 2. Extension to add to blast output file. (ie. bln )

This is the suffix which will be added to the end of your blast
output file. You can use this option to set different extensions
for different types of blast. For example *.bln for blastn
output and *.blx for blastx output.

=back

An example config file:

 #-----------------------------+
 # BLASTN: TIGR GIs            |
 #-----------------------------+
 blastn	bln	8	1e-5	TaGI_10	-a 2 -U
 blastn	bln	8	1e-5	AtGI_13	-a 2 -U
 blastn	bln	8	1e-5	ZmGI_17	-a 2 -U
 #-----------------------------+
 # TBLASTX: TIGR GIs           |
 #-----------------------------+
 tblastx	blx	8	1e-5	TaGI_10	-a 2 -U
 tblastx	blx	8	1e-5	AtGI_13	-a 2 -U
 tblastx	blx	8	1e-5	ZmGI_17	-a 2 -U

=head2 Environment

This program does not make use of variables in the user environment.

=head1 DEPENDENCIES

=head2 Required Software

=over

=item * Software Name

Any required software will be listed here.

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

=item * Known Limitation

If this program has known limitations they will be listed here.

=back

=head1 SEE ALSO

The program is part of the DAWG-PAWS package of genome
annotation programs. See the DAWG-PAWS web page 
( http://dawgpaws.sourceforge.net/ )
or the Sourceforge project page 
( http://sourceforge.net/projects/dawgpaws ) 
for additional information about this package.

=head1 LICENSE

GNU General Public License, Version 3

L<http://www.gnu.org/licenses/gpl.html>

THIS SOFTWARE COMES AS IS, WITHOUT ANY EXPRESS OR IMPLIED
WARRANTY. USE AT YOUR OWN RISK.

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED:

UPDATED:

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 04/08/2010 
# - The current version of Apollo balks at putting exons
#   as childrn of genes. I now need to add a fake transcript
#   model.
