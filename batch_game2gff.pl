#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_game2gff.pl - Batch convert gamexml to gff format   |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 11/01/2007                                       |
# UPDATED: 12/17/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Convert game.xml annotated feature file format to        |
#  Goal is to only export the curated results.              | 
#                                                           |
# REQUIREMENTS:                                             |
#  This requires the apollo genome annotation curation      |
#  program.                                                 |
#                                                           |
# VERSION: $Rev: 948 $                                            |
#                                                           |
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
my ($VERSION) = q$Rev: 948 $ =~ /(\d+)/;

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $indir;
my $outdir;

# Booleans
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;

# Options with default values
my $ap_path = $ENV{DP_APOLLO_BIN} || "apollo";

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|indir=s"   => \$indir,
                    "o|outdir=s"  => \$outdir,
		    # ADDITIONAL OPTIONS
		    "ap-path=s"   => \$ap_path,
		    "q|quiet"     => \$quiet,
		    "verbose"     => \$verbose,
		    # ADDITIONAL INFORMATION
		    "usage"       => \$show_usage,
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

if ($show_version) {
    print "\n$0:\nVersion: $VERSION\n\n";
    exit;
}

if ($show_man) {
    # User perldoc to generate the man documentation.
    system("perldoc $0");
    exit($ok ? 0 : 2);
}

#-----------------------------+
# CHECK REQUIRED ARGS         |
#-----------------------------+
if ( (!$indir) || (!$outdir) ) {
    print "\a";
    print "ERROR: Input directory must be specified\n" if !$indir;
    print "ERROR: Output directory must be specified\n" if !$outdir;
    print_help ("usage", $0 );
}

#-----------------------------+
# MAIN PROGRAM BODY           |
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
# Get the seq files from the  |
# directory provided by the   |
# var $indir                  |
#-----------------------------+
opendir( DIR, $indir ) || 
    die "Can't open directory:\n$indir"; 
my @seq_files = grep /\.game.xml$|\.xml$/, readdir DIR ;
closedir( DIR );

my $count_files = @seq_files;


if ($count_files == 0) {
    print "\a";
    print "ERROR: No game.xml files were found in the intput directory\n";
    print "$indir\n";
}

print "$count_files to process\n";

for my $ind_file (@seq_files)
{
    print "Processing: $ind_file\n";

    my $in_seq_path = $indir.$ind_file;
    my $out_gff_path = $outdir.$ind_file.".gff";
    my $out_fasta_path = $outdir.$ind_file.".fasta";
 
    #  Convert each the game xml to gff using the apollo_convert command
    &apollo_convert ($in_seq_path, "game", $out_gff_path, 
		     "gff", "NULL", "NULL");
    

}

exit;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+


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


sub apollo_convert
{
#-----------------------------+
# CONVERT AMONG FILE FORMATS  |
# USING THE APOLLO PROGRAM    |
#-----------------------------+
# Converts among the various data formats that can be used 
# from the command line in tbe Apollo program. For example
# can convert GFF format files into the game XML format.
# NOTES:
#  - Currently assumes that the input file is in the correct
#    coordinate system.
#  - GFF files will require a sequence file
#  - ChadoDB format will require a db password


    # ApPath - the path of dir with the Apollo binary
    #          Specifying the path will allow for cases
    #          where the program is not in the PATHS
    # ApCmd  - the apollo commands to run

    my ($InFile,$InForm,$OutFile,$OutForm,$SeqFile,$DbPass) = @_;

    #InFile = $_[0];        # Input file path
    #InForm = $_[1];        # Output file format:
    #                         # game|gff|gb|chadoxml|backup
    #$OutFile = $_[2];       # Output file path
    #$OutForm = $_[3];       # Ouput file foramt
    #                           # chadoDB|game|chadoxml|genbank|gff|backup
    #$SeqFile = $_[4];       # The path of the sequence file
    #                           # This is only required for GFF foramt files
    #                           # When not required this can be passed as na
    #$DbPass = $_[5];        # Database password for logging on to the 
    #                           # chado database for reading or writing.
    my $ApCmd;

    # Set the base command line. More may need to be added for different
    # formats. For example, GFF in requires a sequence file and the CHADO
    # format will require a database password.
    $ApCmd = $ap_path." -i ".$InForm." -f ".$InFile.
	" -o ".$OutForm." -w ".$OutFile;

    # Make sure that that input output formats are in lowercase
    # may need to add something here to avoid converting chadoDB
    $InForm = lc($InForm);
    $OutForm = lc($OutForm);
    
    # Determine the proper command to use based on the input format
    # since GFF file also require a sequence file
    if ($InForm =~ "gff" )
    {
	$ApCmd = $ApCmd." -s ".$SeqFile;
    }
    
    if ($InForm =~ "chadodb")
    {
	$ApCmd = $ApCmd." -D ".$DbPass;
    }

    # Do the apollo command
    system ( $ApCmd );

}


=head1 NAME

batch_game2gff.pl - Convert game.xml annotations to gff format

=head1 VERSION

This documentation refers to program version $Rev: 948 $

=head1 SYNOPSIS

=head2 Usage

    batch_game2gff.pl -i InDir -o OutDir

=head2 Required Variables

    --indir         # Path to the input directory
    --outdir        # Path to the output directory

=head1 DESCRIPTION

This is what the program does

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path to the directory that contains the game.xml files to be converted
from game.xml format to gff format.

=item -o,--outfile

Path to the output directory where the GFF format files will be saved to.

=back

=head1 OPTIONS

=over 2

=item --ap-path

The path to the local installation of the Apollo Genome Annotation program.
This can also be defined with the GP_APOLLO_BIN environment variable as 
discussed in the CONFIGURATION AND ENVIRONMENT section of the 
batch_game2gff.pl program manual.

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

=item --verbose

Run the program in verbose mode. This produces the maximum amount
of program status information and can be useful for diagnosing problems.

=back

=head1 DIAGNOSTICS

Error messages generated by this program and possible solutions are listed
below.

=over 2

=item ERROR: No game.xml files were found in the input directory

The input directory does not appear to contain game.xml files.
This could happen because you gave an incorrect path or because your annotation
files do not have the expected *.game.xml extension in the file name.

=item ERROR: Could not create the output directory

The output directory could not be created at the path you specified. 
This could be do to the fact that the directory that you are trying
to place your base directory in does not exist, or because you do not
have write permission to the directory you want to place your file in.

=back

=head1 CONFIGURATION AND ENVIRONMENT

This program does not require external configuration files but it can make 
use of variables defined in the user's environment:

=head2 User Environment

=over 2

=item DP_APOLLO_BIN

This is the path to the Apollo binary file. This can also be defined at
the command line by the --ap-path variable.

=back

An example of setting the Apollo binary variable in the bash shell is:

    export DP_APOLLO_BIN='$HOME/apps/Apollo_1.6.5/apollo/bin/apollo'

=head1 DEPENDENCIES

=head2 Required Software

=over

=item * Apollo Genome Annotation Curation Tool

This program relies on the Apollo Genome Annotation Curation tool to
convert the game.xml format files to the GFF format. This program can
be downloaded at:
http://www.fruitfly.org/annot/apollo/

=back

=head2 Required Perl Modules

=over

=item * File::Copy

This module is required to copy the BLAST results.

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

=item * Apollo GUI Required

Despite using the command line syntax for Apollo to conver the game files
to xml, this program requires that you have GUI access to the Apollo program.

=back

=head1 SEE ALSO

The batch_game2gff.pl program is part of the DAWG-PAWS package of genome
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

STARTED: 11/01/2007

UPDATED: 12/17/2007

VERSION: $Rev: 948 $

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 11/01/2007
# - Main body of program written
#
# 12/17/2007
# - Added SVN tracking of Rev
# - Updated POD documentation
# - Added print_help subfunction that extracts help message
#   from the POD documentation
