#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_repseek2gff.pl - Convert repseek output to gff format |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 11/24/2008                                       |
# UPDATED: 03/02/2010                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Convert repseek output to a GFF format compatible with   |
#  the Apollo genome annotation program.                    |
#                                                           |
# USAGE:                                                    |
#  cnv_repseek2gff.pl -i repseek_out.txt -o repseek.gff     |
#                                                           |
# VERSION: $Rev: 948 $                                      |
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
# Get GFF version from environment, GFF2 is DEFAULT
my $gff_ver = uc($ENV{DP_GFF}) || "GFF2";

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $infile;
my $outfile;
my $seq_id = "seq";        # Default seqname is seq
my $parameter_set = 0;     # The parameter set used, set to default

my $program = "repseek";

# BOOLEANS
my $test = 0;
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
my $ok = GetOptions(# REQUIRED ARGUMENT
		    "i|infile=s"  => \$infile,
                    "o|outfile=s" => \$outfile,
		    # OPTIONS
		    "gff-ver=s"     => \$gff_ver,
		    "s|seqname=s" => \$seq_id,
		    "program=s"   => \$program,
		    "p|param=s"   => \$parameter_set,
		    "q|quiet"     => \$quiet,
		    "verbose"     => \$verbose,
		    # ADDITIONAL INFORMATION
		    "usage"       => \$show_usage,
		    "test"        => \$test,
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
    print "\ncnv_repseek2gff.pl:\n".
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
# MAIN PROGRAM BODY           |
#-----------------------------+
# Just use the repseek2gff subfunction to do the conversion

&repseek2gff ($seq_id, $parameter_set, $infile, $outfile, $program);

exit 0;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub repseek2gff {

    # Need to convert this to using hashes and use STDIN,STOUT
    # when file path arguments are not passed
    
    # $seq_id is the id of the query sequence
    # $repin is the path to the file to convert to gff
    # $repout is the path to the gff output file

    my ($seqname, $param_set, $repin, $repout, $source) = @_;

    #-----------------------------+
    # OPEN FILE HANDLES           |
    #-----------------------------+
    # DEFAULT TO STDIN
    if ($repin) {
	open (REPIN, "<$repin") ||
	    die "Can not open input file:\n$repin\n";
    } 
    else {
	print STDERR "Expecting input from STDIN\n";
	open (REPIN, "<&STDIN") ||
	    die "Can not accept input from standard input.\n";
    }

    # DEFAULT TO STDOUT
    if ($repout) {
	open (GFFOUT, ">$repout") ||
	    die "Can not open output file:\n$repout\n";
    }
    else {
	open (GFFOUT, ">&STDOUT") ||
	    die "Can not print to STDOUT\n";
    }

    my $palindrome_id = 0; # Running palindrome id number
    my $tandem_id = 0;
    my $close_id = 0;
    my $overlap_id = 0;
    my $interseq_id = 0;
    my $repseek_id = 0;

    my $inv_id = 0;        # Counter/id for inverted repeats
    my $dir_id = 0;        # Counter/id for direct repeats

#    my $source;            # The data for the source column
#    my $source = "repseek";

    if ($param_set) {
	$source = $source.":".$param_set;
    }


    # MAY CONSIDER JUST RETURNING THE TANDEM REPEATS ?
    while (<REPIN>) {
	chomp;

	$repseek_id = $repseek_id + 1;
	my @rep_parts = split;
	my @type_parts = split(/\./, $rep_parts[0] );

	my $repeat_type = $type_parts[0];
	my $repeat_direction = $type_parts[1];

	my $copy1_start = $rep_parts[1];
	my $copy2_start = $rep_parts[2];
	my $copy1_len = $rep_parts[3];
	my $copy2_len = $rep_parts[4];
	my $copy1_end = $copy1_start + $rep_parts[3];
	my $copy2_end = $copy2_start + $rep_parts[4];
	my $copy_distance = $rep_parts[5];            # distance between the repeats
	my $percent_identity = $rep_parts[7];   #  (matches / alignment_length)\n
	my $alignment_score = $rep_parts[8];

	# Will split the repseek source into inverted repeats 
	# vs direct repeats .. putative TEs

	# Source sets
	# - Overlap
        # - Palindromes
	my $feature;

	# Translate repeat direction
	if ($repeat_direction =~ "inv") {
	    $repeat_direction = "inverted";
	}
	elsif ($repeat_direction =~ "dir") {
	    $repeat_direction = "direct";
	}

	if ($repeat_type =~ "Overlap") {
	    #$feature = "overlapping_"."$repeat_direction"."_repeat";
	    # overlapping repeats not recognized in SO
	    $feature = $repeat_direction."_repeat";
	}
	else {
	    $feature = "$repeat_direction"."_repeat";
	}

	if ($repeat_type =~ "Palindrome") {
	    # palindromic repeats are currently not recognized
	    # in the sequence ontology, will have to tag this
	    # feature with additional information in the
	    # attribute field.
	    #$feature = "palindromic_repeat";
	    $feature = "inverted_repeat";
	}

	my $attribute = "repseek".$repseek_id."_".$repeat_direction;

	#-----------------------------+
	# PRINT OUTPUT TO GFF         |
	#-----------------------------+
	my $parent_id;
	if ($gff_ver =~ "GFF3") {

	    #-----------------------------+
	    # PRINT PARENT                |
	    #-----------------------------+
	    $parent_id = $attribute;
	    my $parent_attribute = "ID=".$parent_id."_repeat";
	    
	    if ($repeat_type =~ "Palindrome") {
		$parent_attribute = $parent_attribute.
		    " ; Alias=palindrome";
	    }
	    elsif ($repeat_type =~ "Overlap") {
		$parent_attribute = $parent_attribute.
		    "; Alias=overlapping_direct_repeat";
	    }
	    elsif ($repeat_type =~ "Tandem") {
		$parent_attribute = $parent_attribute.
		    "; Alias=tandem_direct_repeat";
	    }
	    # Add not about other attribute
	    $parent_attribute = $parent_attribute.
		"; Note= ".$percent_identity.
		" ".$copy1_len.
		" ".$copy2_len.
		" ".$copy_distance ;


	    # GET SPAN LOCATIONS          |
	    my @locations = ( int($copy1_start), 
			      int($copy1_end), 
			      int($copy2_start), 
			      int($copy2_end) );
	    @locations = sort @locations;
	    my $span_start = $locations[0];
	    my $span_end = $locations[3];

	    # PRINT 
	    print GFFOUT "$seqname\t".       # Seqname
		"$source\t".                 # Source
		"$feature\t".                # Feature
		"$span_start\t".             # Start
		"$span_end\t".               # End
		"$alignment_score\t".        # Score
		"+\t".                       # Strand
		".\t".                       # Frame
		"$parent_attribute".         # Attribute
		"\n";

	}
	else {
	    $attribute = $attribute
	}


	if ($gff_ver =~ "GFF3") {
	    $attribute = "ID=".$parent_id."_copy1".
		"; Parent=".$parent_id."_repeat";
	    $feature = "repeat_unit";
	}
	print GFFOUT "$seqname\t".       # Seqname
	    "$source\t".                 # Source
	    "$feature\t".                # Feature
	    "$copy1_start\t".            # Start
	    "$copy1_end\t".              # End
	    "$alignment_score\t".        # Score
	    "+\t".                       # Strand
	    ".\t".                       # Frame
	    "$attribute".                # Attribute
	    "\n";

	if ($gff_ver =~ "GFF3") {
	    $attribute = "ID=".$parent_id."_copy2".
		"; Parent=".$parent_id."_repeat";
	    $feature = "repeat_unit";
	}
	print GFFOUT "$seqname\t".       # Seqname
	    "$source\t".                 # Source
	    "$feature\t".                # Feature
	    "$copy2_start\t".            # Start
	    "$copy2_end\t".              # End
	    "$alignment_score\t".        # Score
	    "+\t".                       # Strand
	    ".\t".                       # Frame
	    "$attribute".                # Attribute
	    "\n";

    }

    close (REPIN);
    close (GFFOUT);


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

cnv_repseek2gff.pl - Convert repseek output to gff format

=head1 VERSION

This documentation refers to program version $Rev: 948 $

=head1 SYNOPSIS

=head2 Usage

    cnv_repseek2gff.pl -i InFile -o OutFile

=head2 Required Options

    --infile        # Path to the repseek output file
                    # Assumes STDIN if not givien
    --outfile       # Path to the output gff file
                    # Assumes STDOUT if not given

=head1 DESCRIPTION

Convert repseek output to a GFF format compatible with
the Apollo genome annotation program.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path of the input file. If an input file is not provided, the program
will expect input from STDIN.

=item -o,--outfile

Path of the output file. If an output path is not provided, the program
will write output to STDOUT.

=back

=head1 OPTIONS

=over 2

=item --gff-ver

The GFF version for the output. This will accept either gff2 or gff3 as the
options. By default the GFF version will be GFF2 unless specified otherwise.
The default GFF version for output can also be set in the user environment
with the DP_GFF option. The command line option will always override the option
defined in the user environment.

=item -s,--seqname

The name of the sequence contig that is being annotated. This will be used
for the first column in the gff file. If this option is not specified
the name will default to 'seq'.

=item -p, --param

The name of the parameter set used in repseek. This allows the user to
define multiple parameter sets in repseek, and then draw them as
separate tracks in annotation visualization programs.

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

=item -v, --verbose

Run the program in verbose mode.

=back

=head1 DIAGNOSTICS

=over

=item Can not open input file

The input file path provided by the -i, --infile switch is not valid.

=item Can not open output file

The output file path provided by the -o, --outfile witch is not valid.

=item Expecting input from STDIN

When an input file path is not specified, the program expects input 
to come from STDIN. This default behaviour allows output from repseek
to be piped directly into the cnv_repseek2gff.pl program.

=back

=head1 CONFIGURATION AND ENVIRONMENT

This program does not make use of configuration files or variables
set in the user environment.

=head1 DEPENDENCIES

=head2 Required Software

=over

=item RepSeek

This program is designed to parse output from the RepSeek program. RepSeek
is available from:
http://wwwabi.snv.jussieu.fr/~public/RepSeek/

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

=item * Limited RepSeek testing

This program has been tested and know to work with RepSeek version 10May2007.
Other versions have not been tested for compatibility.

=back

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

STARTED: 11/25/2008

UPDATED: 03/02/2010

VERSION: $Rev: 948 $

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 11/24/2008
# - The earlier version of this program was lost so
#   starting this program de novo.
# - This program accepts input from STDIN and will output 
#   to STDOUT.
#
# 01/20/2009
# - Added svn revision tracking
# 
# 03/30/2009
# - changed source to not include repeat direction
# - changed feature name from exon to sequence ontology
#   complient names for repeats
#      -inv -> inverted_repeat
#      -dir -> direct_repeat
#      -overlapping_direct_repeat   
#         -- not SeqOntology complient
#      -overlapping_inverted_repeat
#         -- not SeqOntology complient
#      -palindromic_repeat
#  - Added support for --param
#  - Added support for --program
# 03/02/2010
#  - Added support for gff3 output
#
# REPSEEK OUTPUT
# Overlap.dir     117     141     48      52      -24     119-143-20-3.05 84.615  33.06   2.48    2       0.52
#Distant.dir     3980    10979   1745    1731    5254    5141-12128-20-2.78      96.167  1593.04 7.03    7       0.67
#Distant.dir     3980    13884   1749    1750    8155    3980-13884-57-5.88      86.969  1254.62 7.14    7       0.57
#Distant.dir     4592    72427   459     471     67376   4764-72599-20-4.47      63.559  241.67  6.82    8       0.36
#Distant.dir     10957   13867   1767    1776    1143    10985-13890-46-6.07     85.271  1206.85 6.92    7       0.65
#Distant.dir     14494   72427   387     388     57546   14666-72599-20-4.47     58.763  187.50  6.90    8       0.38
#Close.dir       18272   18347   35      35      40      18274-18349-33-2.00     97.143  32.56   2.00    2       1.00
#Distant.dir     28066   41668   149     149     13453   28099-41701-56-2.99     97.987  142.99  2.51    2       0.51
#Overlap.dir     29924   29954   99      99      -69     29949-29979-23-2.33     87.879  77.74   2.70    3       0.70
#Distant.dir     32892   39731   1747    1749    5092    33404-40245-96-2.07     97.086  1633.51 6.60    7       0.53
#Distant.dir     32892   79038   364     363     45782   32905-79051-30-4.20     84.110  241.70  6.96    7       0.98
#Distant.dir     32894   71717   825     810     37998   32906-71729-70-5.90     70.192  423.29  6.79    7       0.54
#Distant.dir     33585   72890   6510    6511    32795   38931-78243-59-2.51     92.172  5330.09 3.97    3       0.57
#                       "\tRepeats are displayed in 12 COLUMNS\n"
#                       "\t 01 - type of the repeat (Tandem|Close|Distant|Overlap|Palindrome|Interseq).(dir|inv)\n"
#                       "\t        'Tandem'     : repeat is a perfect tandem (spacer=0)\n"
#                       "\t        'Close'      : spacer is smaller than 1 kb\n"
#                       "\t        'Distant'    : spacer is larger  than 1 kb\n"
#                       "\t        'Overlap'    : repeat is composed of two overlapping copies (spacer<0)\n"
#                       "\t        'Palindrome' : repeat is a perfect palindrome (spacer=0)\n"
#                       "\t        'Interseq'   : repeat is shared by two sequences\n"
#                       "\t 02 - position of the first copy\n"
#                       "\t 03 - position of the second copy\n"
#                       "\t 04 - length of the first copy\n"
#                       "\t 05 - length of the second copy\n"
#                       "\t 06 - spacer of the repeats (corrected when the sequence is circular)\n"
#                       "\t 07 - characteristics of the seed that gave that repeat (pos1-pos2-len-meanR)\n"
#                       "\t 08 - %%identity between the two copies (matches / alignment_length)\n"
#                       "\t 09 - score of the alignement\n"
#                       "\t 10 - mean-R of the repeat\n"
#                       "\t 11 - mode-R of the repeat\n"
#                       "\t 12 - fraction of its mode-R\n"
#
#                                );
#
