#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# gff_seg.pl - Segmentation of a large gff file             |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 02/17/2008                                       |
# UPDATED: 04/01/2010                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Given a gff file that contains point or segement data    |
#  will extract segments that exceed a threshold value or   |
#  array of threshold values. Will create a gff segment     |
#  file as well as a gff parse file. The segment file       |
#  converts the vals to segments that meet the threshold    |
#  criteria while the parse file returns all points or      |
#  segments in the input file that exceed the threshold     |
#  value.                                                   |
#                                                           |
# USAGE:                                                    |
#  gff_seg.pl -i infile.gff -s outfile.gff - p parseout.gff |
#             -t [integer] --min [int] --max [int]          |
#                                                           |
# VERSION: $Rev: 874 $                                      |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+
# TO DO -- accept t as array of values 
#          ie -t 1,10,100,1000,1000
# with the bins
# for integers this is 
# 0 ot 1   0<x<2
# 2 to 10   1<x<11
# 11 to 100   10<x<11
# will run through set
# see if cur_t_bin is different then previous t bin
# and draw output according.
# this would also be a good place to 

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
my ($VERSION) = q$Rev: 874 $ =~ /(\d+)/;
# Get GFF version from environment, GFF2 is DEFAULT
my $gff_ver = uc($ENV{DP_GFF}) || "GFF2";

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
# Vars that take values from command line
my $infile;                   # Gff infile
my $thresh;                   # Threshold value for parsing, segmenting
my $outfile_seg;              # Path to the segmentation output file
my $outfile_parse;            # Path to the parsed output file
#my $outdir;                  # Base outdir for storing output when $thresh is an array
my $min_len;                  # Minimum segment length to report
#my $max_len;                  # Maximum segment length to report

# Booleans
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $source;                  # The program name
my $param;                    # The parameter tag
my $seqid;                    # The sequence ID

# THRESHOLD MATRIX
my @tmatrix;                  # Threshold value matrix
my $col;                      # index col in tmatrix
my $col_max;                  # (Num of threshold vals)
                              # index starts at one (0 is labels)
my $row;                      # index base in tmatrix
my $row_max;                  # (Num of bases)                  
                              # index starts at one (0 are labels)
my $dep;                      # The depth of sampling of individual base
my $dep_max;                  # (length of feature segments)
my $i;                        # Temp matrix index
my $j;                        # Temp matrix index

# Summary of tmatrix
my @tsummary;


# GFF IN VALS
my $in_seq_name;       # 0 - Sequence name
my $in_source;         # 1 - Source of the feature
my $in_feat;           # 2 - feature
my $in_start;          # 3 - Start of the feature
my $in_end;            # 4 - End of the feature
my $in_score;          # 5 - Score, This is oligo count too
my $in_strand;         # 6 - Strand of the feature
my $in_frame;          # 7 - Frame of the feature
my $in_attribute;      # 8 - Attribute info

my $base_attribute;    # Base name for new segment attribute
my $attribute;         # Full attribute name

# Feature name to use in output
# ie mathmatetically_defined_repeat
my $feature;

# Allow threshold as array value:
#     --thresh | -t for threshold values
# Use --lables | -l for score labels in output
# then test to see if -t is equal to l in length
# default
#my @thresh_vals = (0,1,10,100,1000,10000,100000);
my @thresh_vals;
#my @thresh_vals = (1,10,100,1000,10000,100000);
my @label_vals;

# Optional output files
my $sgr_out;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED VARIABLES
		    "i|infile=s"       => \$infile,
                    "s|seg-out=s"      => \$outfile_seg,
		    "p|parse-out=s"    => \$outfile_parse,
#		    "t|thresh=s"       => \$thresh,
		    "t|thresh=s@"    => \@thresh_vals,
		    # OPTIONS
		    "sgr=s"            => \$sgr_out,
		    "gff-ver=s"        => \$gff_ver,
		    "attribute=s"      => \$base_attribute,
		    "feature=s"        => \$feature,
		    "seqname=s"        => \$seqid,
		    "program|source=s" => \$source,
		    "param=s"          => \$param,
		    "min-len=s"        => \$min_len,
		    "q|quiet"          => \$quiet,
		    "verbose"          => \$verbose,
		    # ADDITIONAL INFORMATION
		    "usage"            => \$show_usage,
		    "version"          => \$show_version,
		    "man"              => \$show_man,
		    "h|help"           => \$show_help,);

if (@thresh_vals) {
    @thresh_vals = split(/,/,join(',',@thresh_vals));
}
else {
    @thresh_vals = (1,10,100,1000,10000,100000);
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
    print "\nbatch_mask.pl:\n".
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
#if  ( !$thresh ) {
#    print STDERR "\a";
#    print STDERR "ERROR: A threshold value must be specified at the".
#	" command line\n" if (!$thresh);
#    exit;
#}

#-----------------------------------------------------------+ 
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+ 

#-----------------------------+
# OPEN FILE HANDLES           |
#-----------------------------+
if ( ($outfile_seg) ) {
    open (SEGOUT, ">$outfile_seg")
	|| die "ERROR: Can not open segment out file:\n$outfile_seg\n";
    if ($gff_ver =~ "GFF3") {
	print SEGOUT "##gff-version 3\n";
    }
}

if ( ($outfile_parse) ) {
    open (PAROUT, ">$outfile_parse")
	|| die "ERROR: Can not open parse out file:\n$outfile_parse\n";
    if ($gff_ver =~ "GFF3") {
	print PAROUT "##gff-version 3\n";
    }
}


# If neither a seg file or par file are given then do seg file
# to STDOUT
if ( (!$outfile_seg) && (!$outfile_parse) ) {
    open (SEGOUT, ">&STDOUT")
	|| die "ERROR: Can not print to STDOUT\n";
    if ($gff_ver =~ "GFF3") {
	print SEGOUT "##gff-version 3\n";
    }
}



if ($infile) {
    open (GFFIN, "<$infile")
	|| die "ERROR: Can not open gff input file:\n$infile\n";
} else {
    print STDERR "Expecting input from STDIN\n";
    open (GFFIN, "<&STDIN") ||
	die "Can not accept input from standard input.\n";
}

# SORT THE THRESHOLD ARRAY

#-----------------------------+
# PROCESS GFF FILE            |
#-----------------------------+

#my $do_print_seg = 0;  # BOOLEAN
#
my $line_count = 0;

# LOAD HEADERS TO TMATRIX

# LOAD DATA TO TMATRIX
if ($verbose) {
    print STDERR "\n----------------------------------\n";
    print STDERR " Loading Threshold Matrix ... \n";
    print STDERR "----------------------------------\n";
}

my $feat_len;    # Length of feature
my $feat_start;  # Start of the feature
my $feat_end;    # End of the feature
my $feat_bases;  # Numbe of bases in feature
my $min_seq;     # The max sequence position
my $max_seq;     # The max sequence position
my $max_len;     # The max sequence length

while (<GFFIN>) {

    # TO do add ignore of comment lines
    $line_count++;

    chomp;
    my @gff_data = split(/\t/, $_);
    my $num_cols = @gff_data;

    # Get sequence id
    unless ($seqid) {
	$seqid = $gff_data[0];
    };

    # Get source
    unless ($source) {
	$source = $gff_data[1];
    }
    
    # Get feature
    unless ($feature) {
	$feature = $gff_data[2];
    }

    # Get base attribute
    unless ($base_attribute) {
	$base_attribute = $gff_data[8];
    }

    my $cur_bin_val = &get_bin_val($gff_data[5], @thresh_vals);
    $line_count++;
    if ($num_cols == 9) {    
	@gff_data = split(/\t/, $_);
    }
    else {
	print STDERR "WARNING: Unexpected number of cols in input line".
	    "$line_count";
	next;
    }

    # Get start and end values, flip names if end is before start
    $feat_start = int($gff_data[3]);
    $feat_end = int($gff_data[4]);
    if ($feat_end < $feat_start) {
	$feat_start = int($gff_data[4]);
	$feat_end = int($gff_data[3]);
    }
    $feat_len = $feat_end - $feat_start;

    # Get binned score value
    my $bin_val = &get_bin_val($gff_data[5], @thresh_vals);

    #-----------------------------+
    # LOAD BIN VALS TO tmatrix    |
    #-----------------------------+
    # This will do terrible for with an index that does not
    # start at zero
#    print STDERR $feat_start."->".$feat_end."\n";
    my $i_end = $feat_end + 1;
#    print STDERR "\t";
    
    # if we can assume all features of same length
    $j = $feat_len + 1;


    # Get sequence bounds
    if ($min_seq) {
	if ($min_seq > $feat_start) {
	    $min_seq = $feat_start;
	}
    }
    else {
	$min_seq = $feat_start;
    }

    if ($max_seq) {
	if ($max_seq < $feat_end) {
	    $max_seq = $feat_end;
	}
    }
    else {
	$max_seq = $feat_end;
    }

    # Get max feature length
    if ($max_len) {
	if ($max_len < $feat_len) {
	    $max_len = $feat_len;
	}
    }
    else {
	$max_len = $feat_len;
    }

    for ($i = $feat_start; $i < $i_end; $i++) {
	#print STDERR "$i:$j|$bin_val-";
	$j--;
	unless ($bin_val =~ "NULL") {
	    $tmatrix[$i][$j] = $bin_val;
	}
    }
#    print STDERR "\n";

}


#-----------------------------+
# LOAD SUMMARY ARRAY          |
#-----------------------------+
# Currently will just get the max value
print STDERR "\nSeq bounds: $min_seq - $max_seq \n\n" if $verbose;
my $seq_end = $max_seq + 1;
my $col_end = $max_len + 1;
for ($i = 1; $i < $seq_end; $i++) {
    for ($j = 1; $j < $col_end; $j++) {
	if ($tmatrix[$i][$j]) {
	    if ($tsummary[$i]) {
		# The following will just give the max value
		if ( $tmatrix[$i][$j] > $tsummary[$i]) {
		    $tsummary[$i] = $tmatrix[$i][$j];
		}
	    }
	    else {
		$tsummary[$i] = $tmatrix[$i][$j];
	    }
	}
    }
}

#-----------------------------+
# REPORT SUMMARY ARRAY        |
#-----------------------------+
my $prev_bin_val;
# Vars that take values from the gff file
my $seg_start;
my $seg_end;
my $in_seg = 0;               # Boolean, TRUE if in SEG
my $seg_count = 1;            # Count of segment


# Output options are GFF, SGR and 
if ($sgr_out) {
    open (SGROUT, ">$sgr_out") ||
	die "ERROR: Can not open sgr file for output $sgr_out\n";
}

# I need GFF output names as
#  Vmatch:20mer_454_100
my $gff_src;

for ($i = 1; $i < $seq_end; $i++) {
    if ($tsummary[$i]) {

	# SGR OUTPUT
	if ($sgr_out) {
	    print SGROUT $seqid."\t".$i."\t".$tsummary[$i]."\n";
	}

	# GFF OUTPUT
	if ($prev_bin_val) {
	    # Temp for the moment
	    #print STDERR $prev_bin_val." vs. ".$tsummary[$i]."\n";
	    if ($prev_bin_val == $tsummary[$i]) {
		$seg_end = $i;
	    }
	    else {
		$gff_src = $source."_".$prev_bin_val;
		$attribute = $base_attribute."-".$seg_count;
		if ($gff_ver =~ "GFF3") {
		    $attribute = "ID=".$attribute;
		}
		print SEGOUT $seqid."\t".
		    $gff_src."\t".
		    $feature."\t".
		    $seg_start."\t".
		    $seg_end."\t".
		    ".\t".
		    ".\t".
		    $prev_bin_val."\t".
		    $attribute.
		    "\n";
		#print STDERR "from diff\n";
		$seg_start = $i;
		$seg_end = $i;
		$prev_bin_val = $tsummary[$i];
		$in_seg = 1;
		$seg_count++;
	    }
	}
	else {
	    # initialize on first record
	    $in_seg = 1;
	    $seg_start = $i;
	    $seg_end = $i;
	    $prev_bin_val = $tsummary[$i];
	}
    }
    else {
	if ($in_seg) {
	    $gff_src = $source."_".$prev_bin_val;
	    $attribute = $base_attribute."-".$seg_count;
	    if ($gff_ver =~ "GFF3") {
		$attribute = "ID=".$attribute;
	    }
	    print SEGOUT $seqid."\t".
		$gff_src."\t".
		$feature."\t".
		$seg_start."\t".
		$seg_end."\t".
		".\t".
		".\t".
		$prev_bin_val."\t".
		$attribute.
		"\n";
	    #print STDERR "from void\n";
	    $in_seg = 0;
	    $prev_bin_val = "";
	    $seg_count++;
	}
    }
}

# If we drop off in segmetn this print last record
if ($in_seg) {
    $gff_src = $source."_".$prev_bin_val;
    $attribute = $base_attribute."-".$seg_count;
    if ($gff_ver =~ "GFF3") {
	$attribute = "ID=".$attribute;
    }
    print SEGOUT $seqid."\t".
	$gff_src."\t".
	$feature."\t".
	$seg_start."\t".
	$seg_end."\t".
	".\t".
	".\t".
	$prev_bin_val."\t".
	$attribute.
	"\n";
    #print STDERR "Last record\n";
}

# CLOSE FILE HANDLES
close (SGROUT) if ($sgr_out);
close (SEGOUT) if ($outfile_seg);
close (PAROUT) if ($outfile_parse);

exit;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub print_seg {
    my ($seq_name,
	$source_val,
	$feature_val,
	$start_val, 
	$end_val,
	$score_val,
	$attribute_val) = @_;

#    print SEGOUT $seq_name."\t".  # seq name
    print STDERR $seq_name."\t".  # seq name
	$source_val."\t".         # program source
	$feature_val."\t".        # feature
	$start_val."\t".          # start
	$end_val."\t".            # end
	$score_val."\t".          # score (bin num)
	".\t".                    # strand
	".\t".                    # frame
	$attribute_val."\n";      # attribute
    
}

sub get_bin_val {
    
    # RETURNS THE TEXT NULL FOR VALUES less than min
    # threshold, this allows you to set 0 as a minimum threshold
    # too.
    # this requres that the bin vals array is sorted
    my ($qry_val, @bin_vals) = @_;
#    $qry_val = int($qry_val);
    my $bin_ans;
    my $i;
    $i = -1;
    my $min_val = int($bin_vals[0]);

    # Ignore values less than the minimum threshold
    if ( $qry_val < $min_val ) {
	return "NULL";
    }

    # cycle through and get bin value then return
    for my $bin_val (@bin_vals) {
#	$bin_val = int($bin_val);
	if ( $qry_val >= $bin_val) {
	    $i++;
	    print STDERR "\t".$i." - is ".$qry_val." >= ".$bin_val."\n"
		if $verbose;
	}
	else {
	    if ($verbose) {
		print STDERR "\t".$i." - is ".$qry_val." >= ".$bin_val."\n";
		print STDERR "\tReturning ".$i."\n";
		print STDERR "\tVal ".$bin_vals[$i]."\n";
	    }
	    return $bin_vals[$i];
	}

    }


    # this will return the index val in the input array, 
    # this index value
    # can be used to look up the bin value in the label array
    # Doing it this way will allow for multiple types of bins to be
    # reported as scores, and easily do simple transformations

    
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

gff_seg.pl - Segment and parse a large gff file

=head1 VERSION

This documentation refers to program version $Rev: 874 $

=head1 SYNOPSIS
    
=head2 Usage

    gff_seg.pl -i infile.gff -s seg_out.gff -p par_out.gff -t integer

=head2 Required Arguments

    -i,--infile         # Path to the input file
    -s,--seg-out        # Path to the segmented output file
    -p,--parse-out      # Path to the parsed output file
    -t,--thresh         # Threshold value

=head1 DESCRIPTION

Given a gff file that contains point or segement data
will extract segments that exceed a threshold value or
array of threshold values. Will create a gff segment
file as well as a gff parse file. The segment file
converts the vals to segments that meet the threshold
criteria while the parse file returns all points or
segments in the input file that exceed the threshold
value.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path of the input file.

=item -s,--seg-out

Gff file of the segmented data

=item -p,--parse-out

Gff file of the parsed data. This will include all rows in the original
file that exceed the threshold value.

=item -t,--thresh

Threshold value. A single integer that represents the threshold value.
The output will return all features that are greater then or equal
to the threshold value.

=back

=head1 OPTIONS

=over 2

=item --program

The program used to generate the gff result. This is the value in the 
second column of the GFF file.
By default, the program name used in the original GFF file will be used.

=item --param

The parameter used to generate the segmentation. For example, 20mer_100x
for 20mer oligos with a threshold value of 100x coverage.

=item --seqname

The name of the sequence file being annotated. This is the first column
of data in the gff output file.

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

The list of error messages that can be generated,
explanation of the problem
one or more causes
suggested remedies
list exit status associated with each error

=head1 CONFIGURATION AND ENVIRONMENT

Names and locations of config files
environmental variables
or properties that can be set.

=head1 DEPENDENCIES

Other modules or software that the program is dependent on.

=head1 BUGS AND LIMITATIONS

Any known bugs and limitations will be listed here.

=head1 REFERENCE

Please refer to the DAWGPAWS manuscript in Plant Methods when describing
your use of this program:

JC Estill and JL Bennetzen. 2009. 
"The DAWGPAWS Pipeline for the Annotation of Genes and Transposable 
Elements in Plant Genomes." Plant Methods. 5:8.

=head1 LICENSE

GNU General Public License, Version 3

L<http://www.gnu.org/licenses/gpl.html>

THIS SOFTWARE COMES AS IS, WITHOUT ANY EXPRESS OR IMPLIED
WARRANTY. USE AT YOUR OWN RISK.

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 02/17/2008

UPDATED: 04/14/2010

VERSION: $Rev: $

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 02/17/2008
# - Starting program. Treating thresh as a single value
#   later updates will need to accept an array here.
# - ASSUMES CONTIGUOUS FEATURES
# - Overlap in segments okay
# - Two output files will be created:
#     - parse_out - All features exceeding the threshold value
#     - seg_out - All segments exceeding the threshold value
# - Added the min-len variable to set the minimum length of
#   a segment to report
#
# 04/14/2010
# - Matrix based, max base report. For visualizing oligomer 
#   counts in Apollo.
