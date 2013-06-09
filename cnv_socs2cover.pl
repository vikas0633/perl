#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_socs2cover.pl - Convert socs map to oligo coverage    |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 12/09/2009                                       |
# UPDATED: 12/09/2009                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Convert socs SOLiD mapping output to a coverage map      |
#  for the input sequence with GFF output.                  |
#                                                           |
# USAGE:                                                    |
#  ShortFasta Infile.fasta Outfile.fasta                    |
#                                                           |
# VERSION: $Rev$                                            |
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
my ($VERSION) = q$Rev$ =~ /(\d+)/;

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $infile;
my $outfile;
my $qry_name;
my $feature_type;

# BOOLEANS
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_test = 0;                  # Run the program in test mode

#-----------------------------+
# SOCS VARIABLES              |
#-----------------------------+
my @socs_coverage = ();        # Coverage

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"  => \$infile,
                    "o|outfile=s" => \$outfile,
		    # ADDITIONAL OPTIONS
		    "n|s|seqname|name=s"  => \$qry_name,
		    "f|feature=s"  => \$feature_type,
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
# MAIN PROGRAM BODY           |
#-----------------------------+

# Set the names of the files that will be opened
my $socs_plus_map = $infile."\.\+\.map";
my $socs_neg_map = $infile."\.\-\.map";
my $socs_plus_amb = $infile."\.\+\.amb";
my $socs_neg_amb = $infile."\.\+\.amb";;


# Testing variables
print STDERR "Searching for files:\n";
print STDERR "\t$socs_plus_map\n";
print STDERR "\t$socs_neg_map\n";
print STDERR "\t$socs_plus_amb\n";
print STDERR "\t$socs_neg_amb\n";

#-----------------------------------------------------------+
# PROCESS SOCS FILES                                        |
#-----------------------------------------------------------+
# The socs files to process
# be default will be all
my @socs_files = ($socs_plus_map,
		  $socs_neg_map,
		  $socs_plus_amb,
		  $socs_neg_amb);

# Process each SOCS output file
my $socs_file = $socs_files[1]; # for testing

foreach my $socs_file (@socs_files) {
    print STDERR "Processing: $socs_file\n";
    open (SOCSMAP, "<$socs_file") ||
	die "ERROR: Can not open plus strand map file:\n$socs_file\n";
    
    # Using $l to index the file line
    my $i = 0;
    while (<SOCSMAP>) {
	chomp;
#	print STDERR $i.":".$_."\n";
	my $cov_val = int($_);
#	print STDERR "\tas int:".$cov_val."\n";
	
	if ($socs_coverage[$i]) {
	    $socs_coverage[$i] = $socs_coverage[$i] + $cov_val;
#	    print STDERR "\tcov_val:".$socs_coverage[$i]."\n";
	} else {
	    $socs_coverage[$i] = $cov_val;
#	    print STDERR "\tcov_val:".$socs_coverage[$i]."\n";
	}
	
	$i++;
    }
    
    close SOCSMAP;
    
}


#-----------------------------------------------------------+
# PROCESS OUTPUT                                            |
#-----------------------------------------------------------+
unless ($qry_name) {
    $qry_name = $infile;
}
my $source = "socs";
my $feature = "SOLiD";
my $strand = ".";

open (GFFOUT, ">$outfile") ||
    die "Can not open file for output\n";
my $i=1;
foreach my $cov_result (@socs_coverage) {
    print GFFOUT "$qry_name\t". # SeqName
	"$source\t".            # Source (Socs program)
	"$feature\t".           # Feature type name
	"$i\t".                 # Start
	"$i\t".                 # End 
	"$cov_result\t".        # Score (Coverage)
	"$strand\t".            # Strand
	".\t".                  # Frame
	"sc_$i".                # Feature Name
	"\n";
    $i++;
}

$i++;
close (GFFOUT);

exit 0;

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

1;
__END__

=head1 NAME

cnv_socs2cover.pl - Convert socs map to oligo coverage

=head1 VERSION

This documentation refers to program version 0.1

=head1 SYNOPSIS

  USAGE:
    Name.pl -i InSeq -o OutFile

    --infile        # The root name of the sequence file 
    --outfie        # Path to the output file in GFF format

=head1 DESCRIPTION

Takes as its input a set of files output from the socs program, and generates
a GFF format file that shows coverage. By default this will combine output
for unambiguous and ambiguous coverage for both strangs.

=head1 COMMAND LINE ARGUMENTS

=head2 Required Arguments

=over 2

=item -i,--infile

Path of the input file.

=item -o,--outfile

Path of the output file.

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

=head1 LICENSE

GNU General Public License, Version 3

L<http://www.gnu.org/licenses/gpl.html>

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
