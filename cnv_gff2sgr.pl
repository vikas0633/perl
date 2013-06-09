#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_gff2sgr.pl - Converts GFF file to SGR format          |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 12/10/2009                                       |
# UPDATED: 12/10/2009                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Converts a GFF file format to the simple graph (SGR)     |
#  format used by the IGB program and compatible with       |
#  the Apollo genome browswer. Includes the ability to      |
#  tranform the data output. Currently supported data       |
#  transformation is limited to a log10 transfrom of        |
#  integers.                                                |
#                                                           |
# USAGE:                                                    |
#  cnv_gff2sgr.pl -i infile.gff -o outfile.sgr              |
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
use POSIX qw(log10);           # For LOG10 transform

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/;

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $infile;
my $outfile;
my $transform;

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
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"    => \$infile,
                    "o|outfile=s"   => \$outfile,
		    # ADDITIONAL OPTIONS
		    "t|transform=s" => \$transform,
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
# FILE IO                     |
#-----------------------------+
if ($infile) {
    open (GFFIN, "<$infile") ||
	die "Can not open input file at:\n$infile\n";
}
else {
    open (GFFIN, "<&STDIN") ||
	die "Can not open STDIN for intput\n";
    print STDERR "Expecting input from STDIN\n";
}

if ($outfile) {
    open (SGROUT, ">$outfile") ||
	die "Can not open output file at::\t$outfile\n";
}
else {
    open (SGROUT, ">&STDOUT") ||
	die "Can not open STDOUT for output\n";
}

#-----------------------------+
#
#-----------------------------+
while (<GFFIN>) {
    chomp;
    my @gff_data = split;
    my $score = $gff_data[5];

    # Transform of score here
    if ($transform) {
	if ($transform =~ "log") {
	    $score = (int($score));
	    unless ($score == 0) {
		$score = log10($score);
	    }
	}
	elsif ($transform =~ "preston") {
	    $score = (int($score));
#	    print STDERR "send".$score."\t";
	    $score = &val2octave($score);
	}
    }

    # The following only reports the start position
    # it would be better to do midpoint if segment data
    print SGROUT $gff_data[0]."\t".
	$gff_data[3]."\t".
	$score.
	"\n";
}


close (GFFIN);
close (SGROUT);

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

sub val2octave {

# Fits a value to an octave.
# Takes a numeric value as its input.
# Ref:
# Preston, F. W. (1948). "The commonness and rarity of species." 
# Ecology 29: 254–283.
# Currently doing this as a lookup, may be better to jsut do the math
# and return an integer	

my ($inval) =  @_;
my $retval = 0;

#$inval = int($inval);
#print STDERR "IN".$inval."\n";
#if ($inval > 10) {
#    print STDERR "YUP";
#}

# Brute force if/then below..should use hashes
# for this
# this is not efficient since most values will be less than 10 or sho
if ($inval > 134217728) {
    $retval = 28;
}
elsif ($inval > 67108864) {
    $retval = 27;
}
elsif ($inval > 33554432) {
    $retval = 26;
}
elsif ($inval > 16777216) {
    $retval = 25;
}
elsif ($inval > 8388608) {
    $retval = 24;
}
elsif ($inval > 4194304) {
    $retval = 23;
}
elsif ($inval > 2097152) {
    $retval = 22;
}
elsif ($inval > 1048576) {
    $retval = 21;
}
elsif ($inval > 524288) {
    $retval = 20;
}
elsif ($inval > 262144) {
    $retval = 19;
}
elsif ($inval > 131072) {
    $retval = 18;
}
elsif ($inval > 65536) {
    $retval = 17;
}
elsif ($inval > 32768) {
    $retval = 16;
}
elsif ($inval > 16384) {
    $retval = 15;
}
elsif ($inval > 8192) {
    $retval = 14;
}
elsif ($inval > 4096) {
    $retval = 13;
}
elsif ($inval > 2048) {
    $retval = 12;
}
elsif ($inval > 1024){
    $retval = 11;
}
elsif ($inval > 512){
    $retval = 10;
}
elsif ($inval > 256) {
    $retval = 9;
}
elsif ($inval > 128){
    $retval = 8;
}
elsif ($inval > 64) {
    $retval = 7;
}
elsif ($inval > 32) {
    $retval = 6;
}
elsif ($inval > 16) {
    $retval = 5;
}
elsif ($inval > 8 ) {
    $retval = 4;
}
elsif ($inval > 4) {
    $retval = 3;
}
elsif ($inval > 2) {
    $retval = 2;
}
elsif ($inval > 1) {
    $retval = 1;
}

#print $inval.":".$retval."\n";

return ($retval);

}

1;
__END__

=head1 NAME

Name.pl - Short program description. 

=head1 VERSION

This documentation refers to program version 0.1

=head1 SYNOPSIS

  USAGE:
    Name.pl -i InFile -o OutFile

    --infile        # Path to the input file
    --outfie        # Path to the output file

=head1 DESCRIPTION

This is what the program does

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
