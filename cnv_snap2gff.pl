#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_snap2gff.pl - Convert SNAP GFF output to GFF3 format  |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 03/29/2010                                       |
# UPDATED: 03/29/2010                                       |
#                                                           |
# DESCRIPTION:                                              |
#  This is a very limited conversion of SNAP GFF output to  |
#  GFF3 format output.                                      |
#                                                           |
# USAGE:                                                    |
#  cnv_snap2gff.pl -i infile_snap.gff -o outfile.gff        |
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
# Get GFF version from environment, GFF2 is DEFAULT
my $gff_ver = uc($ENV{DP_GFF}) || "GFF2";

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $infile;
my $outfile;
my $prev_gene_name;

# BOOLEANS
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_test = 0;                  # Run the program in test mode

my $program = "SNAP";
my $param;
my $seqname;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED VARIABLES
		    "i|infile=s"  => \$infile,
                    "o|outfile=s" => \$outfile,
		    # OPTIONS
		    "s|seqname=s" => \$seqname,
		    "p|param=s"   => \$param,
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


#-----------------------------+
# MAIN PROGRAM BODY           |
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


# DO CONVERSION
snap2gff ($program, $infile, $outfile, $seqname, $param, 0);

exit 0;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub snap2gff {
    
    my ($source, $snap_in, $gffout, $seq_id, $src_suffix, $do_append ) = @_;
    # Array to hold all snap results for a single contig
    my @snap_results;
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
		"ID=".$parent_id."\t".      # attribute
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

=head2 Usage

    Name.pl -i InFile -o OutFile

=head2 Required Arguments

    --infile        # Path to the input file
    --outfie        # Path to the output file

=head1 DESCRIPTION

This is what the program does

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path of the input file.

=item -o,--outfile

Path of the output file.

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

=item -q,--quiet

Run the program with minimal output.

=back

=head1 EXAMPLES

The following are examples of how to use this script

=head2 Typical Use

This is a typcial use case.

=head1 DIAGNOSTICS

=over 2

=item * Expecting input from STDIN

If you see this message, it may indicate that you did not properly specify
the input sequence with -i or --infile flag. 

=back

=head1 CONFIGURATION AND ENVIRONMENT

Names and locations of config files
environmental variables
or properties that can be set.

=head1 DEPENDENCIES

Other modules or software that the program is dependent on.

=head1 BUGS AND LIMITATIONS

Any known bugs and limitations will be listed here.

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
