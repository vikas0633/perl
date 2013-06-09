#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_eugene2gff.pl - Convert eugene output to GFF format   |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 03/16/2010                                       |
# UPDATED: 03/17/2010                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Short Program Description                                |
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
# Get GFF version from environment, GFF2 is DEFAULT
my $gff_ver = uc($ENV{DP_GFF}) || "GFF2";

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $infile;
my $outfile;

# BOOLEANS
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_test = 0;                  # Run the program in test mode

my $param;                    # Suffix appended to the end of the gff 
my $seqname;                  #
my $program = "Eugene";       # The source program

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"  => \$infile,
                    "o|outfile=s" => \$outfile,
		    # ADDITIONAL OPTIONS
		    "gff-ver=s"   => \$gff_ver,
		    "p|param=s"   => \$param,
		    "program=s"   => \$program,
		    "s|seqname=s" => \$seqname,
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

# CONVERT EUGENE OUTPUT TO GFF FILE
cnv_eugene2gff ($seqname, $program, $infile, $outfile);

exit 0;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+
sub cnv_eugene2gff {

    # CONVERTS THE EUGENE FORMAT GFF FILE TO THE 
    # DAWGPAWS COMPATIBLE GFF FILE
    my ($seq_name, $gff_source, $eug_gff, $gffout) =  @_;
    my @eugene_results;
    my $start;
    my $end;

    #-----------------------------+
    # EUGENE INPUT                |
    #-----------------------------+
    if ($eug_gff) {
	open (EUGIN, "<$eug_gff") ||
	    die "Can not open EuGene gff file $eug_gff";
    }
    else {
	print STDERR "Expecting input from STDIN\n";
	open (EUGIN, "<&STDIN") || 
	    die "Can not STDIN for input\n";
    }
    
    #-----------------------------+
    # GFF OUTPUT                  |
    #-----------------------------+
    if ($gffout) {
	open (GFFOUT, ">$gffout") ||
	    die "Can not open gff output\n";
    }    
    else {
	open (GFFOUT, ">&STDOUT") ||
	    die "Can not print to STDOUT\n";
    }

    # PRINT GFF3 Header
    if ($gff_ver =~ "GFF3") {
	print GFFOUT "##gff-version 3\n";
    }

    my $cur_gene_name;
    my $prev_gene_name = "null";
    my $i = -1;                     # i indexes gene count
    my $j = -1;                     # j indexes exon count

    while (<EUGIN>) {

	my @gff_parts = split;
	my $num_gff_parts = @gff_parts;
	my $exon_count = 0;
	
	# The middle part of this is gene model number
	# Name as HEX20.9.0
	# Prefix is first five characters of the fasta name
	# Middle is gene model number
	# Last part is 0 for UTRS, exon order number otherwise
	my @name_parts = split ( /\./, $gff_parts[0] );
	my $model_num = $name_parts[1];

	my $cur_gene_name = $model_num;
	
	#-----------------------------+
	# MAKE START < END            |
	#-----------------------------+
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
	    
	    $eugene_results[$i]{gene_strand} = $gff_parts[6];
	    $eugene_results[$i]{gene_name} = $gff_source."_".
		"gene_".
		$cur_gene_name;
	    $eugene_results[$i]{gene_start} = $start;
	    
	}
	
	$prev_gene_name = $cur_gene_name;

	# GET GENE END
	$eugene_results[$i]{gene_end} = $end;

	# LOAD EXON INFORMATION TO ARRAY
	$eugene_results[$i]{exon}[$j]{exon_id} = $name_parts[2];
	$eugene_results[$i]{exon}[$j]{start} = $start;
	$eugene_results[$i]{exon}[$j]{end} = $end;
	$eugene_results[$i]{exon}[$j]{score} = $gff_parts[5];
	$eugene_results[$i]{exon}[$j]{strand} = $gff_parts[6];
	$eugene_results[$i]{exon}[$j]{frame} = $gff_parts[7];
	$eugene_results[$i]{exon}[$j]{exon_type} = $gff_parts[2];
	# convert exon type to SONG compatible format

	print STDERR "MODEL NUMBER: $model_num\n" if $verbose;

	# Positions are listed as
	# UTR5
	# UTR3
	# E.Init
	# E.Term
	# E.Intr
	# E.Sngl
	# How to translate this to Apollo
	# For now ignore 5 and 3 UTR label all as exon
	# otherwise add 5 and three to either end
	unless ($gff_parts[2] =~ "UTR5" ||
		$gff_parts[2] =~ "UTR3") {
	    

	    
#	    print GFFOUT "eugene_simp_".$model_num;

	    #-----------------------------+
	    # PRINT GFFOUTPUT AS PARSE    |
	    #-----------------------------+
	    my $attribute = "eugene_simp_".$model_num;;
	    
	    # Need to replace the 
	    # Addiing simp to eugene model name
	    # to indicate that external information
	    # was not used besides the rice matrix
	    print GFFOUT "$seq_name\t";
	    print GFFOUT "eugene\t";
	    print GFFOUT "exon\t";
	    print GFFOUT $start."\t";
	    print GFFOUT $end."\t";
	    print GFFOUT $gff_parts[5]."\t";
	    print GFFOUT $gff_parts[6]."\t";
	    print GFFOUT $gff_parts[7]."\t";
	    print GFFOUT $attribute."\t";
	    print GFFOUT "\n";
	    
	    #print "\t$gff_parts[2]\n" if $verbose;
	}

    }
    

    #-----------------------------+
    # PRINT GFFOUT FROM ARRAY     |
    #-----------------------------+
    my $parent_id;
    for my $href ( @eugene_results ) {
	
	# If GFF3 need to print the parent gene span
	if ($gff_ver =~ "GFF3") {
	    $parent_id = $href->{gene_name};
	    
	    print GFFOUT $seq_name."\t".                # seq id
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
	# EXONS                       |
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
	    
	    # Currently not reporting UTRs
	    # May want to exclude UTRs from gene span reported above
	    unless ($ex->{exon_type} =~ "UTR5" ||
		    $ex->{exon_type} =~ "UTR3") {
		
		print GFFOUT $seq_name."\t".
		    $gff_source."\t".
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

    }


    close (EUGIN);
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
