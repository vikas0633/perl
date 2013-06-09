#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# Name.pl                                                   |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 00/00/2007                                       |
# UPDATED: 00/00/2007                                       |
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

# Include the existing batch_genscan
#require "/Users/jestill/code/dawgpaws/trunk/scripts/batch_genscan.pl";

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


my $seq_id = "test"; 
my $param;
my $program;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"    => \$infile,
                    "o|outfile=s"   => \$outfile,
		    "s|seqname=s"   => \$seq_id,
		    # ADDITIONAL OPTIONS
		    "param=s"       => \$param,
		    "program=s"       => \$program,
		    "gff-ver=s"     => \$gff_ver,
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


# DO the conversion
genscan2gff($infile, $outfile, $seq_id, $program, $param);

exit 0;

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
