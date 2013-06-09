#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_findmite2gff.pl - Convert findmite results to gff     |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 08/31/2007                                       |
# UPDATED: 02/26/2010                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Parses findmite results to a gff file. Assigns unique    |
#  name to all of the results                               |
#                                                           |
# VERSION: $Rev: 948 $                                            |
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
my $infile;                    # Findmite outupt file to convert
my $outfile;                   # The gff output file produced
my $fasta_outfile;             # FASTA output file of putative mites
my $seqname;                   # Name of the sequence being annotated
my $param;                     # Name of the parameter set used
my $prog = "findmite";         # Name of program

# BOOLEANS
my $help = 0;
my $quiet = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $verbose = 0;
my $gff_append = 0;
my $fasta_append = 0;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED ARGUMENTS
		    "i|infile=s"    => \$infile,
                    "o|outfile=s"   => \$outfile,
		    # OPTIONS
		    "gff-ver=s"     => \$gff_ver,
		    "n|name=s"      => \$seqname,
		    "p|param=s"     => \$param,
		    "program=s"     => \$prog,
		    "f|fasta=s"     => \$fasta_outfile,
		    "q|quiet"       => \$quiet,
		    "verbose"       => \$verbose,
		    "append-gff"    => \$gff_append,
		    "append-fasta"  => \$fasta_append,
		    # Additional information
		    "usage"         => \$show_usage,
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

#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+
print STDERR "$0 has started\n" if $verbose;

if ($fasta_outfile) {
    if ($fasta_append) {
	findmite2gff ($prog, $infile, $outfile, $param, $gff_append, 
		      $seqname, $fasta_outfile, $fasta_append);
    } 
    else {
	findmite2gff ($prog, $infile, $outfile, $param, $gff_append, 
		      $seqname, $fasta_outfile, $fasta_append);	
    }
}
else {
    findmite2gff ($prog, $infile, $outfile, $param, $gff_append,
		  $seqname, 0, 0);
}

exit;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub findmite2gff {
    
    #-----------------------------+
    # SUBFUNCTION VARS            |
    #-----------------------------+
    my ( $source, $findmite_in, $gff_out, $src_suffix, $append_gff, 
	 $fm_seqname, $fasta_out, $append_fasta) = @_;

    # GFF OUT VARS
    my $gff_seq_id;                 # Seq id for use in gff
    my $gff_start;                  # Start of the feature
    my $gff_end;                    # End of the feature
    my $gff_source;                 # Source ie findmite_at_11
    my $gff_name;                   # Name of the feature 

    # FINDMITE HEADER VARS
    my $fm_dir_repeats;             # Direct repeats
    my $fm_tir_len;                 # Lenght of the TIR
    my $fm_mismatch_max;            # Max number of mistmatches
    my $fm_filter_at;               # Were A/T strings filtered 
    my $fm_filter_cg;               # Were C/G strings filtered
    my $fm_filter_atta;             # Were AT/TA strings filtered
    my $fm_filter_2_base;           # Percent (0 to 100)
    my $fm_min_dist;                # Minimum distance used
    my $fm_max_dist;                # Maximum distance used
    my $fm_file_an;                 # Name of the file analyzed

    # FINDMATE INDIVIDUAL MITE VARS
    my $fm_pattern;                 # Pattern
    my $fm_seq_id;                  # Id of the query sequence
    my $fm_num_mismatch;            # Num of mismatches
    my $fm_seq = "";                # Sequence string as parsed from findmite
    my $mite_seq_string;            # Sequence string of the putatitve mite
    my $mite_context_5;             # 5' Context of the mite (40bp)
    my $mite_context_3;             # 3' Context of the mite (40bp)
    my $mite_start;                 # Start of the mite as parsed from FM
    my $mite_end;                   # End of the mite as parsed from FM

    # Counter
    my $mite_num = 0;               # Incremented ID number for the mite

    # BOOLEANS
    my $in_seq = 0;                 # Boolean, in seq data (past header info)

    #-----------------------------+
    # OPEN INFILE                 |
    #-----------------------------+
    if ($findmite_in) {
	open (INFILE, "<$findmite_in") ||
	    die "Can not open input file:\n$findmite_in\n";
    }
    else {
	print STDERR "Expecting input from STDIN\n";
	open (INFILE, "<&STDIN") ||
	    die "Could not open STDIN for input.\n"
    }
    
    #-----------------------------+
    # OPEN THE GFF OUTFILE        |
    #-----------------------------+
    if ($gff_out) {
	if ($append_gff) {
	    open (GFFOUT, ">>$gff_out") ||
		die "Can not open output file:\n$gff_out\n";	
	}
	else {
	    open (GFFOUT, ">$gff_out") ||
		die "Can not open gff output file:\n$gff_out\n";
	}
    }
    else {
	open (GFFOUT, ">&STDOUT") ||
	    die "Can not print to STDOUT\n";
    }

    #-----------------------------+
    # FASTA OUTPUT FILE           |
    #-----------------------------+
    if ($fasta_out) {
	if ($append_fasta) {
	    open (FASTAOUT, ">>$fasta_out") ||
		die "Can not FASTA output file:\n$fasta_out\n";
	} 
	else {
	    open (FASTAOUT, ">$fasta_out") ||
	    die "Can not FASTA output file:\n$fasta_out\n";
	}
    } 

    #-----------------------------+
    # PROCESS FINDMITE FILE       |
    #-----------------------------+
    while (<INFILE>) {
	chomp;
 	#print $_."\n";

	if(m/^Pattern  : (.*)/) {

	    # If we have previously loaded seq data then
	    # futher parse the sequence string and and
	    # print the results to the gff output
	    if ($in_seq) {
		                   
		if ($fm_seq =~ m/(.*)\((.*)\)\-{5}(.*)\-{5}\((.*)\)(.*)/) {
		    $mite_context_5 = $1;
		    $mite_start = $2;
		    $mite_seq_string = $3;
		    $mite_end = $4;
		    $mite_context_3 = $5;

		    if ($verbose) {
			print STDERR "\t5CON: $mite_context_5\n";
			print STDERR "\tSTAR: $mite_start\n";
			print STDERR "\tMITE: $mite_seq_string\n";
			print STDERR "\tEND : $mite_end\n";
			print STDERR "\t3CON: $mite_context_3\n";
			print STDERR "\n\n";
		    }

		    #-----------------------------+
		    # SET SEQUENCE ID             |
		    #-----------------------------+
		    if ($fm_seqname) {
			$gff_seq_id = $fm_seqname;
		    }
		    else {
			if ($fm_seq_id =~ m/^>(\S*)\s./) {
			    $gff_seq_id = $1;
			}
			elsif ($fm_seq_id =~ m/^>(.*)/) {
			    $gff_seq_id = $1;
			}
			else {
			    $gff_seq_id = $fm_seq_id;
			}
		    }
		    
		    #-----------------------------+
		    # SET PROGRAM SOURCE           |
		    #-----------------------------+
		    $gff_source = $source;
		    if ($src_suffix) {
			$gff_source = $gff_source.":".$src_suffix;
		    }
#		    else {
#			$gff_source = $gff_source.":".
#			    $fm_dir_repeats."_".$fm_tir_len;
#		    }
		    
		    #-----------------------------+
		    # SET FEATURE NAME            |
		    #-----------------------------+
		    if ($src_suffix) {
			$gff_name = "findmite_".$src_suffix.
			    "_".$mite_num;
		    }
		    else {
			$gff_name = "findmite_".$mite_num;
		    }

		    $gff_start = $mite_start;
		    $gff_end = $mite_end;

		    #-----------------------------+
		    # PRINT TO GFF FILE           |
		    #-----------------------------+
		    my $attribute;
		    if ($gff_ver =~ "GFF3") {
			$attribute = "ID=".$gff_name;
		    }
		    else {
			$attribute = $gff_name
		    }

		    print GFFOUT "$gff_seq_id\t".   # Seq name
			"$gff_source\t".            # Source
			"MITE\t".                   # Feature type
			"$gff_start\t".             # Start
			"$gff_end\t".               # End
			".\t".                      # Score
			"+\t".                      # Strand
			".\t".                      # Frame
			"$attribute\n";              # Feature name

		    #-----------------------------+
		    # PRINT MITE TO FASTA FILE    |
		    #-----------------------------+
		    if ($fasta_out) {
			print FASTAOUT ">$gff_name\n";
			print FASTAOUT "$mite_seq_string\n";
		    }

		}

		# Reset vals to null
		$fm_seq = "";
	    }

	    $in_seq = 1;
	    $mite_num++;
	    $fm_pattern = $1;
	    print STDERR "$fm_pattern\n" if $verbose;
	}
	elsif(m/^Sequence : (.*)/) {
	    $fm_seq_id = $1;
	    print STDERR "\t$fm_seq_id\n" if $verbose;
	}
	elsif(m/^Mismatch : (.*)/){
	    $fm_num_mismatch = $1;
	    print STDERR "\t$fm_num_mismatch\n" if $verbose;
	}
	elsif($in_seq) {
	    $fm_seq = $fm_seq.$_;
	}
	#-----------------------------+
	# HEADER INFORMATION          | 
	#-----------------------------+
	elsif(m/Direct repeats.{17}(.*)/) {
	    $fm_dir_repeats = $1;
	}
	elsif(m/Length of TIR.{18}(.*)/) {
	    $fm_tir_len = $1;
	}
	elsif(m/Number of mis.{18}(.*)/) {
	    $fm_num_mismatch = $1;
	}
	elsif(m/Filtering A\/T.{18}(.*)/) {
	    $fm_filter_at = $1;
	}
	elsif(m/Filtering C\/G.{18}(.*)/) {
	    $fm_filter_cg = $1;
	}
	elsif(m/Filtering AT\/.{18}(.*)/) {
	    $fm_filter_atta = $1;
	}
	elsif(m/Filtering 2.{20}(.*)/) {
	    $fm_filter_2_base = $1;
	}
	elsif(m/Minimum dist.{19}(.*)/) {
	    $fm_min_dist = $1;
	}
	elsif(m/Maximum dist.{19}(.*)/) {
	    $fm_max_dist = $1;
	}
	elsif(m/The results from the input.{17}(.*)/) {
	    $fm_file_an = $1;
	}

    }

    
    #-----------------------------+
    # PRINT OUT OF DATA IN VARS   |
    #-----------------------------+
    if ($fm_seq =~ m/(.*)\((.*)\)\-{5}(.*)\-{5}\((.*)\)(.*)/) {
    #if ($fm_seq =~ m/(.*)\((.*)\)\-\-\-\-\-(.*)\-\-\-\-\-\((.*)\)(.*)/) {
	$mite_context_5 = $1;
	$mite_start = $2;
	$mite_seq_string = $3;
	$mite_end = $4;
	$mite_context_3 = $5;
	

	if ($verbose) {
	    print STDERR "\t5CON: $mite_context_5\n";
	    print STDERR "\tSTAR: $mite_start\n";
	    print STDERR "\tMITE: $mite_seq_string\n";
	    print STDERR "\tEND : $mite_end\n";
	    print STDERR "\t3CON: $mite_context_3\n";
	    print STDERR "\n\n";
	}	

    }



    #////////////////////////////////
    # THE FOLLOWING CODE SHOULD NOT BE THIS REDUNDANT CUT AND COPY
    #-----------------------------+
    # SET SEQUENCE ID             |
    #-----------------------------+
    if ($fm_seqname) {
	$gff_seq_id = $fm_seqname;
    }
    else {
	if ($fm_seq_id =~ m/^>(\S*)\s./) {
	    $gff_seq_id = $1;
	}
	elsif ($fm_seq_id =~ m/^>(.*)/) {
	    $gff_seq_id = $1;
	}
	else {
	    $gff_seq_id = $fm_seq_id;
	}
    }
    
    #-----------------------------+
    # SET PROGRAM SOURCE           |
    #-----------------------------+
    $gff_source = $source;
    if ($src_suffix) {
	$gff_source = $gff_source.":".$src_suffix;
    }
#    else {
#	$gff_source = $gff_source.":".
#	    $fm_dir_repeats."_".$fm_tir_len;
#    }
    
    #-----------------------------+
    # SET FEATURE NAME            |
    #-----------------------------+
    if ($src_suffix) {
	$gff_name = "findmite_".$src_suffix.
	    "_".$mite_num;
    }
    else {
	$gff_name = "findmite_".$mite_num;
    }
    
    $gff_start = $mite_start;
    $gff_end = $mite_end;
    


    my $attribute;
    if ($gff_ver =~ "GFF3") {
	$attribute = "ID=".$gff_name;
    }
    else {
	$attribute = $gff_name
    }
    
    print GFFOUT "$gff_seq_id\t".   # Seq name
	"$gff_source\t".            # Source
	"MITE\t".                   # Feature type
	"$gff_start\t".             # Start
	"$gff_end\t".               # End
	".\t".                      # Score
	"+\t".                      # Strand
	".\t".                      # Frame
	"$attribute\n";              # Feature name


#    print GFFOUT "$gff_seq_id\t".   # Seq name
#	"$gff_source\t".            # Source
#	"MITE\t".                   # Feature type
#	"$gff_start\t".             # Start
#	"$gff_end\t".               # End
#	".\t".                      # Score
#	"+\t".                      # Strand
#	".\t".                      # Frame
#	"$gff_name\n";              # Feature name
#
#    # END REDUNDANT CUT AND COPY
#    #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

    #-----------------------------+
    # PRINT MITE TO FASTA FILE    |
    #-----------------------------+
    if ($fasta_out) {
	print FASTAOUT ">$gff_name\n";
	print FASTAOUT "$mite_seq_string\n";
    }

    # Print while debu
    if ($verbose) {
	print STDERR "FILE  : $fm_file_an\n";
	print STDERR "DIRREP: $fm_dir_repeats\n";
	print STDERR "TIRLEN: $fm_tir_len\n";
	print STDERR "NUMMIS: $fm_num_mismatch\n";
	print STDERR "F_AT  : $fm_filter_at\n";
	print STDERR "F_CG  : $fm_filter_cg\n";
	print STDERR "F_ATTA: $fm_filter_atta\n";
	print STDERR "F_2BAS: $fm_filter_2_base\n";
	print STDERR "MIN   : $fm_min_dist\n";
	print STDERR "MAX   : $fm_max_dist\n";
    }

    close INFILE;
    close GFFOUT;
    close FASTAOUT if ($fasta_out);

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

cnv_findmite2gff.pl - Convert FINDMITE output to gff format

=head1 VERSION

This documentation refers to program $Rev: 948 $

=head1 SYNOPSIS

=head2 Usage

    cnv_findmite2gff -i InFile -o GffOutFile [-f FastaOutFile]

=head2 Required Options

    -i   # Path to the Findmite result to convert
         # If not specified, the program will expect input from STDIN
    -o   # Path to the gff format output file
         # If not specified the program will write output to STDOUT

=head1 DESCRIPTION

Converts FINDMITE output to GFF format. May also be used to extract the
MITEs from to a FASTA file by using the --fasta option.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path of the input file that contains the results of the FINDMITE program.
If this option is not specified, the program
will expect input from STDIN.

=item -o,--outfile

Path to the gff output file. If and outfile file path is not specified, the
program will write the output to STDERR.

=back

=head1 OPTIONS

=over 2

=item -p,--param

The parameter name to append to the source program name. This information will
be appended to the second column of the gff output file. This is used to 
specify the parameter set used to generate the FINDMITE results.

=item --program

The program name to use. This is the data in the second column of the 
gff output file. Be default, this is set to 'findmite'. This option
allows you to specify other program names if desired.

=item -n,--name

Identifier for the sequence file that was masked with repeatmasker. The
out file from repeatmasker may have truncated your original file name,
and this option allows you to use the full sequence name.

=item -f, --fasta

Genreate a fasta format output file. This file will contain the putative MITEs.

=item --append-gff

Append the gff data to the any existing data in the file specified
by -o,--outfile.

=item --append-fasta

Append the MITE fasta file data to any existing data in the file
specified by -f,fasta.

=item -q,--quiet

Run the program with minimal output.

=item --verbose

Run the program in verbose mode.

=item --usage

Short overview of how to use program from command line.

=item --help

Show program usage with summary of options.

=item --version

Show program version.

=item --man

Show the full program manual. This uses the perldoc command to print the 
POD documentation for the program.

=back

=head1 DIAGNOSTICS

The error messages that can be generated will be listed here.

=over 2

=item * Expecting Input from STDIN

You will see this message if you did not specify an input file with the -i
or --input options.

=back

=head1 CONFIGURATION AND ENVIRONMENT

This program does not make use of a configuration file or variables set 
in the user environment.

=head1 DEPENDENCIES

Other modules or software that the program is dependent on.

=head2 Software

=over

=item * FINDMITE.

This program requires the FINDMITE program. This program is available as
an executable for Unix based operating systems at:
at http://jaketu.biochem.vt.edu/dl_software.htm

=back

=head2 Perl Modules

This program does not make use of Perl modules outside of the normal suite
of modules present in a typical installation of perl.

=head1 BUGS AND LIMITATIONS

=head2 Bugs

=over 2

=item * No bugs currently known 

If you find a bug with this software, file a bug report on the DAWG-PAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=204962

=back

=head2 Limitations

=over

=item * No Major Limitations Currently Known

If you discover limitations with your use of this software, please
file a bug 
report on the DAWG-PAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=204962

=back

=head1 LICENSE

GNU General Public License, Version 3

L<http://www.gnu.org/licenses/gpl.html>

THIS SOFTWARE COMES AS IS, WITHOUT ANY EXPRESS OR IMPLIED
WARRANTY. USE AT YOUR OWN RISK.

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 08/30/2007

UPDATED: 02/26/2010

VERSION: $Rev: 948 $

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
# 
# 08/30/2007
# - Program started
# 08/31/2007
# - Added fasta output
# - Updated POD documentation
# 02/03/2009
# - Switched to print_help subfunction to extract help
#   from POD documentation
# - Updated POD documentation
# - Added support for specifying sequence source with
#   the -n,--name option
# - Added support for specifying parameter set with the
#   -p, --param option
# - Added support for specifying program source with the
#   --program option
# 02/26/2010
# - Added option for GFF3 format
# - Changed type from mite to SO compliant MITE
