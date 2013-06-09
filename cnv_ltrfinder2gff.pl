#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_ltrfinder2gff.pl - Converts ltr_finder to gff         |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 09/14/2007                                       |
# UPDATED: 01/29/2010                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Converts the LTR_FINDER results to gff format.           |
#                                                           |
# VERSION: $Rev: 948 $                                      |
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
use Cwd;                       # Get the current working directory
use File::Copy;                # Copy files

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev: 948 $ =~ /(\d+)/;

# Get GFF version from environment, GFF2 is DEFAULT
my $gff_ver = uc($ENV{DP_GFF}) || "GFF2";

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $infile;                    # Infile. textfile result from LTR_FINDER
my $outfile;                   # Outfile.

# Booleans
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $append = 0;                # Append gff output to existing file
my $param;                    # Suffix appended to the end of the gff 
my $seqname;          #
my $program = "ltr_finder";   # The source program

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
		    "append"      => \$append,
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

if ($show_man) {
    # User perldoc to generate the man documentation.
    system ("perldoc $0");
    exit($ok ? 0 : 2);
}

if ($show_version) {
    print "\ncnv_ltrfinder2gff.pl:\n".
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

ltrfinder2gff ($program, $seqname, $infile, $outfile, $append, $param);

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

sub ltrfinder2gff {
    
    # seq_id could be extracted from the ltrfinder result
    # but passing it to the subfunction directly allows for cases
    # where the id assigned by ltr_struc differs from the way
    # the user is referring to the assembly
    # MAY WANT TO ALLOW FOR USING THE $ls_seq_id
    my ($gff_src, $seq_id, $lf_infile, $gffout, $do_append, $gff_suffix) = @_;


    #////////////////////////////////////////
    #////////////////////////////////////////
    #////////////////////////////////////////
    # ltr finder results as an array
    # new 01/28/2010
    my @ltr_results;
    # Starting i at -1 so that increments start at 0
    my $i = -1;
    my $j = -1;
    #////////////////////////////////////////
    #////////////////////////////////////////
    #////////////////////////////////////////
    
    # The gff src id
    #my $gff_src = "ltr_finder";
    if ($gff_suffix) {
	$gff_src = $gff_src.":".$gff_suffix;
    }

    my $print_gff_out = 0;  # Boolean to print out gff data

    my $gff_str_out;        # A single out string line of gff out

    # LTF FINDER COUTNERS/INDICES
    my $lf_id_num = 0;       # Incremented for every predicted model

    # LTR_FINDER BOOLEANS
    my $in_emp_details = 0;  # Exact match pairs details
    my $in_ltr_align = 0;    # Details of the LTR alignment
    my $in_pbs_align = 0;
    my $in_ppt;
    
    my $lf_prog_name;             # LTR Finder program name
    my $lf_seq_id;                # Seq id
    my $lf_seq_len;               # Sequence length
    my $lf_version;               # Version of LTR finder being parsed
    my $lf_trna_db;               # tRNA database used
    
    
    # Status strings
    my $has_5ltr_tg;                 # TG in 5' END of 5' LTR
    my $has_5ltr_ca;                 # CA in 3' END of 5' LTR
    my $has_3ltr_tg;                 # TG in 5' END of 3' LTR
    my $has_3ltr_ca;                 # CA in 3' END of 3' LTR
    my $has_tsr;                     # Has Target Site Replication
    my $has_pbs;                     # Has Primer Binding Site
    my $has_ppt;                     # Has Poly Purine Tract
    my $has_rt;                      # Has Reverse Transcriptase
    my $has_in_core;                 # Has Integrase Core
    my $has_in_cterm;                # Has Integrase C-term
    my $has_rh;                      # Has RNAseH


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

    #-----------------------------+
    # OPEN INPUT FILE             |
    #-----------------------------+
    if ($lf_infile) {
	open (INFILE, "<$lf_infile") ||
	    die "ERROR: Can not open LTR_FINDER result file\n $lf_infile\n";
	
    }
    else {
	print STDERR "Expecting input from STDIN\n";
	open (INFILE, "<&STDIN") ||
	    die "Can not accept input from standard input.\n";
    }


    
    while (<INFILE>) {
	chomp;
	# 
	if (m/Nothing Header(.*)/) {
	    
	}
	

	# IN NEW REC, GET ID
	elsif (m/^\[(.*)\]/) {

	    #-----------------------------+
	    # PRINT STORED GFF OUTPUT     |
	    # IF VALUES ARE PRESENT       |
	    #-----------------------------+
	    
	    # Increment $I and initialize j for protein domains
	    $i++;
	    $j = -1;

	    # override seq id if one passed 
	    if ($seqname) {
		$lf_seq_id = $seqname;
	    }
	    
#	    # FULL SPAN
#	    if ($lf_span_start) {
#		$gff_str_out = "$lf_seq_id\t". # Seq ID
#		    "$gff_src\t".                # Source
#		    "LTR_retrotransposon\t".     # Data type
#		    "$lf_span_start\t".          # Start
#		    "$lf_span_end\t".            # End
#		    "$lf_score\t".               # Score
#		    "$lf_strand\t".              # Strand
#		    ".\t".                       # Frame
#		    "ltr_finder_$lf_ltr_id\n";   # Retro ID
#		print GFFOUT $gff_str_out;
#	    }
#	    
#	    if ($lf_5ltr_start) {
#		$gff_str_out = "$lf_seq_id\t".  # Seq ID
#		    "$gff_src\t".               # Source
#		    "five_prime_LTR\t".         # Data type
#		    "$lf_5ltr_start\t".         # Start
#		    "$lf_5ltr_end\t".           # End
#		    "$lf_score\t".              # Score
#		    "$lf_strand\t".             # Strand
#		    ".\t".                      # Frame
#		    "ltr_finder_$lf_ltr_id\n";  # Retro ID
#		print GFFOUT $gff_str_out;
#	    }
#	    
#	    if ($lf_3ltr_start) {
#		$gff_str_out = "$lf_seq_id\t".  # Seq ID
#		    "$gff_src\t".               # Source
#		    "three_prime_LTR\t".        # Data type
#		    "$lf_3ltr_start\t".         # Start
#		    "$lf_3ltr_end\t".           # End
#		    "$lf_score\t".              # Score
#		    "$lf_strand\t".             # Strand
#		    ".\t".                      # Frame
#		    "ltr_finder_$lf_ltr_id\n";  # Retro ID
#		print GFFOUT $gff_str_out;
#	    }
#	    
#	    if ($has_pbs) {
#		$gff_str_out = "$lf_seq_id\t".  # Seq ID
#		    "$gff_src\t".               # Source
#		    "primer_binding_site\t".    # Data type
#		    "$lf_pbs_start\t" .         # Start
#		    "$lf_pbs_end\t".            # End
#		    "$lf_score\t".              # Score
#		    "$lf_strand\t".             # Strand
#		    ".\t".                      # Frame
#		    "ltr_finder_$lf_ltr_id\n";  # Retro ID
#		print GFFOUT $gff_str_out;
#	    }
#	    
#	    
#	    if ($has_ppt) {
#		$gff_str_out = "$lf_seq_id\t".  # Seq ID
#		    "$gff_src\t".               # Source
#		    "RR_tract\t".               # Data type
#		    "$lf_ppt_start\t".          # Start
#		    "$lf_ppt_end\t".            # End
#		    "$lf_score\t".              # Score
#		    "$lf_strand\t".             # Strand
#		    ".\t".                      # Frame
#		    "ltr_finder_$lf_ltr_id\n";  # Retro ID
#		print GFFOUT $gff_str_out;
#	    }
#	    
#	    
#	    if ($has_tsr) {
#		
#		$gff_str_out = "$lf_seq_id\t".   # Seq ID
#		    "$gff_src\t".                # Source
#		    "target_site_duplication\t". # Data type
#		    "$lf_5tsr_start\t".          # Start
#		    "$lf_5tsr_end\t".            # End
#		    "$lf_score\t".               # Score
#		    "$lf_strand\t".              # Strand
#		    ".\t".                       # Frame
#		    "ltr_finder_$lf_ltr_id\n";   # Retro ID
#		print GFFOUT $gff_str_out;	    
#		
#		$gff_str_out = "$lf_seq_id\t".   # Seq ID
#		    "$gff_src\t".                # Source
#		    "target_site_duplication\t". # Data type
#		    "$lf_3tsr_start\t".          # Start
#		    "$lf_3tsr_end\t".            # End
#		    "$lf_score\t".               # Score
#		    "$lf_strand\t".              # Strand
#		    ".\t".                       # Frame
#		    "ltr_finder_$lf_ltr_id\n";   # Retro ID
#		print GFFOUT $gff_str_out;
#		
#	    }
#	    
#	    
#	    # Integrase Core
#	    if ($has_in_core) {
#		#/////////
#		# NOT SONG
#		#\\\\\\\\\
#		$gff_str_out = "$lf_seq_id\t".     # Seq ID
#		    "$gff_src\t".              # Source
#		    "integrase_core_domain\t".  # Data type
#		    "$lf_in_core_dom_start\t".   # Start
#		    "$lf_in_core_dom_end\t".     # End
#		    "$lf_score\t".               # Score
#		    "$lf_strand\t".              # Strand
#		    ".\t".                       # Frame
#		    "ltr_finder_$lf_ltr_id\n";   # Retro ID
#		print GFFOUT $gff_str_out;
#		
#		#/////////
#		# NOT SONG
#		#\\\\\\\\\
#		$gff_str_out = "$lf_seq_id\t".     # Seq ID
#		    "$gff_src\t".              # Source
#		    "integrase_core_orf\t".      # Data type
#		    "$lf_in_core_orf_start\t".   # Start
#		    "$lf_in_core_orf_end\t".     # End
#		    "$lf_score\t".               # Score
#		    "$lf_strand\t".              # Strand
#		    ".\t".                       # Frame
#		    "ltr_finder_$lf_ltr_id\n";   # Retro ID
#		print GFFOUT $gff_str_out;
#		
#	    } # End of has in_core
#	    
#	    
#	    if ($has_in_cterm) {
#		
#		#/////////
#		# NOT SONG
#		#\\\\\\\\\
#		$gff_str_out = "$lf_seq_id\t".     # Seq ID
#		    "$gff_src\t".                # Source
#		    "integrase_cterm_domain\t".  # Data type
#		    "$lf_in_cterm_dom_start\t".  # Start
#		    "$lf_in_cterm_dom_end\t".    # End
#		    "$lf_score\t".               # Score
#		    "$lf_strand\t".              # Strand
#		    ".\t".                       # Frame
#		    "ltr_finder_$lf_ltr_id\n";   # Retro ID
#		print GFFOUT $gff_str_out;
#		
#		$gff_str_out = "$lf_seq_id\t".     # Seq ID
#		    "$gff_src\t".                # Source
#		    "integrase_cterm_orf\t".     # Data type
#		    "$lf_in_cterm_orf_start\t".  # Start
#		    "$lf_in_cterm_orf_end\t".    # End
#		    "$lf_score\t".               # Score
#		    "$lf_strand\t".              # Strand
#		    ".\t".                       # Frame
#		    "ltr_finder_$lf_ltr_id\n";   # Retro ID
#		print GFFOUT $gff_str_out;
#		    
#	    } # End of has_in_cterm
#	    
#	    if ($has_rh) {
#		
#		#/////////
#		# NOT SONG
#		#\\\\\\\\\
#		$gff_str_out = "$lf_seq_id\t".     # Seq ID
#		    "$gff_src\t".                # Source
#		    "rnaseh_domain\t".           # Data type
#		    "$lf_rh_dom_start\t".        # Start
#		    "$lf_rh_dom_end\t".          # End
#		    "$lf_score\t".               # Score
#		    "$lf_strand\t".              # Strand
#		    ".\t".                       # Frame
#		    "ltr_finder_$lf_ltr_id\n";   # Retro ID
#		print GFFOUT $gff_str_out;
#		
#		$gff_str_out = "$lf_seq_id\t".     # Seq ID
#		    "$gff_src\t".                # Source
#		    "rnaseh_orf\t".              # Data type
#		    "$lf_rh_orf_start\t".        # Start
#		    "$lf_rh_orf_end\t".          # End
#		    "$lf_score\t".               # Score
#		    "$lf_strand\t".              # Strand
#		    ".\t".                       # Frame
#		    "ltr_finder_$lf_ltr_id\n";   # Retro ID
#		print GFFOUT $gff_str_out;
#		
#	    }
#	    
#	    if ($has_rt) {
#		
#		#/////////
#		# NOT SONG
#		#\\\\\\\\\
#		$gff_str_out = "$lf_seq_id\t".     # Seq ID
#		    "$gff_src\t".                # Source
#		    "rt_domain\t".               # Data type
#		    "$lf_rt_dom_start\t".        # Start
#		    "$lf_rt_dom_end\t".          # End
#		    "$lf_score\t".               # Score
#		    "$lf_strand\t".              # Strand
#		    ".\t".                       # Frame
#		    "ltr_finder_$lf_ltr_id\n";   # Retro ID
#		print GFFOUT $gff_str_out;
#		
#		$gff_str_out = "$lf_seq_id\t".     # Seq ID
#		    "$gff_src\t".                # Source
#		    "rt_orf\t".                  # Data type
#		    "$lf_rt_orf_start\t".        # Start
#		    "$lf_rt_orf_end\t".          # End
#		    "$lf_score\t".               # Score
#		    "$lf_strand\t".              # Strand
#		    ".\t".                       # Frame
#		    "ltr_finder_$lf_ltr_id\n";   # Retro ID
#		print GFFOUT $gff_str_out;
#		
#	    } # End of has_reverse_transcriptase
	    
	    #-----------------------------+
	    # LOAD ID VAR                 |
	    #-----------------------------+
#	    $lf_ltr_id = $1;
	    $ltr_results[$i]{lf_ltr_id} = $1;


	}
	
	# SEQ ID AND LENGTH
	elsif (m/>Sequence: (.*) Len:(.*)/){
	    # Incrmenting here would be the number of contigs that are being processed

	    $lf_seq_id = $1;
	    $lf_seq_len = $2;

	}
	
	# SPAN LOCATION, LENGTH, AND STRAND
	elsif (m/^Location : (\d*) - (\d*) Len: (\d*) Strand:(.)/){
#	    $lf_span_start = $1;
#	    $lf_span_end = $2;
#	    $lf_length = $3;
#	    $lf_strand = $4;

	    $ltr_results[$i]{lf_span_start} = $1;
	    $ltr_results[$i]{lf_span_end} = $2;
	    $ltr_results[$i]{lf_length} = $3;
	    $ltr_results[$i]{lf_strand} = $4;
	    $ltr_results[$i]{lf_seq_id} = $lf_seq_id;

	}
	
	# SCORE SIMILARITY
	elsif (m/^Score    : (.*) \[LTR region similarity:(.*)\]/){
	    $ltr_results[$i]{lf_score} = $1;            # Score
	    $ltr_results[$i]{lf_ltr_similarity} = $2;   # Similarity of LTRs
	}
	
	# STATUS SET
	elsif (m/^Status   : (\d)(\d)(\d)(\d)(\d)(\d)(\d)(\d)(\d)(\d)(\d)/){
	    # Since this is a binary string, it can be split as digits
	    # and used to load the $has_* booleans
	    $has_5ltr_tg = $1;
	    $has_5ltr_ca = $2;
	    $has_3ltr_tg = $3;
	    $has_3ltr_ca = $4;
	    $has_tsr = $5;
	    $has_pbs = $6;
	    $has_ppt = $7;
	    $has_rt = $8;
	    $has_in_core = $9;
	    $has_in_cterm = $10;
	    $has_rh = $11;
	}
	
	# 5' LTR
	elsif (m/^5\'-LTR   : (\d*) - (\d*) Len: (\d*)/){
	    $ltr_results[$i]{lf_5ltr_start} = $1;
	    $ltr_results[$i]{lf_5ltr_end} = $2;
	    $ltr_results[$i]{lf_5ltr_len} = $3;
	}
	
	# 3' LTR
	elsif (m/^3\'-LTR   : (\d*) - (\d*) Len: (\d*)/){
	    $ltr_results[$i]{lf_3ltr_start} = $1;
	    $ltr_results[$i]{lf_3ltr_end} = $2;
	    $ltr_results[$i]{lf_3ltr_len} = $3;
	}
	
    # TARGET SITE REPLICATION
	elsif (m/TSR      : (\d*) - (\d*) , (\d*) - (\d*) \[(.*)\]/){
	    $ltr_results[$i]{lf_5tsr_start} = $1;
	    $ltr_results[$i]{lf_5tsr_end} = $2;
	    $ltr_results[$i]{lf_3tsr_start} = $3;
	    $ltr_results[$i]{lf_3tsr_end} = $4;
	    $ltr_results[$i]{lf_tsr_string} = $5;
	}
	
	# SHARPNESS METRIC
	elsif (m/^Sharpness: (.*),(.*)/){
	    $ltr_results[$i]{lf_sharp5} = $1;
	    $ltr_results[$i]{lf_sharp_3} = $2;
	}
	
	# PBS
	elsif (m/PBS   : \[(\d*)\/(\d*)\] (\d*) - (\d*) \((.*)\)/) {
	    $ltr_results[$i]{lf_pbs_num_match} = $1;  # Number of matching bases
	    $ltr_results[$i]{lf_pbs_aln_len} = $2;    # PBS alignment length
	    $ltr_results[$i]{lf_pbs_start} = $3;      # Start of PBS signal
	    $ltr_results[$i]{lf_pbs_end} = $4;        # End of PBS signal
	    $ltr_results[$i]{lf_pbs_trna} = $5;       # PBS tRNA type and anti-codon
	}
	
	# PPT
	elsif (m/PPT   : \[(\d*)\/(\d*)\] (\d*) - (\d*)/) {
	    $ltr_results[$i]{lf_ppt_num_match} = $1;
	    $ltr_results[$i]{lf_ppt_aln_len} = $2;
	    $ltr_results[$i]{lf_ppt_start} = $3;
	    $ltr_results[$i]{lf_ppt_end} = $4;
	}
	
	# PROTEIN DOMAINS
	# This will need to be modified and checked after another run
	# using ps_scan to get the additional domains
	#
	#Domain: 56796 - 57326 [possible ORF:56259-59144, (IN (core))]
	elsif (m/Domain: (\d*) - (\d*) \[possible ORF:(\d*)-(\d*), \((.*)\)\]/) {
	    

	    my $lf_domain_name = $5;

	    #///////////////////////////
	    # PUSH ALL DOMAIN DATA TO HREF
	    #//////////////////////////
	    # Increment j for index of protein domains
	    $j++;

	    $ltr_results[$i]{lf_dom}[$j]{lf_dom_start} = $1;
	    $ltr_results[$i]{lf_dom}[$j]{lf_dom_end} = $2;
	    # 3 and 5 are potential ORFS
	    $ltr_results[$i]{lf_dom}[$j]{lf_dom_name} = $5;

	    #-----------------------------+
	    # NORMALIZE DOMAIN NAME       |
	    #-----------------------------+
	    
	    if ($lf_domain_name =~ 'IN \(core\)') {
		$ltr_results[$i]{lf_dom}[$j]{lf_dom_name_norm} = "Integrase";
		$ltr_results[$i]{lf_dom}[$j]{lf_dom_name_short} = "int";
	    }
	    elsif ($lf_domain_name =~ 'IN \(c-term\)') {
		$ltr_results[$i]{lf_dom}[$j]{lf_dom_name_norm} = "C-terminal";
		$ltr_results[$i]{lf_dom}[$j]{lf_dom_name_short} = "cterm";
	    }
	    elsif ($lf_domain_name =~ 'RH') {
		$ltr_results[$i]{lf_dom}[$j]{lf_dom_name_norm} = "RnaseH";
		$ltr_results[$i]{lf_dom}[$j]{lf_dom_name_short} = "rh";
	    }
	    elsif ($lf_domain_name =~ 'RT') {
		$ltr_results[$i]{lf_dom}[$j]{lf_dom_name_norm} = "Reverse Transcriptase";
		$ltr_results[$i]{lf_dom}[$j]{lf_dom_name_short} = "rvt";
	    }
	    else {
#		print "\a";
#		print STDERR "Unknown domain type: $lf_domain_name\n";
		$ltr_results[$i]{lf_dom}[$j]{lf_dom_name_norm} = $lf_domain_name;
		$ltr_results[$i]{lf_dom}[$j]{lf_dom_name_short} = $lf_domain_name;
	    }
	    
	    
	    
	} # End of elsif Domain
	
	#-----------------------------+
	# FILE HEADER INFORMATION     |
	#-----------------------------+
	
	# PROGRAM NAME
	elsif (m/^Program    : (.*)/) {
	    $lf_prog_name = $1;
#	    $ltr_results[$i]{program} = $1;
	}
	
	# PROGRAM VERSION
	elsif (m/^Version    : (.*)/) {
	    $lf_version = $1;
#	    $ltr_results[$i]{lf_version} = $1;
	}



	# MOD HERE
	if ($print_gff_out) {

	}


    }
    
    close INFILE;


    #-----------------------------------------------------------+
    # PRINT OUTPUT FROM THE ARRAY OF HASHES                     |
    #-----------------------------------------------------------+
    for my $href ( @ltr_results ) {

	
	# For GFF2 format attribute is same for all parts
	my $attribute = "ltr_finder_".$href->{lf_ltr_id};
	my $model_num = $href->{lf_ltr_id};
	$model_num = sprintf("%04d", $model_num);

	# parent id
	my $parent_id;

	# If GFF3 format encode seq_ids to follow specifications
	if ($gff_ver =~ "GFF3") {
	    $href->{lf_seq_id} = seqid_encode( $href->{lf_seq_id} );
	    if ($param) {
		$parent_id = "ltr_finder".
		    "_par".$param.
		    "_model".$model_num;
	    }
	    else {
		$parent_id = "ltr_finder".
		    "_model".$model_num;
	    }
	}

	#-----------------------------+
	# LTR retrotransposon span    |
	#-----------------------------+
	if ($gff_ver =~ "GFF3") {
	    $attribute = "ID=".$parent_id;
	}
	$gff_str_out = $href->{lf_seq_id}."\t".    # Seq ID
	    "$gff_src\t".                          # Source
	    "LTR_retrotransposon\t".               # Data type
	    $href->{lf_span_start}."\t".           # Start
	    $href->{lf_span_end}."\t".             # End
	    $href->{lf_score}."\t".                # Score
	    $href->{lf_strand}."\t".               # Strand
	    ".\t".                                 # Frame
	    $attribute.
	    "\n";
	print GFFOUT $gff_str_out;

	#-----------------------------+
	# 5 Prime LTR                 |
	#-----------------------------+
	if ($gff_ver =~ "GFF3") {
	    $attribute = "ID=".$parent_id."_five_prime_LTR".
		";Name=Five Prime LTR".
		";Parent=".$parent_id;
	}
	$gff_str_out = $href->{lf_seq_id}."\t".    # Seq ID
	    "$gff_src\t".                          # Source
	    "five_prime_LTR\t".                    # Data type
	    $href->{lf_5ltr_start}."\t".           # Start
	    $href->{lf_5ltr_end}."\t".             # End
	    $href->{lf_score}."\t".                # Score
	    $href->{lf_strand}."\t".               # Strand
	    ".\t".                                 # Frame
	    $attribute.
	    "\n";
	print GFFOUT $gff_str_out;


	#-----------------------------+
	# 3 Prime LTR                 |
	#-----------------------------+
	if ($gff_ver =~ "GFF3") {
	    $attribute = "ID=".$parent_id."_three_prime_LTR".
		";Name=Three Prime LTR".
		";Parent=".$parent_id;
	}
	$gff_str_out = $href->{lf_seq_id}."\t".    # Seq ID
	    "$gff_src\t".                          # Source
	    "three_prime_LTR\t".                   # Data type
	    $href->{lf_3ltr_start}."\t".           # Start
	    $href->{lf_3ltr_end}."\t".             # End
	    $href->{lf_score}."\t".                # Score
	    $href->{lf_strand}."\t".               # Strand
	    ".\t".                                 # Frame
	    $attribute.
	    "\n";
	print GFFOUT $gff_str_out;

	#-----------------------------+
	# PRIMER BINDING SITE         |
	#-----------------------------+
	if ($gff_ver =~ "GFF3") {
	    $attribute = "ID=".$parent_id."_primer_binding_site".
		";Name=PBS".
		";Parent=".$parent_id;
	}
	$gff_str_out = $href->{lf_seq_id}."\t".    # Seq ID
	    "$gff_src\t".                          # Source
	    "primer_binding_site\t".                   # Data type
	    $href->{lf_pbs_start}."\t".           # Start
	    $href->{lf_pbs_end}."\t".             # End
	    $href->{lf_score}."\t".                # Score
	    $href->{lf_strand}."\t".               # Strand
	    ".\t".                                 # Frame
	    $attribute.
	    "\n";
	print GFFOUT $gff_str_out;

	#-----------------------------+
	# RR_tract                    |
	#-----------------------------+
	if ($gff_ver =~ "GFF3") {
	    $attribute = "ID=".$parent_id."_RR_tract".
		";Name=PPT".
		";Parent=".$parent_id;
	}
	$gff_str_out = $href->{lf_seq_id}."\t".    # Seq ID
	    "$gff_src\t".                          # Source
	    "RR_tract\t".                         # Data type
	    $href->{lf_ppt_start}."\t".           # Start
	    $href->{lf_ppt_end}."\t".             # End
	    $href->{lf_score}."\t".                # Score
	    $href->{lf_strand}."\t".               # Strand
	    ".\t".                                 # Frame
	    $attribute.
	    "\n";
	print GFFOUT $gff_str_out;

	#-----------------------------+
	# 5' TSD                      |
	#-----------------------------+
	if ( $href->{lf_5tsr_start}) {
	    if ($gff_ver =~ "GFF3") {
		$attribute = "ID=".$parent_id."_tsd5".
		    ";Name=Target Site Duplication".
		    ";Parent=".$parent_id;
	    }
	    $gff_str_out = $href->{lf_seq_id}."\t".   # Seq ID
		"$gff_src\t".                         # Source
		"target_site_duplication\t".          # Data type
		$href->{lf_5tsr_start}."\t".           # Start
		$href->{lf_5tsr_end}."\t".             # End
		$href->{lf_score}."\t".               # Score
		$href->{lf_strand}."\t".              # Strand
		".\t".                                # Frame
		$attribute.
		"\n";
	    print GFFOUT $gff_str_out;
	}

	#-----------------------------+
	# 3' TSD                      |
	#-----------------------------+
	if ( $href->{lf_3tsr_start}) {
	    if ($gff_ver =~ "GFF3") {
		$attribute = "ID=".$parent_id."_tsd3".
		    ";Name=Target Site Duplication".
		    ";Parent=".$parent_id;
	    }
	    $gff_str_out = $href->{lf_seq_id}."\t".   # Seq ID
		"$gff_src\t".                         # Source
		"target_site_duplication\t".          # Data type
		$href->{lf_3tsr_start}."\t".          # Start
		$href->{lf_3tsr_end}."\t".            # End
		$href->{lf_score}."\t".               # Score
		$href->{lf_strand}."\t".              # Strand
		".\t".                                # Frame
		$attribute.
		"\n";
	    print GFFOUT $gff_str_out;
	}

	#-----------------------------+
	# PROTEIN DOMAINS             |
	#-----------------------------+
	for my $dom ( @{ $href->{lf_dom} } ) {
	    # The following will just print the name

	    my $type = $dom->{lf_dom_name_norm};

	    if ($gff_ver =~ "GFF3") {
		$attribute = "ID=".$parent_id.$dom->{lf_dom_name_short}.
		    ";Name=".$dom->{lf_dom_name_norm}.
		    ";Parent=".$parent_id;
		#$type = "polypeptide_domain";
		$type = "transposable_element_gene";
		# It may be necessary to first define a "gene" and
		# then set domains..
	    }

	    $gff_str_out = $href->{lf_seq_id}."\t".   # Seq ID
		"$gff_src\t".                         # Source
		$type."\t".
		$dom->{lf_dom_start}."\t".          # Start
		$dom->{lf_dom_end}."\t".            # End
		$href->{lf_score}."\t".               # Score
		$href->{lf_strand}."\t".              # Strand
		".\t".                                # Frame
		$attribute.
		"\n";
	    print GFFOUT $gff_str_out;

	}


    }

    close GFFOUT;


} # End of ltrfinder2gff subfunction


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

cnv_ltrfinder2gff.pl - Converts LTR_Finder output to gff format

=head1 VERSION

This documentation refers to program version $Rev: 948 $

=head1 SYNOPSIS

=head2 Usage

    cnv_ltrfinder2gff.pl -i lf_result.txt -o lf_result.gff

=head2 Required Arguments

    --infile        # Path to the input file
                    # Result from a single record fasta file
    --outdir        # Base output dir

=head1 DESCRIPTION

Convert the ltrfinder output to gff format. This assumes that the output
from ltr_finder correspons to a single BAC. A suffix can be passed with the
--suffix option to provide a suffix for the source column. For example
run default ltr_finder results with --suffix def to create 
ltr_finder:def. T

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path of the input file. If an input file is not provided, the program
will expect input from STDIN.

=item -o,--outfile

Path of the output file.  If an output path is not provided,
the program will write output to STDOUT.

=back

=head1 OPTIONS

=over 2

=item --gff-ver

The GFF version for the output. This will accept either gff2 or gff3 as the
options. By default the GFF version will be GFF2 unless specified otherwise.
The default GFF version for output can also be set in the user environment
with the DP_GFF option. The command line option will always override the option
defined in the user environment.

=item -n,--name

The sequence name to use in the GFF output file. Otherwise, this will
just use 'seq' as the sequence name.

=item -p, --param

The name of the paramter set used. This will be appened to the data in the
second column of the GFF file, and can be used to distinguish among 
parameter combinations
for multiple applications of ltrfinder to the same sequence file.

=item --apend

Append the GFF output to the gff file indicate by the --outfile option. This
allows you to append results from multiple programs to a single GFF file, but
it is generally recommended to create separate GFF files and concatenate
them at a later stage.

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

Error messages that you may encounter and possible solutions are listed below:

=over 2

=item Expecting input from STDIN

If a file is not specified by the -i or --infile option, the program will
expect to receive intput from standard input.

=back

=head1 CONFIGURATION AND ENVIRONMENT

=over 2

=item DP_GFF

The DP_GFF variable can be defined in the user environment to set
the default GFF version output. Valid settings are 'gff2' or
'gff3'.

=back

=head1 DEPENDENCIES

=head2 Software

This program requires the following software:

=over

=item * LTR Finder

This program parses output from the LTR finder program. It is possible to
obtain a linux binary by contacting the authors : xuzh <at> fudan.edu.cn.
It is also possible to obtain these results using the LTR_FINDER web page:
http://tlife.fudan.edu.cn/ltr_finder/

=back

=head2 Perl Modules

This program does not make use of perl modules beyond those installed
with the basic Perl package. If you discover a dependency that is not
documented here, please email the author or file a bug report.

=head1 BUGS AND LIMITATIONS

Any known bugs and limitations will be listed here.

=head2 Bugs

=over 2

=item * No bugs currently known 

If you find a bug with this software, file a bug report on the DAWG-PAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=204962

=back

=head2 Limitations

=over 2

=item * Limited Testing with LTR_FINDER versions

This program is known to parse the results from LTR_FINDER v 1.0.2. This 
program has not been tested with the results from the LTR_FINDER web site.

=back

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

STARTED: 09/14/2007

UPDATED: 01/14/2010

VERSION: $Rev: 948 $

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 09/14/2007
# - Program started
# 10/01/2007
# - Moved main body of program to subfunction ltrfinder2gff
#   This subfunction operates under the assumption that all
#   results in the ltrfinder output file are the results for
#   a single assembly. This will facilitate copying the results
#   to an individual directory for the contig being annotated.
# - Making gff output, making this song complient
# 10/02/2007
# - Finishing gff output, now saving to string to write to
#   both gffout and stdout.
# 01/27/2009
# - Modifying program to use the print_help subfunction that
#   extracts the help message from the POD documentation
# - Updated POD documentation
# - Modified to accept intput from STDIN when --infile not
#   specified at the command line
# - Modified to write output to STOUT when --outfile not
#   specified at the command line
# 01/28/2009
# - Finished update of POD documentation
#
# 03/27/2009
# - Renamed --name to --seqname
# - Fixed program to accept seqname to override parsed name
# - Added program
#
# 04/27/2009
# - modified the parser to work with a multiple record fasta
#   file
# - This was done by removing the statement where I did not
#   print output when the ltr_finder id was set to one
# - I added a logical to check for values, and this will output
#   GFF only when values are present
# 01/14/2010
# - Added sufunctions for encoding GFF3
#   *seqid_encode
#   *gff3_encode
# 01/28/2010
# - Loading variabls to array of hashes
# 01/29/2010
# - Reporting data from array of hashes
