#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_ltrstruc2ann.pl - Convert ltrstruc rpt to annotation  |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 04/14/2008                                       |
# UPDATED: 06/02/2009                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Convert *rpt.txt output files from LTR_STRUC to          |
#  annotation coordinates mapped on the LTR Ret             |
#                                                           |
# VERSION: $Rev$                                            |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+
#
# TO DO: * Add percent divergence of the LTRs to the annotation
#        * Decide of this should be append to an existing gff file or write
#          a new gff for every sequence
#        * Add support for GFF3 output
#

package DAWGPAWS;

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;
use File::Copy;
use Getopt::Long;
# The following needed for printing help
use Pod::Select;               # Print subsections of POD documentation
use Pod::Text;                 # Print POD doc as formatted text file
use IO::Scalar;                # For print_help subfunction
use IO::Pipe;                  # Pipe for STDIN, STDOUT for POD docs
use File::Spec;                # Convert a relative path to an abosolute path
use Bio::Seq;
use Bio::Tools::Run::StandAloneBlast;

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/;

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $e_val = 0.00001;

# Required variables
my $outfile;
my $gff_outfile;
my $repdir;
my $fs_outfile;               # Feature summary outfile
my $outdir;

# Booleans
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_copy = 0;
my $do_seq_data = 0;          # Create files in outdir with sequence data

my $name_root;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
#		    "i|indir=s"     => \$indir,
                    "o|outfile=s"   => \$outfile,
		    "g|gff-out=s"   => \$gff_outfile,
		    "f|feat-sum=s"  => \$fs_outfile,
		    "r|results=s"   => \$repdir,
		    # ADDITIONAL OPTIONS
		    "outdir=s"      => \$outdir,
		    "e|e-val=s"     => \$e_val,
		    "c|copy"        => \$do_copy,
		    "q|quiet"       => \$quiet,
		    "s|seq-data"    => \$do_seq_data,
		    "verbose"       => \$verbose,
		    # ADDITIONAL INFORMATION
		    "usage"         => \$show_usage,
		    "version"       => \$show_version,
		    "man"           => \$show_man,
		    "h|help"        => \$show_help,);

#-----------------------------+
# SHOW REQUESTED HELP         |
#-----------------------------+
if ( ($show_usage) ) {
    print_help ("usage", $0 );
}

if ( ($show_help) || (!$ok) ) {
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
#if ( (!$outfile) || (!$repdir)  ) {
#    print "ERROR: Featsummary file must be specified\n" if !$fs_outfile;
#    print "ERROR: Output file must be specified\n" if !$outfile;
#    print "ERROR: Reports directory must be specified\n" if !$repdir;
#    print "\n";
##    print_help ("usage", $0 );
#}

#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+

#-----------------------------+
# CHECK FOR SLASH IN DIR      |
# VARIABLES                   |
#-----------------------------+
# If the indir does not end in a slash then append one
# TO DO: Allow for backslash
unless ($repdir =~ /\/$/ ) {
    $repdir = $repdir."/";
}

#-----------------------------+
# OPEN FEATURE SUMMARY OUTFILE|
# PRINT HEADER                |
#-----------------------------+
open (FEATSUM, ">$fs_outfile") ||
    die "Can not open the feature summary file\n";
print FEATSUM "#seq_id\t".  # 1
    "fl_rpn_start\t".      # 2
    "fl_rpn_end\t".        # 3
    "ltr_5_start\t".       # 4
    "ltr_5_end\t".         # 5 
    "ltr_5_dn_start\t".    # 6
    "ltr_5_dn_end\t".      # 7
    "ltr_3_start\t".       # 8
    "ltr_3_end\t".         # 9
    "ltr_3_dn_start\t".    # 10
    "ltr_3_dn_end\t".      # 11
    "pbs_start\t".         # 12
    "pbs_end\t".           # 13
    "ppt_start\t".         # 14
    "ppt_end\t".           # 15
    "tsd_5_start\t".       # 16
    "tsd_5_end\t".         # 17
    "tsd_5_seq\t".         # 18
    "tsd_3_start\t".       # 19
    "tsd_3_end\t".         # 20
    "tsd_3_seq\t".         # 21
    "gag_s\tgag_e\t".      # 22 GAG
    "zf_s\tzf_e\t".        # 24 Zinc finger, gag component
    "rvp_s\trvp_e\t".      # 26 Protease/AP
    "rve_s\trve_e\t".      # 28 INTEGRASE
    "rvt_s\trvt_e\t".      # 30 Reverse Transcriptase
    "rh_s\trh_e\t".        # 32 RnaseH
    "chrom_s\tchrom_e\t".  # 34 Chromodomain
    "env_s\tenv_e\t".      # 36Envelope
    "flank_5_seq\t".       # 38 5' Flanking Sequence
    "flank_3_seq\t".       # 39 3' Flanking Sequence
    "strand".              # 40 Strand as reported from ltr_struc
    "\n";        

#-----------------------------+
# Get the LTR_STRUC report    |
# files                       |
#-----------------------------+
opendir(REPDIR,$repdir) 
    || die "Can not open results direcotry:\n$repdir\n";
my @report_files = grep /rprt\.txt$/, readdir REPDIR;
@report_files = sort(@report_files);
closedir (REPDIR); 
# Sort the array of reports

my $num_report_files = @report_files;

if ($verbose) {
    print STDERR "\n-----------------------------------------------\n";
    print STDERR " Report files to process: $num_report_files\n";
    print STDERR "-----------------------------------------------\n\n";
}

my $ind_report_num=0;
for my $ind_report (@report_files) {
    
    
    # This is the id of the sequence the report is for
    #my $ind_report_id = substr ($ind_report,0,$name_root_len);    
    my $ind_report_id = $ind_report;
 
    print STDERR "\tReport: $ind_report\n" if $verbose;
    $ind_report_num++;
    print STDERR "\tProcessing report $ind_report_num".
	" of $num_report_files\n" unless $quiet;
    
    my $seq_id = $ind_report;
    my $gff_out_path = $gff_outfile;
    my $report_file_path = $repdir.$ind_report;

    if ($ind_report_num == 1) {
    # Generall if first record start new gff file
    # 04/14/2008
    # Changed the following to append to existing gff file
	ltrstruc2ann ( $seq_id, $report_file_path,
		       1, $ind_report_num, 
		       $do_seq_data, $gff_out_path);
    }
    else {
	ltrstruc2ann ( $seq_id, $report_file_path,
		       1, $ind_report_num,
		       $do_seq_data, $gff_out_path);
    }
     
} # End for for each individual report


exit;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

# Dropping fasta_in
sub ltrstruc2ann {
    
    # Not renaming the vars ... but gff_out and gff_append both
    # refer to the annotation file
    # VARS PASSED TO THE SUBFUNCTION

#    my ($report_in, $gff_out, $gff_append, $ls_id_num, 
#	$print_seq_data) = @_;

    my ($gff_seq_id, $report_in, $gff_append, $ls_id_num, 
	$print_seq_data, $gff_out) = @_;

# Trying to make $gff_out an option
#    my ($report_in, $gff_append, $ls_id_num, 
#	$print_seq_data) = @_;
#    my ($gff_out) = shift;

    # FASTA RELATED VARS
    my $qry_seq;
    
    # Counters
    my $ltr5_blank_count = 0;

    # LTR STRUC VARS
    my $ls_score;        # Score assigned by LTR_STRUC
    my $ls_contig_len;   # Length of the source contig
    my $ls_orientation;  # Orientation 
    my $ls_retro_len;    # Overall length of the retrotransposon
    my $ls_longest_orf;  # Length of the longest ORF
    my $ls_rt_frame1;
    my $ls_rt_frame2;
    my $ls_rt_frame3;
    my $ls_5ltr_len;     # Length of the 5' LTR
    my $ls_3ltr_len;     # Length of the 3' LTR
    my $ls_ltr_homology; # Percent ID of LTRs
    my $ls_5dinuc;       # 5' dinucleotide sequence
    my $ls_3dinuc;       # 3' dinucleotide sequence
    my $ls_5tsr_seq;     # 5' Target site rep sequence
    my $ls_3tsr_seq;     # 3' Target site rep sequence
    my $ls_5flank_seq;   # 5' Flanking sequence
    my $ls_3flank_seq;   # 3' Flanking sequence
    my $ls_pbs_seq;      # Primer Binding Site Sequence
    my $ls_ppt_seq;      # Polypuring Tract sequence
    my $ls_5id_seq;
    my $ls_3id_seq;
    my $ls_5ltr_seq;     # Sequence of the 5' LTR
    my $ls_3ltr_seq;
    my $ls_full_retro_seq; 
    my $par_5ltr_len;      # Length of the 5' LTR as parsed
    my $par_3ltr_len;      # Length of the 3' LTR as parsed

    # BOOLEANS
    my $in_rt_frame_1 = 0;
    my $in_rt_frame_2 = 0;
    my $in_rt_frame_3 = 0;
    my $in_5id_seq = 0;
    my $in_3id_seq = 0;
    my $in_5ltr = 0;
    my $in_3ltr = 0;
    my $in_complete_seq = 0;
    my $in_aligned_ltrs = 0;

    # Coordinate values
    my $full_retro_start;
    my $full_retro_end;
    my $ltr5_start;
    my $ltr5_end;
    my $ltr3_start;
    my $ltr3_end;
    my $pbs_start;
    my $pbs_end;
    my $ppt_start;
    my $ppt_end;

    # GFF Coordinates values
    # These are from absolute start of the query sequence string
    # starting the index value at one
    my $gff_full_retro_start;
    my $gff_full_retro_end;
    my $gff_ltr5_start;
    my $gff_ltr5_end;
    my $gff_ltr3_start;
    my $gff_ltr3_end;
    my $gff_pbs_start;
    my $gff_pbs_end;
    my $gff_ppt_start;
    my $gff_ppt_end;
    my $gff_5tsr_start;
    my $gff_5tsr_end;
    my $gff_3tsr_start;
    my $gff_3tsr_end;

    # Coordinate substring tests
    my $ppt_from_full_retro;
    my $ppt_from_query_seq;
    my $pbs_from_full_retro;
    my $pbs_from_query_seq;
    my $tsd5_from_query_seq;
    my $tsd3_from_query_seq;

    #-----------------------------+
    # OPEN GFF OUTPUT FILE        |
    #-----------------------------+
    if ($gff_append) {
	open (GFFOUT, ">>$gff_out") ||
	    die "ERROR: Could not open gff output file:\b$gff_out\n";
    }
    else 
    {
	open (GFFOUT, ">$gff_out") ||
	    die "ERROR: Could not open gff output file:\n$gff_out\n";
    }
    
    #-----------------------------+
    # GET DATA FROM REPORT FILE   |
    #-----------------------------+
    open (REPIN, "<$report_in")
	|| die "ERROR: Can not open report file:\n$report_in\n";

    while (<REPIN>) {
	# Remove windows line endings
	# Since LTR_Struc work
	s/\r//g;
	chomp;
	if ($in_rt_frame_1) {

	}
	elsif (m/COMPLETE SEQUENCE OF PUTATIVE TRANSPOSON/) {
	    $in_3ltr = 0;
	    $in_complete_seq = 1;
	}
	elsif (m/^ALIGNED LTRS:/) {
	    $in_complete_seq = 0;
	    $in_aligned_ltrs = 1;
	}
	elsif ($in_rt_frame_2) {

	}
	elsif ($in_rt_frame_3) {

	}
	elsif ($in_5id_seq) {

	}
	elsif ($in_3id_seq) {

	}
	elsif ($in_complete_seq) {
	    $ls_full_retro_seq = $ls_full_retro_seq.$_;
	}
	elsif(/CUT-OFF SCORE:\s+(\d\.\d+)/){
	    $ls_score = $1;
	}
	elsif(m/TRANSPOSON IS IN (.*) ORIENTATION/) {
	    $ls_orientation = $1;
	    if ($ls_orientation =~ "NEGATIVE") {
		$ls_orientation = "-";
	    }
	    elsif ($ls_orientation =~ "POSITIVE") {
		$ls_orientation = "+";
	    }
	    else {
		# If return can not be parsed just use dot
		# this indicates unknown orientation if gff
		$ls_orientation = "."; 
	    }
	}
	#-----------------------------+
	# SEQUENCE DATA               |
	#-----------------------------+
	elsif ($in_5ltr) {
	    $ls_5ltr_seq = $ls_5ltr_seq.$_;
	    my $len_inline = length ($_);

	    # The following for debug
	    #print STDERR "\tLEN: $len_inline\n";

#	    if ($len_inline == 0) {
#           Changed the above to the following
#            .. this may cause some problems if single character for correct carriage
#               returns
	    if ($len_inline == 0 ) {
		$ltr5_blank_count++;
		if ($ltr5_blank_count == 2) {
		    # Set in_5ltr boolean to false
		    $in_5ltr = 0;
		    $in_3ltr = 1;
		}
	    }

	}
	elsif ($in_3ltr) {
	    $ls_3ltr_seq = $ls_3ltr_seq.$_;

	    
	}
	elsif (m/DINUCLEOTIDES: (.*)\/(.*)/) {
	    $ls_5dinuc = $1;
	    $ls_3dinuc = $2;
	}
	elsif (m/DIRECT REPEATS: (.*)\/(.*)/) {
	    $ls_5tsr_seq = $1;
	    $ls_3tsr_seq = $2;
	}
	elsif (m/PBS: (.*)/) {
	    $ls_pbs_seq = $1;
	}
	elsif (m/POLYPURINE TRACT: (.*)/){
	    $ls_ppt_seq = $1;
	}
	elsif (m/5\' FLANK: (.*)/) {
	    $ls_5flank_seq = $1;
	}
	elsif (m/3\' FLANK: (.*)/) {
	    $ls_3flank_seq = $1;
	}
	#-----------------------------+
	# OTHER DATA                  |
	#-----------------------------+
	elsif(m/OVERALL LENGTH OF TRANSPOSON:\s+(\d+)/){
	    $ls_retro_len = $1;
	}
	elsif(m/LENGTH OF LONGEST ORF:\s+(\d+)/){
	    $ls_longest_orf = $1;
	}
	elsif(m/LENGTH OF PUTATIVE 3\' LTR:\s+(\d+)/){
	    $ls_3ltr_len = $1;
	}
	elsif(m/LENGTH OF PUTATIVE 5\' LTR:\s+(\d+)/){
	    $ls_5ltr_len = $1;
	}
	elsif(m/LTR PAIR HOMOLOGY:\s+(\S+)%/){
	    $ls_ltr_homology = $1;
	}
	#-----------------------------+
	# SET BOOLEAN FLAGS           |
	#-----------------------------+
	elsif (m/LTRS:/) {
	    $in_5ltr = 1;
	}
    }

    close (REPIN);


    #-----------------------------+
    # GET COORDINATES             |
    #-----------------------------+
    $par_5ltr_len = length($ls_5ltr_seq);
    $par_3ltr_len = length($ls_3ltr_seq);

    print STDERR "5 LTR Len: $par_5ltr_len\t" if $verbose;
    print STDERR "$ls_5ltr_len\n" if $verbose;
    print STDERR "3 LTR Len: $par_3ltr_len\t" if $verbose;
    print STDERR "$ls_3ltr_len\n" if $verbose;


#    $full_retro_start = index($qry_seq,$ls_full_retro_seq) + 1;
    # Changed the above to the following to make the full retro start at zero
    $full_retro_start = 0;
    $full_retro_end = $full_retro_start + $ls_retro_len;
    
    # The following will have a problem on 100% identical LTRs
    # however, telling the search to start at the end of the
    # 5' LTR will solve this problem since the index function
    # will accept an offset at the third argument
    $ltr5_start = index ($ls_full_retro_seq, $ls_5ltr_seq) + 1;
    $ltr5_end = $ltr5_start + $ls_5ltr_len;
    $ltr3_start = index ($ls_full_retro_seq, $ls_3ltr_seq) + 1;
    $ltr3_end = $ltr3_start + $ls_3ltr_len;
    $pbs_start = index ($ls_full_retro_seq, $ls_pbs_seq) + 1 ;
    $pbs_end = $pbs_start + length($ls_pbs_seq);
    $ppt_start = index ($ls_full_retro_seq, $ls_ppt_seq) + 1;
    $ppt_end = $ppt_start + length($ls_ppt_seq);

    #-----------------------------+
    # GET EXTRACTED SEQS          |
    #-----------------------------+
    # This is to test if the coordinates I am returning matches the
    # observations that LTR_STRUC is reporting
    $ppt_from_full_retro = substr ($ls_full_retro_seq, $ppt_start - 1,
				   length($ls_ppt_seq) );
    $pbs_from_full_retro = substr ($ls_full_retro_seq, $pbs_start - 1,
				   length($ls_pbs_seq) );

    #-----------------------------+
    # GFF COORDINATES             |
    #-----------------------------+
    # Full retro (not includings tsds)
    $gff_full_retro_start = $full_retro_start;                      # OK
    $gff_full_retro_end = $gff_full_retro_start + $ls_retro_len - 1;# OK
    # Subunits (including the putative tsds)
    $gff_pbs_start = $pbs_start - 1 + $full_retro_start;            # OK
    $gff_pbs_end = $gff_pbs_start + length($ls_pbs_seq) - 1;        # OK
    $gff_ppt_start = $ppt_start - 1 + $full_retro_start;            # OK
    $gff_ppt_end = $gff_ppt_start + length($ls_ppt_seq) - 1;        # OK
    # LTRs
    $gff_ltr5_start = $gff_full_retro_start;                        # OK
    $gff_ltr5_end = $gff_ltr5_start + $ls_5ltr_len - 1;             # OK
    $gff_ltr3_end = $gff_full_retro_end;                            # OK
    $gff_ltr3_start = $gff_ltr3_end - $ls_3ltr_len + 1;             # OK
    # TSRs - Currently returning correct on positive strand
    $gff_5tsr_start = $gff_full_retro_start - length($ls_5tsr_seq); # OK
    $gff_5tsr_end = $gff_5tsr_start + length($ls_5tsr_seq) - 1;     # OK
    $gff_3tsr_start = $gff_full_retro_start + $ls_retro_len;        # OK
    $gff_3tsr_end = $gff_3tsr_start + length($ls_3tsr_seq) - 1;     # OK


    # Annotate the internal sequence features using local databases
    my @blast_dbs = ("/db/jlblab/pfam/gag_poaceae",
		     "/db/jlblab/pfam/zf_cchc_poaceae",
		     "/db/jlblab/pfam/rvp_all_poa",
		     "/db/jlblab/pfam/rve_poaceae",
		     "/db/jlblab/pfam/rvt_poaceae", # both rvt_1 and rvt_2
		     "/db/jlblab/pfam/rh_poaceae",
		     "/db/jlblab/pfam/chromo_viridiplant",
		     "/db/jlblab/pfam/dros_env",
		 );


    # SEQID
    my $name_root = $gff_seq_id;
    
    

    # PRINT GFF OUTPUT            |
    #-----------------------------+
    # Data type follows SONG
    # http://song.cvs.sourceforge.net/*checkout*/song/ontology/so.obo
#    my $name_root = "ltr_retro_model";
    my $gff_result_id = "ltr_struc_".$ls_id_num;
    print GFFOUT "$name_root\t".   # Seq ID
	"ltr_struc\t".             # Source
	"LTR_retrotransposon\t".        # Data type, has to be exon for APOLLO
	"$gff_full_retro_start\t". # Start
	"$gff_full_retro_end\t".   # End
	"$ls_score\t".             # Score
	"$ls_orientation\t".       # Strand
	".\t".                     # Frame
	"$gff_result_id\n";        # Retro Id
    print FEATSUM "$name_root\t";
    print FEATSUM "$gff_full_retro_start\t";
    print FEATSUM "$gff_full_retro_end\t";

    print GFFOUT "$name_root\t".   # Seq ID
	"ltr_struc\t".             # Source
	"five_prime_LTR\t".        # SO:0000425
	"$gff_ltr5_start\t".       # Start
	"$gff_ltr5_end\t".         # End
	"$ls_score\t".             # Score
	"$ls_orientation\t".       # Strand
	".\t".                     # Frame
	"$gff_result_id\n";        # Retro Id
    print FEATSUM "$gff_ltr5_start\t";
    print FEATSUM "$gff_ltr5_end\t";

    # DINUCLEOTIDES
    my $ltr_5_dn_start = substr($ls_5ltr_seq, 0, 2);
    my $ltr_5_dn_end = substr($ls_5ltr_seq, 
			      length($ls_5ltr_seq) - 2 , 2);
    print FEATSUM "$ltr_5_dn_start\t";
    print FEATSUM "$ltr_5_dn_end\t";


    print GFFOUT "$name_root\t".   # Seq ID
	"ltr_struc\t".             # Source
	"three_prime_LTR\t".       # SO:0000426  
	"$gff_ltr3_start\t".       # Start
	"$gff_ltr3_end\t".         # End
	"$ls_score\t".             # Score
	"$ls_orientation\t".       # Strand
	".\t".                     # Frame
	"$gff_result_id\n";        # Retro Id
    print FEATSUM "$gff_ltr3_start\t";
    print FEATSUM "$gff_ltr3_end\t";
    
    # SEQ HAS TG/CA
#    #print "\t\t3\'LTR: $ls_3ltr_seq\n";
    my $ltr_3_dn_start = substr($ls_3ltr_seq,0,2);
    my $ltr_3_dn_end = substr($ls_3ltr_seq,
			      length($ls_3ltr_seq) - 2, 2);
    print FEATSUM "$ltr_3_dn_start\t";
    print FEATSUM "$ltr_3_dn_end\t";

    print GFFOUT "$name_root\t".   # Seq ID
	"ltr_struc\t".             # Source
	"primer_binding_site\t".   # SO:0005850
	"$gff_pbs_start\t".        # Start
	"$gff_pbs_end\t".          # End
	"$ls_score\t".             # Score
	"$ls_orientation\t".       # Strand
	".\t".                     # Frame
	"$gff_result_id\n";        # Retro Id
    print FEATSUM "$gff_pbs_start\t";
    print FEATSUM "$gff_pbs_end\t";
    
    print GFFOUT "$name_root\t".   # Seq ID
	"ltr_struc\t".             # Source
	"RR_tract\t".              # SO:0000435 
	"$gff_ppt_start\t".        # Start
	"$gff_ppt_end\t".          # End
	"$ls_score\t".             # Score
	"$ls_orientation\t".       # Strand
	".\t".                     # Frame
	"$gff_result_id\n";        # Retro Id
    print FEATSUM "$gff_ppt_start\t";
    print FEATSUM "$gff_ppt_end\t";

    print GFFOUT "$name_root\t".     # Seq ID
	"ltr_struc\t".               # Source
	"target_site_duplication\t". # SO:0000434
	"$gff_5tsr_start\t".         # Start
	"$gff_5tsr_end\t".           # End
	"$ls_score\t".               # Score
	"$ls_orientation\t".         # Strand
	".\t".                       # Frame
	"$gff_result_id\n";          # Retro Id
    print FEATSUM "$gff_5tsr_start\t";
    print FEATSUM "$gff_5tsr_end\t";
    print FEATSUM "$ls_5tsr_seq\t";

    print GFFOUT "$name_root\t".     # Seq ID
	"ltr_struc\t".               # Source
	"target_site_duplication\t". # SO:0000434
	"$gff_3tsr_start\t".         # Start
	"$gff_3tsr_end\t".           # End
	"$ls_score\t".               # Score
	"$ls_orientation\t".         # Strand
	".\t".                       # Frame
	"$gff_result_id\n";          # Retro Id
    print FEATSUM "$gff_3tsr_start\t";
    print FEATSUM "$gff_3tsr_end\t";
    print FEATSUM "$ls_3tsr_seq\t";
    
    # Extract BLAST OUTPUT
    my $max_e_val = $e_val;
    my $seq_id = $name_root;
    my $seq_string = $ls_full_retro_seq;

    #///////////////////////
    foreach my $blastdb (@blast_dbs) {
	my $annotated_rpn = seq_annotate ($seq_id, $seq_string, 
					  $blastdb, $max_e_val );

	if ($annotated_rpn) {

	    print GFFOUT "$name_root\t".     # Seq ID
		"blastx\t".               # Source
		$annotated_rpn->{ feat_name }."\t".
		$annotated_rpn->{ feat_start }."\t".
		$annotated_rpn->{ feat_end }."\t".
		"$ls_score\t".               # Score
		"$ls_orientation\t".         # Strand
		".\t".                       # Frame
		"$gff_result_id\n";          # Retro Id

	    print FEATSUM $annotated_rpn->{ feat_start }."\t";
	    print FEATSUM $annotated_rpn->{ feat_end }."\t";
	    
	    if ($verbose) {
		print STDERR "$name_root\t".     # Seq ID
		    "blastx\t".               # Source
		    $annotated_rpn->{ feat_name }."\t".
		    $annotated_rpn->{ feat_start }."\t".
		    $annotated_rpn->{ feat_end }."\t".
		    "$ls_score\t".               # Score
		    "$ls_orientation\t".         # Strand
		    ".\t".                       # Frame
		    "$gff_result_id\n";          # Retro Id
	    }

	}
	else {
	    print FEATSUM "N\tN\t";
	}
    }


    # PRINT FLANKING SEQUENCES
    print FEATSUM "$ls_5flank_seq\t";
    print FEATSUM "$ls_3flank_seq\t";

    # PRINT ORIENTATION FROM LTR STRUC
    print FEATSUM "$ls_orientation";
    #my $ls_5flank_seq;   # 5' Flanking sequence
    #my $ls_3flank_seq;   # 3' Flanking sequence

    # END OF ALL FEATURES, PRINT NEWLINE TO FEATSUM
    print FEATSUM "\n";

    
    if ($verbose)  {
	print STDOUT "$name_root\t".   # Seq ID
	    "ltr_struc\t".             # Source
	    "LTR_retrotransposon\t".        # Data type, has to be exon for APOLLO
	    "$gff_full_retro_start\t". # Start
	    "$gff_full_retro_end\t".   # End
	    "$ls_score\t".             # Score
	    "$ls_orientation\t".       # Strand
	    ".\t".                     # Frame
	    "$gff_result_id\n";        # Retro Id
	
	print STDOUT "$name_root\t".   # Seq ID
	    "ltr_struc\t".             # Source
	    "primer_binding_site\t".   # SO:0005850
	    "$gff_pbs_start\t".        # Start
	    "$gff_pbs_end\t".          # End
	    "$ls_score\t".             # Score
	    "$ls_orientation\t".       # Strand
	    ".\t".                     # Frame
	    "$gff_result_id\n";        # Retro Id
	
	print STDOUT "$name_root\t".   # Seq ID
	    "ltr_struc\t".             # Source
	    "RR_tract\t".              # SO:0000435 
	    "$gff_ppt_start\t".        # Start
	    "$gff_ppt_end\t".          # End
	    "$ls_score\t".             # Score
	    "$ls_orientation\t".       # Strand
	    ".\t".                     # Frame
	    "$gff_result_id\n";        # Retro Id
	
	print STDOUT "$name_root\t".   # Seq ID
	    "ltr_struc\t".             # Source
	    "five_prime_LTR\t".        # SO:0000425
	    "$gff_ltr5_start\t".       # Start
	    "$gff_ltr5_end\t".         # End
	    "$ls_score\t".             # Score
	    "$ls_orientation\t".       # Strand
	    ".\t".                     # Frame
	    "$gff_result_id\n";        # Retro Id
	
	print STDOUT "$name_root\t".   # Seq ID
	    "ltr_struc\t".             # Source
	    "three_prime_LTR\t".       # SO:0000426  
	    "$gff_ltr3_start\t".       # Start
	    "$gff_ltr3_end\t".         # End
	    "$ls_score\t".             # Score
	    "$ls_orientation\t".       # Strand
	    ".\t".                     # Frame
	    "$gff_result_id\n";        # Retro Id
	
	print STDOUT "$name_root\t".     # Seq ID
	    "ltr_struc\t".               # Source
	    "target_site_duplication\t". # SO:0000434
	    "$gff_5tsr_start\t".         # Start
	    "$gff_5tsr_end\t".           # End
	    "$ls_score\t".               # Score
	    "$ls_orientation\t".         # Strand
	    ".\t".                       # Frame
	    "$gff_result_id\n";          # Retro Id
	
	print STDOUT "$name_root\t".     # Seq ID
	    "ltr_struc\t".               # Source
	    "target_site_duplication\t". # SO:0000434
	    "$gff_3tsr_start\t".         # Start
	    "$gff_3tsr_end\t".           # End
	    "$ls_score\t".               # Score
	    "$ls_orientation\t".         # Strand
	    ".\t".                       # Frame
	    "$gff_result_id\n";          # Retro Id
	
    }
    
    #-----------------------------+
    # PRINT SUMMARY OUTPUT        |
    #-----------------------------+
    if ($verbose) {
	print STDERR "\t\tScore: $ls_score\n";
	print STDERR "\t\tOrientation: $ls_orientation\n";
	print STDERR "\t\tRetro Len: $ls_retro_len\n";
	print STDERR "\t\tLongest Orf: $ls_longest_orf\n";
	print STDERR "\t\tPair Homology: $ls_ltr_homology\n";
	print STDERR "\t\t5LTR Len LST: $ls_5ltr_len\n";
	print STDERR "\t\t5LTR Lem PAR: $par_5ltr_len\n";
	print STDERR "\t\t3LTR Len LST: $ls_3ltr_len\n";
	print STDERR "\t\t3LTR Len PAR: $par_3ltr_len\n";
	
	print STDERR "\n\t\tEXTRACTED COORDINATES:\n";
	print STDERR "\t\t5\'LTR Start: $ltr5_start\n";
	print STDERR "\t\t5\'LTR End: $ltr5_end\n";
	print STDERR "\t\t3\'LTR Start: $ltr3_start\n";
	print STDERR "\t\t3\'LTR End: $ltr3_end\n";
	print STDERR "\t\tPBS Start: $pbs_start\n";
	print STDERR "\t\tPBS end : $pbs_end\n";
	print STDERR "\t\tPPT Start: $ppt_start\n";
	print STDERR "\t\tPPT End: $ppt_end\n";
	
	print STDERR "\t\tRetro Start: $full_retro_start\n";
	print STDERR "\t\tRetro End: $full_retro_end\n";
	
	# SEQUENCE DATA
	print STDERR "\n\t\tSEQUENCE DATA:\n";
	print STDERR "\t\t5Dinuc: $ls_5dinuc\n";
	print STDERR "\t\t3Dinuc: $ls_3dinuc\n";
	
	print STDERR "\n\t\t=====================\n";
	print STDERR "\t\t5\' TSD TEST:\n";
	print STDERR "\t\tLTR_STRUC: $ls_5tsr_seq\n";
	print STDERR "\t\t=====================\n";
	
	print STDERR "\n\t\t=====================\n";
	print STDERR "\t\t3\' TSD TEST:\n";
	print STDERR "\t\tLTR_STRUC: $ls_3tsr_seq\n";
	print STDERR "\t\t=====================\n";
	
	print STDERR "\n\t\t=====================\n";
	print STDERR "\t\tPBS STRING TEST:\n";
	print STDERR "\t\tLTR_STRUC:  $ls_pbs_seq\n";
	print STDERR "\t\tFROM RETRO: $pbs_from_full_retro\n";
	print STDERR "\t\t=====================\n";
	
	print STDERR "\n\t\t=====================\n";
	print STDERR "\t\tPPT STRING TEST:\n";
	print STDERR "\t\tLTR_STRUC:  $ls_ppt_seq\n";
	print STDERR "\t\tFROM RETRO: $ppt_from_full_retro\n";
	print STDERR "\t\t=====================\n";
	
	print STDERR "\n";
	print STDERR "\t\t5\'Flank: $ls_5flank_seq\n";
	print STDERR "\t\t3\'Flank: $ls_3flank_seq\n";
	#print "\t\t5\'LTR: $ls_5ltr_seq\n";
	#print "\t\t3\'LTR: $ls_3ltr_seq\n";
	#print "$ls_full_retro_seq\n";
    }



    # 05/16/2008 - Added the following back into the code
    

    #-----------------------------+
    # PRINT SEQ DATA              |
    #-----------------------------+
    if ($print_seq_data) {
	
	# This uses the global $outdir variable

	#-----------------------------+
	# 5' LTR                      |
	#-----------------------------+
	my $ltr5_seq_path = $outdir."ltr5_ltr_struc.fasta";
	open (LTR5OUT, ">>$ltr5_seq_path") ||
	    die "Can not open 5\'LTR sequence output file:\n$ltr5_seq_path\n";
	print LTR5OUT ">".$name_root."_".$gff_result_id.
	    "|five_prime_ltr\n";
	print LTR5OUT "$ls_5ltr_seq\n";
	close (LTR5OUT);
	
	#-----------------------------+
	# 3' LTR                      |
	#-----------------------------+
	my $ltr3_seq_path = $outdir."ltr3_ltr_struc.fasta";
	open (LTR3OUT, ">>$ltr3_seq_path") ||
	    die "Can not open 3\'LTR sequence output file:\n$ltr3_seq_path\n";
	print LTR3OUT ">".$name_root."_".$gff_result_id.
	    "|three_prime_ltr\n";
	print LTR3OUT "$ls_3ltr_seq\n";
	close (LTR3OUT);

	#-----------------------------+
	# PBS                         |
	#-----------------------------+
	my $pbs_seq_path = $outdir."pbs_ltr_struc.fasta";
	open (PBSOUT, ">>$pbs_seq_path") ||
	    die "Can not open PBS sequence output file:\n$pbs_seq_path\n";
	print PBSOUT ">".$name_root."_".$gff_result_id.
	    "|primer_binding_site\n";
	print PBSOUT "$ls_pbs_seq\n";
	close (PBSOUT);

	#-----------------------------+
	# PPT                         |
	#-----------------------------+
	my $ppt_seq_path = $outdir."ppt_ltr_struc.fasta";
	open (PBSOUT, ">>$ppt_seq_path") ||
	    die "Can not open PPT sequence output file:\n$ppt_seq_path\n";
	print PBSOUT ">".$name_root."_".$gff_result_id.
	    "|RR_tract\n";
	print PBSOUT "$ls_pbs_seq\n";
	close (PBSOUT);
	
	#-----------------------------+
	# FULL RETRO MODEL            |
	#-----------------------------+ 
	# NOTE THIS INCLUDES NESTED LTR RETROS
	# LTR_retrotransposon
	my $ltr_seq_out = $outdir."full_ltr_retro.fasta";
	open (LTROUT, ">>$ltr_seq_out") ||
	    die "Can not open full ltr retro sequence outfile\n";
	print LTROUT ">".$name_root."_".$gff_result_id.
	    "|LTR_retrotransposon\n";
	print LTROUT "$ls_full_retro_seq\n";
	close (LTROUT);
	    

    }

    close (GFFOUT);

}


sub seq_annotate {

    my $max_num_hits = 1;
    my @ans;              # Answer
    my $vals;             # Hash Reference to values
    my $feat_name;

    # annotate a seq from a datbase
    my ($seq_id, $seq_string, $blastdb, $max_e_val ) = @_;
    # dbh is the database handle where the data are to be store
    # seq_id is the id of the query sequence
    # seq_string is the
    # b should set the number of alignments returned
    my @bl_params = ('b'       => 1,
		     'e-value' => $max_e_val,
		     'program' => 'blastx', 
		     'database' => $blastdb);
    my $factory = Bio::Tools::Run::StandAloneBlast->new(@bl_params);
    
    my $qry_input = Bio::Seq->new(-id=> $seq_id,
				  -seq=>$seq_string );
    
    my $blast_report = $factory->blastall($qry_input);

    # This currently assumes a single query sequence was used
    my $hit_count = 0;
    
    while (my $blast_result = $blast_report->next_result()) {

	while (my $blast_hit = $blast_result->next_hit()) {
	    
	    while (my $blast_hsp = $blast_hit->next_hsp())
	    {

		if ($hit_count < $max_num_hits) {
		    my ($feat_start, $feat_end) = $blast_hsp->range('query');

		    # Print to STDERR IF VERBOSE
#		    print STDERR $feat_name."\n" if $verbose;    
#		    print $feat_start."\n" if $verbose;
#		    print $feat_end."\n" if $verbose;
#		    print STDERR
#			$blast_result->query_name."\n" if $verbose;
#		    print STDERR 
#			$blast_hsp->query_string."\n" if $verbose;
		    $hit_count++;
		    #my @range = $blast_hsp->range('query');

		    # Load values to the hash reference
		    #$vals->{'feat_name'} = $feat_name;
		    $vals->{'feat_name'} = $blast_result->database_name;
		    $vals->{'feat_start'} = $feat_start;
		    $vals->{'feat_end'} = $feat_end;
		    $vals->{'qry_id'} = $blast_result->query_name;
		    $vals->{'feat_seq'} = $blast_hsp->query_string;

		} # End of hit_count
	    }
	}
    } # End of while blast_result

    # Return the feature annotation
    if ($hit_count > 0) {
	return $vals;
    }
    else {
	return 0;
    }

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

# Old print_help
sub print_help {

    # Print requested help or exit.
    # Options are to just print the full 
    my ($opt) = @_;

    my $usage = "USAGE:\n". 
	"MyProg.pl -i InFile -o OutFile";
    my $args = "REQUIRED ARGUMENTS:\n".
	"  --infile       # Path to the input file\n".
	"  --outfile      # Path to the output file\n".
	"\n".
	"OPTIONS::\n".
	"  --version      # Show the program version\n".     
	"  --usage        # Show program usage\n".
	"  --help         # Show this help message\n".
	"  --man          # Open full program manual\n".
	"  --quiet        # Run program with minimal output\n";
	
    if ($opt =~ "full") {
	print "\n$usage\n\n";
	print "$args\n\n";
    }
    else {
	print "\n$usage\n\n";
    }
    
    exit;
}

=head1 NAME

cnv_ltrstruc2ann.pl - Convert LTR_STRUC report output files to gff

=head1 VERSION

This documentation refers to program version $Rev: 343 $

=head1 SYNOPSIS

=head2 Usage

    cnv_ltrstruc2ann.pl -i InDir -o OutDir -r LStrucOut

=head2 Required Arguments

    --indir         # Directory with the fasta files
    --outdir        # Directory for the base output dir
    --results       # Directory containing the LTR_STRUC results

=head1 DESCRIPTION

Given a directory containing the output from LTR_STRUC and a directory
containing the fasta files that the structural predictions were made
from,

=head1 REQUIRED ARGUMENTS

=over 2

=item -o,--outfile

Path of the output file. -- THIS DOES NOT GET USED

=item -g,--gff-out

Path to the gff format output file.

=item -f,--feat-sum

Path to the feature summary output file.

=item -r,--results

Path of the directory containing the results from LTR_STRUC. It is 
expected that these file names will end with rprt.txt.

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

=head1 DIAGNOSTICS

Error messages generated by this program and possible solutions are listed
below.

=over 2

=item ERROR: No fasta files were found in the input directory

The input directory does not contain fasta files in the expected format.
This could happen because you gave an incorrect path or because your sequence 
files do not have the expected *.fasta extension in the file name.

=item ERROR: Could not create the output directory

The output directory could not be created at the path you specified. 
This could be do to the fact that the directory that you are trying
to place your base directory in does not exist, or because you do not
have write permission to the directory you want to place your file in.

=back

=head1 CONFIGURATION AND ENVIRONMENT

The program cnv_ltrstruc2gff.pl does not currently require an 
external configuration file or make use of variables defined in the
user's environment.

=head1 DEPENDENCIES

=head2 Required Software

=over

=item * LTR_Struc

This program is designed to work with output files generated by
the LTR_Struc
program. This program is available for download from:
http://www.genetics.uga.edu/retrolab/data/LTR_Struc.html

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

=head1 SEE ALSO

The cnv_ltrstruc2gff.pl program is part of the DAWG-PAWS package of genome
annotation programs. See the DAWG-PAWS web page 
( http://dawgpaws.sourceforge.net/ )
or the Sourceforge project page 
( http://sourceforge.net/projects/dawgpaws ) 
for additional information about this package.

=head1 REFERENCE

Your use of the DAWGPAWS programs should reference the following manuscript:

JC Estill and JL Bennetzen. 2009. 
"The DAWGPAWS Pipeline for the Annotation of Genes and Transposable 
Elements in Plant Genomes." Plant Methods. 5:8.

The use of the cnv_ltrstruc2ann.pl program should also reference the 
LTR_STRUC manuscript:

EM McCarthy and JF McDonald. 2003. "LTR_STRUC: a novel search and 
identification program for LTR retrotransposons." Bioinformatics 19(3): 362-7.

=head1 LICENSE

GNU GENERAL PUBLIC LICENSE, VERSION 3

http://www.gnu.org/licenses/gpl.html

THIS SOFTWARE COMES AS IS, WITHOUT ANY EXPRESS OR IMPLIED
WARRANTY. USE AT YOUR OWN RISK.

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 09/25/2007

UPDATED: 04/15/2010

VERSION: $Rev: 343 $

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
# 09/25/2007
# - Basic input output started 
# 09/26/2007
# - ltrstruc2gff subfunction added
# 09/27/2007
# - Fixed coordinates in the gff output
# - Added code to copy ltr_struc related files to the
#   outdir.
# 09/28/2007
# - Working with gff files produced and wheat.tiers file to
#   to display the LTR-Retro subcomponent. Not really working out.
#   May have to define features of the ltr_retrospans
v#   and do annotation of the LTR retromodels in external process
# - Added code to generate fasta files for the following components
#    - Full Retro Sequence
#    - 5' LTR
#    - 3' LTR
#    - PBS
#    - PPT
# 12/12/2007
# - Changed print_help subfunction to extract help and 
#   usage messages from the POD documentation
# - Updated POD documentation
#
# 04/14/2008
# - Above documentation refers to cnv_ltrstruc2gff
# - Created cnv_ltrstruc2ann.pl from cnv_ltstruc2gff.
#   Unlike the earlier program this converts an ltrstruc
#   report to an annotation file mapped on the 
# - Added code to LTR_Struc Report parser to remove
#   windows line endings
# - Added seq_annotate subfunction
# - Can now extract feature locations from LTR_Struc and combine this
#   with BLASTX hits to an array of hard coded databases.
# - Adding feature summary file
#   This currently includes 39 feature elements
# 04/17/2008
# - Added e-value as an option
# 04/30/2008
# - Added the orientation to the feature summary file
#   This uses the orientation as reported by ltr_struc
# 05/16/2008
# - Added outdir back into code as the outdir option
#   at the commadn nline
