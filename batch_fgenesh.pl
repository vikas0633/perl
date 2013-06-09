#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_fgenesh.pl - Run Fgenesh in batch mode              |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 03/27/2010                                       |
# UPDATED: 03/29/2010                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Runs the Fgenesh gene annotation program in batch mode.  |
#                                                           |
# USAGE:                                                    |
#  batch_fgenesh.pl -i indir -o outdir -c fgen_config.cfg   |
#                                                           |
# VERSION: $Rev$                                            |
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
my $indir;                     # Dir for input
my $outdir;                    # Dir for output
my $config_file;               # Path to the config file
my @fg_params = ();            # Params for running Fgenesh
my $program = "FGenesh";       # Name used for source in GFF output
my $param;                     # Parameter name

# BOOLEANS
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_test = 0;                    # Run the program in test mode

# Path to the FGENESH program, assume in PATH otherwise
my $fg_path = $ENV{FGENESH_PATH} ||
    "fgenesh";

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED PARAMS
		    "i|indir=s"    => \$indir,
                    "o|outdir=s"   => \$outdir,
		    "c|config=s"   => \$config_file,
		    # OPTIONS
		    "program"      => \$program,
#		    "param"       => \$param,    # Set param
		    "q|quiet"      => \$quiet,
		    "fgenesh-bin", => \$fg_path, # FGenesh Path
		    "gff-ver=s"    => \$gff_ver,
		    "verbose"      => \$verbose,
		    # ADDITIONAL INFORMATION
		    "usage"        => \$show_usage,
		    "test"         => \$do_test,
		    "version"      => \$show_version,
		    "man"          => \$show_man,
		    "h|help"       => \$show_help,);

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
# CHECK REQUIRED ARGS         |
#-----------------------------+
if ( (!$indir) || (!$outdir) || (!$config_file) ) {
    print "\a";
    print STDERR "\n";
    print STDERR "ERROR: An input directory was not specified at the".
	" command line\n" if (!$indir);
    print STDERR "ERROR: An output directory was specified at the".
	" command line\n" if (!$outdir);
    print STDERR "ERROR: A config file was not specified at the".
	" command line\n" if (!$config_file);
    print_help ("usage", $0 );
}

#-----------------------------+
# CHECK FOR SLASH IN DIR      |
# VARIABLES                   |
#-----------------------------+
# If the indir does not end in a slash then append one
unless ($indir =~ /\/$/ ) {
    $indir = $indir."/";
}

unless ($outdir =~ /\/$/ ) {
    $outdir = $outdir."/";
}

#-----------------------------+
# Get the FASTA files from the|
# directory provided by the   |
# var $indir                  |
#-----------------------------+
opendir( DIR, $indir ) || 
    die "Can't open directory:\n$indir"; 
my @fasta_files = grep /\.fasta$|\.fa$/, readdir DIR ;
closedir( DIR );

my $count_files = @fasta_files;

#-----------------------------+
# SHOW ERROR IF NO FILES      |
# WERE FOUND IN THE INPUT DIR |
#-----------------------------+
if ($count_files == 0) {
    print STDERR "\a";
    print STDERR "\nERROR: No fasta files were found in the input directory\n".
	"$indir\n".
	"Fasta files must have the fasta or fa extension.\n\n";
    exit;
}

#-----------------------------+
# LOAD THE CONFIG FILE        |
#-----------------------------+
my $i=0;
my $config_line_num=0;

open (CONFIG, "<$config_file") ||
    die "ERROR Can not open the config file:\n $config_file";

while (<CONFIG>) {
    chomp;
    $config_line_num++;
    unless (m/^\#/) {
       	my @in_line = split;           # Implicit split of $_ by whitespace
	my $num_in_line = @in_line; 
	
# HOW TO COUNT EXPECTED INPUT 
	if ( $num_in_line == 2 ||
	     $num_in_line == 3) { 
	    # May want name boolean then path
	    $fg_params[$i][0] = $in_line[0] || "NULL";  # Name
	    $fg_params[$i][1] = $in_line[1] || "NULL";  # Path
	    # If using HML_seq this path must be included in 3rd column
	    $fg_params[$i][3] = $in_line[3] || "NULL";  # Other OPTIONS
	    $i++;
	} # End of if $num_in_line is 10
	else {
	    print "\a";
	    print STDERR "WARNING: Unexpected number of options in config".
		" file line $config_line_num\n$config_file\n";
	}

   } # End of unless comment line
} # End of while CONFIG file
close CONFIG;


# Number of parameter sets specified in the config file
my $num_par_sets = $i;

if ($num_par_sets == 0) {
    print "\a";
    print STDERR "ERROR: No parameter sets were found in the config file:\n".
	"$config_file\n";
    exit;
}

my $num_proc_total = $count_files * $num_par_sets;


print STDERR "NUMBER OF FILES TO PROCESS: $count_files\n" if $verbose;

#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+
# For each fasts file in the input directory

# For each fasta
# TO DO : Consider pulling from fasta header for input name
my $file_num = 0;
my $proc_num = 0;
for my $ind_file (@fasta_files) {
    
    $file_num++;
    my $name_root;
    
    #-----------------------------+
    # Get the root name of the    |
    # file to predict                |
    #-----------------------------+
    if ($ind_file =~ m/(.*)\.masked\.fasta$/) {
	# file ends in .masked.fasta
	$name_root = "$1";
    }
    elsif ($ind_file =~ m/(.*)\.fasta$/ ) {	    
	# file ends in .fasta
	$name_root = "$1";
    }  
    elsif ($ind_file =~ m/(.*)\.fa$/ ) {	    
	# file ends in .fa
	$name_root = "$1";
    } 
    else {
	$name_root = $ind_file;
    }
    
    #-----------------------------+
    # Create parent dir if it     |
    # does not already exist      |
    #-----------------------------+
    my $dir_parent = $outdir.$name_root."/";
    unless (-e $dir_parent) {
	print STDERR "creating dir: $dir_parent\n" if $verbose;
	mkdir $dir_parent ||
	    die "Could not creat the output dir:\n$dir_parent\n";
    }
    
    
    #-----------------------------+
    # Create the dir to hold the  |
    # FGENESH output              |
    #-----------------------------+
    my $dir_fgenesh_out = $outdir.$name_root."/fgenesh/";
    unless (-e $dir_fgenesh_out ) {
	print STDERR "Creating output dir\n: $dir_fgenesh_out\n" 
	    if $verbose;
	mkdir $dir_fgenesh_out ||
	    die "Could not create the output directory:\n$dir_fgenesh_out";
    }

    #-----------------------------+
    # CREATE THE DIR TO HOLD      |
    # THE GFF OUTPUT              |
    #-----------------------------+
    my $dir_gff_out = $outdir.$name_root."/gff/";
    unless (-e $dir_gff_out) {
	print STDERR "Creating output dir\n:$dir_gff_out"
	    if $verbose;
	mkdir $dir_gff_out ||
	    die "Could not create the gff output directory:\n$dir_gff_out";
    }

    #-----------------------------+
    # FOR EACH PARAM SET  IN THE  |
    # CONFIG FILE                 |
    #-----------------------------+
    my $gff_count=0;               # Count of the number of gff files
#    for my $ind_param (@fg_params) {
    for ($i=0; $i<$num_par_sets; $i++) {	
	$proc_num++;
	my $fg_param_name = $fg_params[$i][0];
	my $fg_matrix = $fg_params[$i][1];
	my $fg_opt = $fg_params[$i][2] || "NULL";

	# FILE PATHS
	my $fg_infile = $indir.$ind_file;
	my $fg_resout = $dir_fgenesh_out.$name_root.
	    "_fg_".$fg_param_name.".txt";
	# GFF Format FGENESH results
	my $fg_gffout = $dir_fgenesh_out.$name_root.
	    "_fg_".$fg_param_name.".gff";

	#-----------------------------+
	# RUN FGENESH                 |
	#-----------------------------+
	my $fg_cmd = $fg_path." ".$fg_matrix." ".$fg_infile;
	unless ($fg_opt =~ "NULL") {
	    $fg_cmd = $fg_cmd." ".$fg_opt;
	}
	$fg_cmd = $fg_cmd." > ".$fg_resout;

	# CONVERT RESULT TO GFF
	print STDERR "CMD: ".$fg_cmd."\n" if $verbose; 


	system ($fg_cmd) || 
	    die "ERROR: Could not run command: $fg_cmd \n";

	#-----------------------------+
	# CONVERT OUTPUT TO GFF       |
	#-----------------------------+
	my $append = 0;
	if (-e $fg_resout) {
	    fgenesh2gff ($program, $fg_resout, $fg_gffout, 
			 $name_root, $fg_param_name, $append);
	}

    }
    

} # End of for each fasta file in the input directory

#-----------------------------+
# RUN FGENESH PROGRAM         |
#-----------------------------+
# Running the fgenesh program is like
# fgenesh MATRIX_NAME infile_seq.fasta > result.txt
#
# Options are placed after the infile_seq option
#
# fgenesh MATRIX_NAME infile_seq.fasta -pmrna > result.txt
# will print the mRNA sequence of the predicted gene
# OPTIONS ALSO:
# Program Par_file Seq_File <Hml_Seq> <-options>
#
# Options:
#   -GC:xxx            -  Use potential GC donor splice sites with score above
#                         xxx. xxx is floating point value corresponded
#                          the number
#                         of matched letters in site consensus (xxx >= 0).
#   -GC                -  Use all potential GC donor splice sites.
#   -p1:xxx            -  Get sequence from position xxx.
#   -p2:xxx            -  Get sequence to position xxx.
#   -c                 -  Use condensed sequence.
#   -min_hml:xxx       -  Minimal considered homology (FGENESH+ and 
#                         FGENESH_C only). Default is 90%
#   -exon_table:file   -  File with table of exons.
#   -exon_bonus:xxx    -  Add bonus xxx for all exons from table.
#   -hml_bonus:xxx     -  Addition multiplier for homology (default - 1.0)
#                         (FGENESH_C only).
#   -send              -  Soft homology termination (FGENESH_C only).
#   -full_start        -  Homologous sequences have a head (FGENESH_C only).
#   -full_end t        -  Homologous sequences have a tail (FGENESH_C only).
#   -ipen:xxx          -  Penalties for all internal exons without homology 
#                         (FGENESH_C only).
#   -pmrna             -  Print mRNA sequences for predicted genes.
#   -pexons            -  Print exons sequences for predicted genes.
#   -t:table           -  Use translation table.
#   -st:table          -  Print selested translation table.
#                         table values are - 
#              1 - Standard. (Default)
#              2 - Vertebrate Mitochondrial.
#              3 - Yeast Mitochondrial.
#              4 - Mold Mitochondria, Protozoan Mitochondrial, 
#                  Colenterate Mitochondrial,
#                  Mycoplasma, Spiroplasma.
#              5 - Invertebrate Mitochondrial.
#              6 - Ciliate Nuclear, Dasycladacean Nuclear, Hexamita Nuclear.
#              9 - Echinoderm Nuclear.
#             10 - Euplotid Nuclear.
#             11 - Bacterial.
#             12 - Alternative Yeast Nuclear.
#             13 - Ascidian Mitochondrial.
#             14 - Flatworm Mitochondrial.
#             15 - Blepharisma Macronuclear.


exit 0;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+
sub fgenesh2gff {
    # fgnesh_in  - path to the fgenesh program
    # gff_out    - path to the gff output file
    # seq_id     - id of the source sequence
    # src_suffix - parameter id for fgenesh run

    my ($source, $fgenesh_in, $gff_out, $seq_id, $src_suffix, $do_append ) = @_;

    my $attribute;

    #-----------------------------+
    # OPEN THE FGENESH INFILE     |
    #-----------------------------+
    my $fgenesh_result;
    if ($fgenesh_in) {
	$fgenesh_result = Bio::Tools::Fgenesh->new(-file => $fgenesh_in);
    }
    else {
	print STDERR "Expecting input from STDIN\n";
	$fgenesh_result = Bio::Tools::Fgenesh->new( -fh  => \*STDIN );
    }

    #-----------------------------+
    # OPEN THE GFF OUTFILE        |
    #-----------------------------+
     # Default to STDOUT if no argument given
    if ($gff_out) {
	if ($do_append) {
	    open (GFFOUT, ">>$gff_out") ||
		die "ERROR: Can not open gff outfile:\n $gff_out\n";
	}
	else {
	    open (GFFOUT,">$gff_out") ||
		die "ERROR: Can not open gff outfile:\n $gff_out\n";
	    if ($gff_ver =~ "GFF3") {
		print GFFOUT "##gff-version 3\n";
	    }
	    
	}
    } 
    else {
	open (GFFOUT, ">&STDOUT") ||
	    die "Can not print to STDOUT\n";
	if ($gff_ver =~ "GFF3") {
	    print GFFOUT "##gff-version 3\n";
	}

    }


    #-----------------------------+
    # SET PROGRAM SOURCE          |
    #-----------------------------+
    unless ($source) {
	$source = "fgenesh";
    }
    if ($src_suffix) {
	$source = $source.":".$src_suffix;
    }

    my $gene_num = 0;

    while (my $gene = $fgenesh_result->next_prediction()) {

	$gene_num++;     
	my $gene_id = sprintf("%05d", $gene_num); # Pad the gene number
	my $gene_name = "Fgenesh_gene_".$gene_id."\n";
	$gene_id = "gene".$gene_id;

	#-----------------------------+
	# SET SEQUENCE ID             |
	#-----------------------------+
	unless ($seq_id) {
	    if ($gene->seq_id()) {
		$seq_id = $gene->seq_id();
	    }
	    else {
		$seq_id = "seq";
	    }
	}

	# NOTE:
	# $gene is an instance of Bio::Tools::Prediction::Gene, which inherits
	# off Bio::SeqFeature::Gene::Transcript.
	#
	# $gene->exons() returns an array of 
	# Bio::Tools::Prediction::Exon objects
	# all exons:

	#-----------------------------+
	# IF GFF3 PRINT PARENT GENE   |
	# AND TRANSCRIPT INFORMATION  |
	#-----------------------------+

	if ($gff_ver =~ "GFF3") {

	    #-----------------------------+
	    # GENE                        |
            #-----------------------------+

	    $attribute = "ID=".$source."_".$gene_id;
	    #$gene_id;
#	    my @exons_ordered = $gene->exons_ordered();

	    # Get gene start and end
	    my $gene_start = $gene->start();
	    my $gene_end = $gene->end();
	    if ($gene_start > $gene_end) {
		$gene_end =  $gene->start();
		$gene_start = $gene->end();
	    }

	    # Get gene score
	    my $gene_score;
	    if ($gene->score()) {
		$gene_score = $gene->score();
	    }
	    else {
		$gene_score = ".";
	    }

	    # Get gene strand
	    my $gene_strand = $gene->strand()."\t";
	    if ($gene_strand =~ "-1") {
		$gene_strand = "-";
	    }
	    else {
		$gene_strand = "+";
	    }

	    print GFFOUT $seq_id."\t".     # Seqname
		$source."\t".              # Source
		"gene\t".                  #feature
		"$gene_start\t".           # start
		"$gene_end\t".             # end
		"$gene_score\t".           # score
		"$gene_strand\t".          # strand
		".\t".                     # frame
		$attribute.                # attribute
		"\n";

#	    #-----------------------------+
#	    # TRANSCRIPT
#	    #-----------------------------+
#	    #my @transcripts = $gene->transcripts();
#	    my @promoters = $gene->promoters();
#	    my @utr_5 = $gene->utrs('utr5prime');
#	    my @utr_3 = $gene->utrs('utr3prime');
#
#	    # Poly a is returning a single base
#	    my $poly_a_site = $gene->poly_A_site;
#	    my $poly_a_start = $poly_a_site->start();
#	    my $poly_a_end = $poly_a_site->end();
#	    if ($poly_a_site) {
#		print STDERR "polya".$poly_a_start."-".$poly_a_end."\n";
#	    }
#	    foreach my $ind_utr_5 (@utr_5) {
#		print STDERR "YUP\n\n";
#	    }
#	    foreach my $promoter (@promoters) {
#		print STDERR "Prom\t".
#		    $promoter->start()."\t".
#		    $promoter->end()."\t".
#		    "\n;"
#	    }
#	    foreach my $poly_a_site (@polya_sites) {
#		print STDERR "polya\n";
#	    }
#	    print STDERR ;

	}


	#-----------------------------+
	# EXONS                       |
	#-----------------------------+
	# Exons will be directly assigned to the gene as the parent
	my @exon_arr = $gene->exons();
	
	my $exon_num = 0;
	foreach my $ind_exon (@exon_arr) {

	    $exon_num++;
	    my $exon_id = sprintf("%05d", $exon_num); # Pad the exon number
	    $exon_id = "exon".$exon_id;
#	    $exon_id = 
	    #print STDERR $ind_exon;
#	    if ($ind_exon->is_coding()) {
#		print STDERR "coding\t";
#	    }
#	    else {
#	    }

	    #-----------------------------+
	    # FORMAT STRAND               |
	    #-----------------------------+
	    my $strand = $ind_exon->strand()."\t";
	    if ($strand =~ "-1") {
		$strand = "-";
	    }
	    else {
		$strand = "+";
	    }

	    #-----------------------------+
	    # GET START AND END           |
	    #-----------------------------+
	    my $start = $ind_exon->start();
	    my $end = $ind_exon->end();
	    if ($start > $end) {
		$end =  $ind_exon->start();
		$start = $ind_exon->end();
	    }

	    #-----------------------------+
	    # START GFF3 WORK HERE        |
	    #-----------------------------+
	    my $attribute;
	    my $feature;

	    if ($gff_ver =~ "GFF3") {
		$attribute = "ID=".$source."_".$gene_id."_".$exon_id.
		    "\;Parent=".$source."_".$gene_id;
		$feature = "exon";
	    }
	    else {
		$feature = "exon";
		$attribute = "gene_".$gene_num;
	    }

	    print GFFOUT $seq_id."\t".     # Seqname
		$source."\t".              # Source
		$feature."\t".             #feature
		$start."\t".               # start
		$end."\t".                 # end
		$ind_exon->score()."\t".   # score
		$strand."\t".              # strand
		".\t".                     # frame
		$attribute.                # attribute
		"\n";


	    # The following does work
	    #print GFFOUT $ind_exon->primary_tag()."\t";


	    #///////////////////////////////
	    # The following do not work
	    #///////////////////////////////
	    #print GFFOUT $gene->cds()."\n";
	    #print GFFOUT $ind_exon->significance()."\t";
	    #print $ind_exon->predicted_cds();
	    #print GFFOUT $ind_exon->coding_signal_score()."\t";
	    #print GFFOUT $ind_exon->seq_id()."\t";
	    #print $ind_exon->significance()."\n";
	    # Get the CDS of the sequence
	    #print STDERR $ind_exon->cds()."\n";

	}
	
#       # initial exons only
#       @init_exons = $gene->exons('Initial');
#       # internal exons only
#       @intrl_exons = $gene->exons('Internal');
#       # terminal exons only
#       @term_exons = $gene->exons('Terminal');
#       # singleton exons: 
#       ($single_exon) = $gene->exons();


   }

    # CLOSE FGENESH
    $fgenesh_result->close();
    

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

batch_fgenesh.pl - Run Fgenesh in batch mode

=head1 VERSION

This documentation refers to program version $Rev$

=head1 SYNOPSIS

  USAGE:
    batch_fgenesh.pl -i indir -o outdir -c fgen_config.cfg

    --indir        # Path to the input directory
    --outdir        # Path to the output directory

=head1 DESCRIPTION

Runs the Fgenesh gene annotation program in batch mode.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path of the directory containing the sequences to process.

=item -o,--outdir

Path of the directory to place the program output.

=item -c, --config

Path to a config file. This is a tab delimited text file
indicating the required information for each of the databases to blast
against. Lines beginning with # are ignored.

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

=item --verbose

Run the program with maximum output.

=item -q,--quiet

Run the program with minimal output.

=item --test

Run the program without doing the system commands. This will
test for the existence of input files.

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

=head2 Configuration File

The location of the configuration file is indicated by the --config option
at the command line.
This is a tab delimited text file
indicating required information for each of the databases to blast
against. Lines beginning with # are ignored, and data are in six 
columns as shown below:

=over 2

=item Col 1. Blast program to use [ tblastx | blastn | blastx ]

The blastall program to use. DAWG-PAWS will support blastn,
tblastx, and blastx format.

=item Col 2. Extension to add to blast output file. (ie. bln )

This is the suffix which will be added to the end of your blast
output file. You can use this option to set different extensions
for different types of blast. For example *.bln for blastn
output and *.blx for blastx output.

=back

An example config file:

 #-----------------------------+
 # BLASTN: TIGR GIs            |
 #-----------------------------+
 blastn	bln	8	1e-5	TaGI_10	-a 2 -U
 blastn	bln	8	1e-5	AtGI_13	-a 2 -U
 blastn	bln	8	1e-5	ZmGI_17	-a 2 -U
 #-----------------------------+
 # TBLASTX: TIGR GIs           |
 #-----------------------------+
 tblastx	blx	8	1e-5	TaGI_10	-a 2 -U
 tblastx	blx	8	1e-5	AtGI_13	-a 2 -U
 tblastx	blx	8	1e-5	ZmGI_17	-a 2 -U

=head2 Environment

This program does not make use of variables in the user environment.

=head1 DEPENDENCIES

=head2 Required Software

=over

=item * Software Name

Any required software will be listed here.

=back

=head2 Required Perl Modules

=over

=item * Getopt::Long

This module is required to accept options at the command line.

=back

=head1 BUGS AND LIMITATIONS

Any known bugs and limitations will be listed here.

=head2 Bugs

=over 2

=item * No bugs currently known 

If you find a bug with this software, file a bug report on the DAWG-PAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=204962

=back

=head2 Limitations

=over

=item * Known Limitation

If this program has known limitations they will be listed here.

=back

=head1 SEE ALSO

The program is part of the DAWG-PAWS package of genome
annotation programs. See the DAWG-PAWS web page 
( http://dawgpaws.sourceforge.net/ )
or the Sourceforge project page 
( http://sourceforge.net/projects/dawgpaws ) 
for additional information about this package.

=head1 REFERENCE

Please refer to the DAWGPAWS manuscript in Plant Methods when describing
your use of this program:

JC Estill and JL Bennetzen. 2009. 
"The DAWGPAWS Pipeline for the Annotation of Genes and Transposable 
Elements in Plant Genomes." Plant Methods. 5:8.

as well as the Fgenesh program description:

Solovyev V.V. (2001) Statistical approaches in Eukaryotic gene
prediction. In Handbook of Statistical genetics (eds. Balding D. et al.),
John Wiley & Sons, Ltd., p. 83-127.

=head1 LICENSE

GNU General Public License, Version 3

L<http://www.gnu.org/licenses/gpl.html>

THIS SOFTWARE COMES AS IS, WITHOUT ANY EXPRESS OR IMPLIED
WARRANTY. USE AT YOUR OWN RISK.

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 03/27/2010

UPDATED: 03/29/2010

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 03/27/2010
# - Program started
