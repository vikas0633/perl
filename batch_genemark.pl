#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_genemark.pl - Run Genemark gene prediction program  |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill at gmail.com                         |
# STARTED: 11/09/2007                                       |
# UPDATED: 03/16/2010                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Run the genmark gene prediction program in batch mode.   |
#  Runs genmark as well as converts output to gff format.   |
#  Requires a config file to specify libraries to use.      |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use File::Copy;
use Getopt::Long;
use Bio::Tools::Genemark; 
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
my @config_opts;               # Configuration file options
my $file_config;               # Path to the configuration file
my $logfile;                   # Path to a logfile to log error info
my $indir;                     # Directory containing the seq files to process
my $outdir;                    # Directory to hold the output
my $msg;                       # Message printed to the log file

my $search_name;               # Name searched for in grep command
my $bac_out_dir;               # Dir for each sequnce being masked
my $name_root;                 # Root name to be used for output etc

# GENMARK PATHS
# genmakr_dir is the dir that contains the genmark binaries
my $genmark_dir = $ENV{GM_BIN_DIR};
# lib_dir is the dir that contains the genmark matrix librarires
my $lib_dir = $ENV{GM_LIB_DIR};

# The name used as the source program
my $src_prog = "GeneMarkHMM";  


# BOOLEANS
my $show_help = 0;             # Show program help
my $show_version = 0;          # Show program version
my $show_man = 0;              # Show program manual page using peldoc
my $show_usage = 0;            # Show program usage command             
my $quiet = 0;                 # Boolean for reduced output to STOUT
my $apollo = 0;                # Path to apollo and apollo variables
my $test = 0;
my $verbose = 0;

# COUNTERS
my $num_proc = 1;              # Number of processes

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(
		    # Required
		    "i|indir=s"      => \$indir,
                    "o|outdir=s"     => \$outdir,
		    "c|config=s"     => \$file_config,
		    # Optional strings
		    "genemark-dir=s" => \$genmark_dir,
		    "lib-dir=s"      => \$lib_dir,
		    "program=s"      => \$src_prog, 
		    #"logfile=s"      => \$logfile,
		    # Booleans
		    #"apollo"         => \$apollo,
		    "verbose"        => \$verbose,
		    "test"           => \$test,
		    "usage"          => \$show_usage,
		    "version"        => \$show_version,
		    "man"            => \$show_man,
		    "h|help"         => \$show_help,
		    "q|quiet"        => \$quiet,);

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

my $proc_num = 0;

#//////////////////////
my $file_num_max = 4;
my $file_num = 0;
#\\\\\\\\\\\\\\\\\\\\\\

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
    print "\nbatch_genmark.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}

#-----------------------------+
# CHECK REQUIRED ARGS         |
#-----------------------------+
if ( (!$indir) || (!$outdir) || (!$file_config) || (!$lib_dir) || 
     (!$genmark_dir) ) {
    print "\a";
    print STDERR "\n";
    print STDERR "ERROR: A paramter file directory was not specified".
	" at the command line" if (!$lib_dir);
    print STDERR "ERROR: A genmark binary directory  was not specified".
	" at the command line" if (!$genmark_dir);
    print STDERR "ERROR: An input directory was not specified at the".
	" command line\n" if (!$indir);
    print STDERR "ERROR: An output directory was not specified at the".
	" command line\n" if (!$outdir);
    print STDERR "ERROR: An configuration file was not specified at the".
	" command line\n" if (!$file_config);
    print_help("full");
}


#-----------------------------+
# CHECK FOR SLASH IN DIR      |
# VARIABLES                   |
#-----------------------------+
# If the indir does not end in a slash then append one
# TO DO: Allow for backslash
unless ($indir =~ /\/$/ ) {
    $indir = $indir."/";
}

unless ($outdir =~ /\/$/ ) {
    $outdir = $outdir."/";
}

#-----------------------------+
# GET INFO FROM CONFIG FILE   |
#-----------------------------+
# Load to the 2d dbs array
open (CONFIG, $file_config) ||
    die "Can't open the config file:\n$file_config";
$i = 0;
my $line_num = 0;
while (<CONFIG>) {
    $line_num ++;
    
    # ignore comment lines
    next if m/^\#/;       

    chomp;
    my @tmpary = split (/\t/);
    my $count_tmp = @tmpary;

    if ($count_tmp == 4) {
	
	#-----------------------------+
	# LOAD CONFIG INFO            |
	#-----------------------------+
	$config_opts[$i]{param_name} = $tmpary[0];
	$config_opts[$i]{param_short} = $tmpary[1];
	$config_opts[$i]{gmhmme_ver} = $tmpary[2];
	$config_opts[$i]{param_file} = $tmpary[3];

	#-----------------------------+
	# SET PROGRAM NUMBER TO       |
	# PROGRAM NAME OF THE GMHMME  |
	# BINARY PROGRAM              |
	#-----------------------------+
	if ($config_opts[$i]{gmhmme_ver} =~ 3 ) {
	    $config_opts[$i]{gmhmme_ver} = "gmhmme3";
	}
	elsif ($config_opts[$i]{gmhmme_ver} =~ 2 ) {
	    $config_opts[$i]{gmhmme_ver} = "gmhmme2";
	}
	else {
	    my $err_msg = "ERROR: The GMHMMER version option in line ".
		$line_num." is set to ".$config_opts[$i]{gmhmme_ver}.
		"Valid options are 2 or 3.";
	    die $err_msg;
	}
	
	# TEST THAT THE PARAM FILE EXISTS
	my $param_path = $lib_dir.$config_opts[$i]{param_file};
	unless (-e $param_path) {
	    my $err_msg = "ERROR: The paramter file could not be located at:\n".
		"$param_path\n";
	    die $err_msg;
	}

	$i++;
	
    } # End of if count_tmp = 4
    else {
	print STDERR "WARNING: Config file line number $line_num\n";
	print STDERR "         Only $line_num variables were found\n" 
    }
    
    
} # End of while CONFIG

my $num_configs = @config_opts;
if ($num_configs == 0) {
    my $err_msg = "ERROR: No configuration options were present in the"; 
    print STDERR "\a";
    print STDERR $err_msg;
    print STDERR 
    die;
}
else {
    print STDERR "Processing $num_configs parameter files\n" if $verbose;
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
my $num_files = @fasta_files;

#-----------------------------+
# SHOW ERROR IF NO FILES      |
# WERE FOUND IN THE INPUT DIR |
#-----------------------------+
if ($num_files == 0) {
    print STDERR "\a";
    print STDERR "\nERROR: No fasta files were found in the input directory\n".
	"$indir\n".
	"Fasta files must have the fasta or fa extension.\n\n";
    exit;
}

#-----------------------------+
# CREATE THE OUT DIR          |
# IF IT DOES NOT EXIST        |
#-----------------------------+
unless (-e $outdir) {
    print STDERR "Creating output dir ...\n" unless $quiet;
    mkdir $outdir ||
	die "Could not create the output directory:\n$outdir";
}


for my $ind_file (@fasta_files) {

    $proc_num++;
    $file_num++;
    #if ($file_num == $file_num_max){exit;}

    #-----------------------------+
    # GET THE ROOT NAME OF THE    |
    # FASTA FILE                  |
    #-----------------------------+
    if ($ind_file =~ m/(.*)\.hard\.fasta$/) {
	# file ends in .hard.fasta
	# This is hard masked fasta files
	$name_root = "$1";
    }
    elsif ($ind_file =~ m/(.*)\.masked\.fasta$/) {
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
    # CREATE ROOT NAME DIR        |
    #-----------------------------+
    my $name_root_dir = $outdir.$name_root;
    unless (-e $name_root_dir) {
	mkdir $name_root_dir ||
	    die "Could not create dir:\n$name_root_dir\n"
    }

    #-----------------------------+
    # CREATE GENMARK OUTDIR       |
    #-----------------------------+
    # Dir to hold gene prediction output from local software
    my $genmark_out_dir = $outdir.$name_root."/genemark/";
    unless (-e $genmark_out_dir) {
	mkdir $genmark_out_dir ||
	    die "Could not create genemark out dir:\n$genmark_out_dir\n";
    }

    #-----------------------------+
    # CREATE GFF OUTDIR           |
    #-----------------------------+
    # Dir to hold gene prediction output from local software
    my $gff_out_dir = $outdir.$name_root."/gff/";
    #print STDERR "$gff_out_dir\n";
    unless (-e $gff_out_dir) {
	mkdir $gff_out_dir ||
	    die "Could not create gff out dir:\n$gff_out_dir\n";
    }

    print STDERR "\n===========================================\n" if $verbose;
    print STDERR " Running GeneMark for $name_root\n" if $verbose;
    print STDERR " File $file_num of $num_files\n" if $verbose;
    print STDERR "===========================================\n" if $verbose;
    
    my $infile_path = $indir.$ind_file;
    my $out_dir = $genmark_out_dir.$name_root.".genmark.out";
    my $gff_path = $genmark_out_dir.$name_root.".genmark.gff";

    #-----------------------------+
    # RUNN FOR EACH PARAMETER     |
    # SET IN THE CONFIG FILE      |
    #-----------------------------+
    for my $config (@config_opts) {
	print STDERR " Processing: ".$name_root." - ".
	    $config->{param_name}."\n" if $verbose;
	
	my $gff_out = $gff_out_dir.$name_root."_genemark_".
	    $config->{param_short}.".gff";

	my $gm_out = $genmark_out_dir.$name_root."_genemark_".
	    $config->{param_short}.".out";

	my $gm_cmd;
	if ($config->{gmhmme_ver} =~ "gmhmme2") {
	    $gm_cmd = $genmark_dir."gmhmme2 -m ";
	}
	elsif ( $config->{gmhmme_ver} =~ "gmhmme3" ) {
	    $gm_cmd = $genmark_dir."gmhmme3 -p -m";
	}
	else {
	    die "ERROR: Unknown gmhmme command:".$config->{gmhmme_ver};
	}
	
	$gm_cmd = $gm_cmd.$lib_dir.$config->{param_file}.
	    " -o ".$gm_out.
	    " ".$infile_path;

	print STDERR "\t$gm_cmd\n" if $verbose;
	system($gm_cmd) unless $test;
	
	# Convert result to gff format
	if (-e $gm_out) {
#	    # OLD
#	    genemark_to_gff($gm_out, $gff_out, 
#			    $name_root, $conig->{param_name});
	    # NEW FUNCTION
	    genemark_to_gff ($name_root, $src_prog, $gm_out, $gff_out, 
			     $config->{param_name});
	}
	
    }

} # End of for each file in the input folder

exit;

#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+


sub genemark_to_gff {

#    my ($gm_in_path, $gff_out_path, $gm_src_seq, $gm_src_prog) = @_;
    
    my ($gm_src_seq, $gm_src_prog, $gm_in_path, $gff_out_path, 
	$parameter) = @_;

    my $attribute;
    my $source = $gm_src_prog;

    #-----------------------------+
    # OPEN THE GFF OUTFILE        |
    #-----------------------------+
    # Default to STDOUT if no argument given
    if ($gff_out_path) {
	open (GFFOUT, ">$gff_out_path") ||
	    die "ERROR: Can not output gff output file:\n$gff_out_path\n"
    }
    else {
	open (GFFOUT, ">&STDOUT") ||
	    die "Can not print to STDOUT\n";
    }
    if ($gff_ver =~ "GFF3") {
	print GFFOUT "##gff-version 3\n";
    }

    #-----------------------------+
    # OPEN THE GENEMARK INFILE    |
    #-----------------------------+
    my $gm_obj;
    if ($gm_in_path) {
	$gm_obj = Bio::Tools::Genemark->new(-file => $gm_in_path);
    }
    else {
	$gm_obj = Bio::Tools::Genemark->new(-fh => \*STDIN);
    }
    
    if ($parameter) {
	$gm_src_prog = $gm_src_prog.":".$parameter;
    }

    my $rna_count = 0;
    my $gene_num = 0;

    while(my $gene = $gm_obj->next_prediction()) {
       
	$gene_num++;     
	my $gene_id = sprintf("%05d", $gene_num); # Pad the gene number
	my $gene_name = $gm_src_prog."_gene_".$gene_id."\n";
	$gene_id = "gene".$gene_id;

	# The following is a duplicate of gene_num
	$rna_count++;
	my $rna_id = sprintf("%05d", $rna_count);
	# End duplicate

	my @exon_ary = $gene->exons();

	my $num_exon = @exon_ary;

	#-----------------------------+
	# Gene span for GFF3          |
	#-----------------------------+
	if ($gff_ver =~ "GFF3") {

	    #-----------------------------+
	    # GENE                        |
            #-----------------------------+
	    $attribute = "ID=".$source."_".$gene_id;

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

	    print GFFOUT  $gm_src_seq."\t".     # Seqname
		$source."\t".                   # Source
		"gene\t".                       #feature
		"$gene_start\t".                # start
		"$gene_end\t".                  # end
		"$gene_score\t".                # score
		"$gene_strand\t".               # strand
		".\t".                          # frame
		$attribute.                     # attribute
		"\n";

	}

	#-----------------------------+
	# EXONS                       |
	#-----------------------------+
	my $exon_num = 0;
	for my $ind_exon (@exon_ary) {
	    
	    $exon_num++;
	    my $exon_id = sprintf("%05d", $exon_num); # Pad the exon number
	    $exon_id = "exon".$exon_id;

	    # Set Feature
	    my $feature = "exon";

	    #-----------------------------+
	    # GET START AND END           |
	    #-----------------------------+
	    my $start = $ind_exon->start;
	    my $end = $ind_exon->end;
	    my $strand = $ind_exon->strand;
	    if ($start > $end) {
		$end =  $ind_exon->start();
		$start = $ind_exon->end();
	    }

	    #-----------------------------+
	    # FORMAT STRAND               |
	    #-----------------------------+
	    if ($strand =~ '-1') {
		$strand = "-"; 
	    }
	    elsif ($strand =~ '1') {
		$strand = "+";
	    }
	    else {
		$strand = ".";
	    }

	    #-----------------------------+
	    # GET SCORE                   |
	    #-----------------------------+
	    # It does not look like a score is assigned
	    # but will put this here in case one shows up
	    # in the future
	    my $score = $ind_exon->score() || ".";

	    #-----------------------------+
	    # SET ATTRIBUTE               |
	    #-----------------------------+
	    if ($gff_ver =~ "GFF3") {
		$attribute = "ID=".$source."_".$gene_id."_".$exon_id.
		    "\;Parent=".$source."_".$gene_id;
	    }
	    else {
		$attribute = $gene_id;
	    }

	    #-----------------------------+
	    # GET FRAME
	    #-----------------------------+
	    #my $frame = $ind_exon->frame();
	    my $frame = ".";

	    #-----------------------------+
	    # PRINT GFF OUTPUT            |
	    #-----------------------------+
	    # The following produces output for Apollo
	    print GFFOUT $gm_src_seq."\t".   # seq name
		$gm_src_prog."\t".           # source
		$feature."\t".               # feature
		$start."\t".                 # start
		$end."\t".                   # end
		$score."\t".                 # score
		$strand."\t".                # strand
		$frame."\t".                 # frame
		$attribute.                  # attribute		
		"\n";      

	} # End of for each exon

    } # End of for each gene model

    $gm_obj->close();

    if ($gff_out_path) {
	close (GFFOUT);
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


=head1 NAME

batch_genemark.pl - Run GenMark.hmm and parse results to a gff format file. 

=head1 VERSION

This documentation refers to batch_genemark.pl version $Rev: 948 $

=head1 SYNOPSIS

=head2 Usage

    batch_genemark.pl -i DirToProcess -o OutDir

=head2 Required Arguments

    -i, --indir    # Directory of fasta files to process
    -o, --outdir   # Path to the base output directory
    -c, --config   # Path to the config file

=head1 DESCRIPTION

Run the GeneMarkHMM gene prediction program in batch mode.
Runs genmark as well as converts output to gff format.
Requires a config file to specify libraries to use.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path of the directory containing the sequences to process.

=item -o,--outdir

Path of the directory to place the program output.

=item -c, --config

Currently this program does NOT make use of a configuration file. This will
be a fairly easy thing to add, and is on my TODO List.

=back

=head1 OPTIONS

=over 2

=item --gff-ver

The GFF version for the output. This will accept either gff2 or gff3 as the
options. By default the GFF version will be GFF2 unless specified otherwise.
The default GFF version for output can also be set in the user environment
with the DP_GFF option. The command line option will always override the option
defined in the user environment. 

=item --genemark-dir

Directory that contains the GeneMark.hmm binaries. This can also be set
with the environment variable GM_BIN_DIR.

=item --lib-dir

The full path to the directory that contains the model libraries
for GeneMarkHMM. This can also be set with the environment varaible
GM_LIB_DIR.

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

=item --test

Run the program without doing the system commands.

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

The batch_genemark.pl program does not currently make use of a configuration
file. 

=head2 User Environment

This program makes use of the following variables defined in the
user's environment.

=over 2

=item GM_BIN_DIR

Directory that contains the GeneMark.hmm binaries.

=item GM_LIB_DIR

The full path to the directory that contains the model libraries
for GeneMarkHMM.

=back

The following example illustrates the ENV options set in the bash
shell.

    export GM_BIN_DIR='$HOME/apps/GenMark/genemark_hmm_euk.linux/'
    export GM_LIB_DIR='$HOME/apps/GenMark/genemark_hmm_euk.linux/'

=head1 DEPENDENCIES

=head2 Required Software

=over

=item * GeneMark.HMM

The GeneMark.HMM program is available for non-commercial Academic use
for a limited time by applying for a license at:
The http://opal.biology.gatech/edu/GeneMark

=back

=head2 Required Perl Modules

=over

=item * File::Copy

This module is required to copy the output results.

=item * Getopt::Long

This module is required to accept options at the command line.

=item * Bio::Tools::Genemark

This module is required to parse the results from the Genemark program
The module is part of the BioPerl package http://www.bioperl.org.

=back

=head1 BUGS AND LIMITATIONS

=head2 Bugs

=over 2

=item * No bugs currently known 

If you find a bug with this software, file a bug report on the DAWG-PAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=204962

=back

=head2 Limitations

=over

=item * Limited gene model supported

This program is currently limited to using the gene models that are
relevant to wheat annotation (Rice|Maize|Wheat|Barley). I will be adding
a config file option that will allow multiple gene models to be used.

=back

=head1 SEE ALSO

The batch_genemark.pl program is part of the DAWG-PAWS package of genome
annotation programs. See the DAWG-PAWS web page 
( http://dawgpaws.sourceforge.net/ )
or the Sourceforge project page 
( http://sourceforge.net/projects/dawgpaws ) 
for additional information about this package.

=head1 REFERENCE

A manuscript is being submitted describing the DAWGPAWS program. 
Until this manuscript is published, please refer to the DAWGPAWS 
SourceForge website when describing your use of this program:

JC Estill and JL Bennetzen. 2009. 
The DAWGPAWS Pipeline for the Annotation of Genes and Transposable 
Elements in Plant Genomes.
http://dawgpaws.sourceforge.net/

=head1 LICENSE

GNU GENERAL PUBLIC LICENSE, VERSION 3

http://www.gnu.org/licenses/gpl.html

THIS SOFTWARE COMES AS IS, WITHOUT ANY EXPRESS OR IMPLIED
WARRANTY. USE AT YOUR OWN RISK.

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 11/09/2007

UPDATED: 03/16/2010

VERSION: $Rev: 948 $

=cut

#-------------------------------------------------------+
# HISTORY                                               |
#-------------------------------------------------------+
#
# 11/09/2007
# - Program started to run GenMark.hmm.euk in batch mode
#   on the local machine.
# 
# 12/14/2007
# - Updated POD documentation
# - Added SVN tracking of Rev for versioning
# - Added print_help subfunction that extracts help 
#   and usage messages from the command line.
# - Added a better developed check for required 
#   arguments
#
# 03/16/2010
# - Added support for config files
# - Added support for GFF3 output
