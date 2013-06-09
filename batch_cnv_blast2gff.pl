#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_blast2gff.pl - Convert BLAST output to gff            |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_at_gmail.com                         |
# STARTED: 07/06/2006                                       |
# UPDATED: 03/02/2010                                       |
#                                                           |  
# DESCRIPTION:                                              | 
# Convert blast output to a Apollo compatible gff file.     |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;                    # Follow the rules
use Getopt::Long;              # Get options from the command line
use Bio::SearchIO;             # Parse BLAST output
use File::Copy;                # Copy the gff output to the gff dir
# The following needed for printing help
use Pod::Select;               # Print subsections of POD documentation
use Pod::Text;                 # Print POD doc as formatted text file
use IO::Scalar;                # For print_help subfunction
use IO::Pipe;                  # Pipe for STDIN, STDOUT for POD docs
use File::Spec;                # To convert a relative path to an abosolute path

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev: 948 $ =~ /(\d+)/;
# Get GFF version from environment, GFF2 is DEFAULT
my $gff_ver = uc($ENV{DP_GFF}) || "GFF2";

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+

# Set variable scope
my $indir;
my $outdir;
my $logfile;
my $name_root;
my $out_gff_path;
my $out_gff_dir;
my $out_gff_copy;
my $msg;
my $ind_blast_dir;
my $blast_file_num;
my $blast_opt = 0;

# Booleans
my $verbose = 0;
my $help = 0;
my $quiet = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_append = 0;

# Arrays
my @blast_files;


#//////////////////////
# file_num_max is the number of seqs to process in test run
my $file_num_max = 1;
my $file_num = 0;
#\\\\\\\\\\\\\\\\\\\\\\

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED
		    "i|indir=s"     => \$indir,
                    "o|outdir=s"    => \$outdir,
		    "b|blast-opt=s" => \$blast_opt,  # Has default
		    # OPTIONS
		    "logfile=s"     => \$logfile,
		    "gff-ver=s"     => \$gff_ver,
		    # BOOLEANS
		    "verbose"       => \$verbose,
		    "append"        => \$do_append,
		    "usage"         => \$show_usage,
		    "version"       => \$show_version,
		    "man"           => \$show_man,
		    "h|help"        => \$show_help,
		    "q|quiet"       => \$quiet,);

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

#-----------------------------+
# CHECK REQUIRED ARGS         |
#-----------------------------+
if ( (!$indir) || (!$outdir) ) {
    print "\a";
    print STDERR "\n";
    print STDERR "ERROR: An input directory was not specified at the".
	" command line\n" if (!$indir);
    print STDERR "ERROR: An output directory was specified at the".
	" command line\n" if (!$outdir);
    print_help("usage", $0);
}

#-----------------------------+
# OPEN THE LOG FILE           |
#-----------------------------+
if ($logfile) {
    # Open file for appending
    open ( LOG, ">>$logfile" ) ||
	die "Can not open logfile:\n$logfile\n";
    my $time_now = time;
    print LOG "==================================\n";
    print LOG "  batch_mask.pl\n";
    print LOG "  JOB: $time_now\n";
    print LOG "==================================\n";
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
    print "\a";
    print "\nERROR: No fasta files were found in the input directory\n".
	"$indir\n".
	"Fasta files must have the fasta or fa extension.\n\n";
    exit;
}

#-----------------------------+
# CREATE THE OUT DIR          |
# IF IT DOES NOT EXIST        |
#-----------------------------+
unless (-e $outdir) {
    print "Creating output dir ...\n" unless $quiet;
    mkdir $outdir
	|| die "Could not create the output directory:\n$outdir";
}


#-----------------------------+
# CONVERT BLAST TO GFF FOR    |
# EACH FILE IN THE DIR        |
#-----------------------------+
for my $ind_file (@fasta_files)
{

    $file_num++;
    $blast_file_num = 0;

    #-----------------------------+
    # GET BASE FILE NAME          |
    #-----------------------------+
    # Working on gettin this fixed
    if ($ind_file =~ m/(.*)\.masked\.fasta$/ ) {	    
	$name_root = "$1";
    }  
    elsif ($ind_file =~ m/(.*)\.fasta$/ ) {	    
	$name_root = "$1";
    }
    elsif ($ind_file =~ m/(.*)\.masked\.fa$/ ) {	    
	$name_root = "$1";
    }
    elsif ($ind_file =~ m/(.*)\.fa$/ ) {	    
	$name_root = "$1";
    } 
    else {
	$name_root = "UNDEFINED";
    }

    # PRINT PROCESS STATUS TO LOG
    print LOG "\n\nProcessing File $file_num of $num_files.\n" if $logfile;

    # PRINT PROCESS STATUS TO TERMINAL
    print STDERR 
	"\n\n+-----------------------------------------------------------+\n"
	if $verbose;
    print STDERR 
	"| Processing File $file_num of $num_files.\n" if $verbose;
    print STDERR 
	"+-----------------------------------------------------------+\n"
	if $verbose;
    
    #-----------------------------+
    # LOAD LIST OF BLAST FILES    |
    # TO PROCESS                  |
    #-----------------------------+ 
    $ind_blast_dir = $outdir.$name_root."/blast/";
    if (-e $ind_blast_dir) {
	opendir( BLASTDIR, $ind_blast_dir ) || 
	    die "Can't open directory:\n$ind_blast_dir"; 
	my @blast_files = grep /\.blo$|\.bln$|\.blx$/, readdir BLASTDIR ;
	closedir( BLASTDIR );


	#-----------------------------+
	# SET PATH TO OUTPUT GFF FILES|
	#-----------------------------+
	$out_gff_path = $ind_blast_dir.$name_root."_all_blast.gff";
	
	$out_gff_dir = $outdir.$name_root."/gff/";
	unless (-e $out_gff_dir) {
	    print "Creating output dir ...\n$out_gff_dir" if $verbose;
	    mkdir $out_gff_dir
		|| die "Could not create the output directory:\n$out_gff_dir";
	}
	
	$out_gff_copy = $out_gff_dir.$name_root."_all_blast.gff";

	
	print STDERR "IN:\t$ind_blast_dir\n";
	print STDERR "OUT:\t$out_gff_path\n";
	print STDERR "COPY:\t$out_gff_copy\n\n";

	#-----------------------------+
	# CONVERT EACH BLAST OUTFILE  |
	# TO GFF                      |
	#-----------------------------+
	$blast_file_num = 0;

	for my $ind_blast_file (@blast_files) {
	    $blast_file_num++; 

	    # For first blast output file overwrite any existing data
	    # This is done by setting the do_append boolean to 0
	    if ($blast_file_num == 1) {
		$do_append = 0;
	    }
	    else {
		$do_append = 1;
	    }
	    
	    my $blast_file_path = $ind_blast_dir.$ind_blast_file;
	    print "Converting: $ind_blast_file\n";

	    blast2gff ( $blast_file_path, $out_gff_path, 
			$do_append, $name_root, $blast_opt );
	}

	
	#-----------------------------+
	# MAKE COPY OF GFF FILE IN    |
	# GFF DIR                     |
	#-----------------------------+
	$msg = "ERROR: Could not copy file from:\n".
	    "\t$out_gff_path\n".
	    "\t$out_gff_copy\n";
	copy ($out_gff_path, $out_gff_copy)
	    || print LOG "$msg";

    }
    else {
	print LOG "ERROR: Could not find BLAST dir:\n$ind_blast_dir\n";
    }

    # If max file number then exit
    # This is for debug tests
    #if ($file_num = $file_num_max) {exit;}

}

close LOG if $logfile;

exit;

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



# REMOVING THE FOLLOWING FOR THE bast2gff function 
# that is GFF3 capable

sub blast2gff {
# CONVERT BLAST TO GFF 
    
    # seqname - ID to use in the seqname field
    # blastin - path to the blast input file
    # gffout  - path to the gff output file
    # append  - boolean append data to existing file at gff out
    # bopt    - The type of blast report to parse (0,8,9)
    my ($blastin, $gffout, $append, $seqname, $bopt ) = @_;
    my $blastprog;        # Name of the blast program used (blastn, blastx)
    my $dbname;           # Name of the database blasted
    my $hitname;          # Name of the hit
    my $start;            # Start of the feature
    my $end;              # End of the feature
    my $strand;           # Strand of the hit
    my $blast_report;     # The handle for the blast report
    my $blast_score;      # The score for the blast hit

    # Temp change feature from exon to the more approprate 'match'
    my $feature = "match";

    # Remove prohibited characters
    if ($gff_ver =~ "GFF3") {
	$seqname = seqid_encode($seqname);
    }

    my $seq_name_len = length($seqname);

    #-----------------------------+
    # GET BLAST REPORT            |
    #-----------------------------+
    if ($bopt == 8 || $bopt == 9) {

	# PARSE m8 or m9 ALIGNMENTS
	$blast_report = new Bio::SearchIO ( '-format' => 'blasttable',
					    '-file'   => $blastin)
	    || die "Could not open BLAST input file:\n$blastin.\n";
	
    }
    else {

	# PARSE DEFAULT ALIGNMNETS (m0)
	$blast_report = new Bio::SearchIO ( '-format' => 'blast',
					    '-file'   => $blastin) 
	    || die "Could not open BLAST input file:\n$blastin.\n";
	
    }

    #-----------------------------+
    # OPEN GFF OUTFILE            |
    #-----------------------------+
    # Default to STDOUT if no argument given
    if ($gffout) {
	if ($append) {
	    open (GFFOUT, ">>$gffout") 
		|| die "Can not open file:\n $gffout\n";
	}
	else {
	    open (GFFOUT, ">$gffout") 
		|| die "Can not open file:\n $gffout\n";
	}
    }
    else {
	open (GFFOUT, ">&STDOUT") ||
	    die "Can not print to STDOUT\n";
    }
    
    while (my $blast_result = $blast_report->next_result())
    {
	$blastprog = $blast_result->algorithm;

	# FOR m8 output the database name does not work
	# so I will use the file name as given by blastin

	if ($bopt == 8 || $bopt == 9) { 
	    $dbname = $blastin;

	    # Since DAWG-PAWS uses a specific naming stragety, I will
	    # try to extract the database name from the file name
	    # Basic pattern is 'path_root/file name'_'database name'
	    # can get the 
	    if ($dbname =~ m/(.*)\/(.*)\.bln$/ ) {
		$dbname = $2;
		$dbname = substr($dbname, $seq_name_len+1 );

	    }
	    elsif ($dbname =~ m/(.*)\/(.*)\.blx$/ ) {
		$dbname = $2;
		$dbname = substr($dbname, $seq_name_len+1);
	    }

	}
	else {
	    $dbname = $blast_result->database_name();
	}

	# remove trailing white space
	$dbname =~ s/\s+$//;
	if ($gff_ver =~ "GFF3") {
	    $dbname = gff3_encode($dbname);
	}

    	while (my $blast_hit = $blast_result->next_hit())
	{

	    while (my $blast_hsp = $blast_hit->next_hsp())
	    {

		my $hitname = $blast_hit->name();
		if ($gff_ver =~ "GFF3") {
		    $hitname = gff3_encode($hitname);
		}
		
		$strand = $blast_hsp->strand('query');
		
		if ($strand =~ "-1") {
		    $strand = "-";
		}
		elsif ($strand =~ "1") {
		    $strand = "+";
		}
		else {
		    die "Error parsing strand\n";
		}

		#-----------------------------+
		# GET QRY START AND END       |
		#-----------------------------+
		# Make certain that start coordinate is
		# less then the end coordinate
		if ( $blast_hsp->start() < $blast_hsp->end() ) {
		    $start = $blast_hsp->start();
		    $end = $blast_hsp->end();
		}
		else {
		    $start = $blast_hsp->end();
		    $end = $blast_hsp->start();
		}
		
		#-----------------------------+
		# GET HIT START AND END       |
		#-----------------------------+
		my $hit_start;
		my $hit_end;
		if ( $blast_hsp->start('hit') < $blast_hsp->end('hit') ) {
		    $hit_start = $blast_hsp->start('hit');
		    $hit_end = $blast_hsp->end('hit');
		}
		else {
		    $hit_start = $blast_hsp->end('hit');
		    $hit_end = $blast_hsp->start('hit');
		}

		#-----------------------------+
		# GET BLAST SCORE             |
		#-----------------------------+
		if ($bopt == 8 || $bopt == 9) {
		    #$blast_score = ".";
		    # trying bits
		    $blast_score = $blast_hsp->bits();
		}
		else {
		    $blast_score = $blast_hsp->score();
		}

		#-----------------------------+
		# PRINT OUTPUT TO GFF FILE    |
		# WITH BAC DATA               |
		#-----------------------------+
		# Changing BLASTN to the Bac Name appears to allow 
		# these to be drawn on different levels.
		#-----------------------------+
		# Set attribute               |
		#-----------------------------+
		my $attribute;
		#my $feature;
		my $source = $blastprog.":".$dbname;   
		if ($gff_ver =~ "GFF3") {
		    $feature = "match_part"; 
		    $source =~ s/\s+$source//;
		    $attribute = "ID=".$source."_".$hitname.
			";".
			"Name=".$hitname.";".
			"Target=".$hitname." ".
			$hit_start." ".$hit_end;
		}
		else {
		    $attribute = $hitname;
		}

		print GFFOUT 
		    "$seqname\t".                            # Seqname
		    $blastprog.":".$dbname."\t".                  # Source
		    "$feature\t".                            # Feature type name
		    "$start\t".                              # Start
		    "$end\t".                                # End
		    $blast_score."\t".                       # Score
		    "$strand\t".                             # Strand
		    ".\t".                                   # Frame
		    "$attribute\n";                          # Feature name

		# PRINT GFF TO STDERR IF IN VERBOSE MODE
		if ($verbose) {
		    print STDERR "\t   SEQ:\t$seqname\n";
		    print STDERR "\t SOURC:\t$blastprog:$dbname\n";
		    print STDERR "\t START:\t$start\n";
		    print STDERR "\t   END:\t$end\n";
		    print STDERR "\t SCORE:\t".$blast_score."\n";
		    print STDERR "\tSTRAND:\t$strand\n";
		    print STDERR "\t   HIT:\t$hitname\n";
		}

	    } # End of while next hsp
	} # End of while next hit
    } # End of while next result
    
    close GFFOUT;
    
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

Batch_cnv_blast2gff.pl - Convert blast output to GFF format

=head1 VERSION

This documentation refers to program version $Rev: 948 $

=head1 SYNOPSIS

=head2 Usage

    batch_cnv_blast2gff.pl -i InDir -o OutDir

=head2 Required Arguments

    -i, --indir    # Directory of fasta files to process
    -o, --outdir   # Path to the base output directory

=head1 DESCRIPTION

Convert NCBI BLAST output to an Apollo compatible GFF file. This can produce
a separate gff file for each BLAST report, or can merge all BLAST results
into a single GFF file.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path of the directory containing the sequences to process.

=item -o,--outdir

Path of the directory to place the program output.

=back

=head1 OPTIONS

=over 2

=item --logfile

Path to a file that will be used to log program status.
If the file already exists, additional information will be concatenated
to the existing file

=item --append

Append results to existing gff file.

=item --usage

Print a short overview of how to use program from the command line.

=item --help

Print a short program usage state with a summary of options.

=item --version

Show program version. This will print the SVN version of the script.

=item --man

Show the full program manual. This uses the perldoc command to print the 
POD documentation for the program.

=item -q,--quiet

Run the program with minimal output.

=item --verbose

Run the program with maximum output to the screen. This will provided a 
detailed status of the progress of the program.

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

This program does not depend on any configuration files or environmental
settings.

=head1 DEPENDENCIES

=head2 Required Perl Modules

=over

=item * Bio::SearchIO;

This module is part of the BioPerl package of programs. It is used to
parse the BLAST output. For information on downloading and installing
BioPerl see http://www.bioperl.org/wiki/Main_Page

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

=head2 Limitations

=over

=item * Limited to NCBI BLAST

The current version is limited to using the NCBI version of BLAST.

=item * Limited file extensions are supported

BLAST output file must currently end with blo, bln, or blx. For example
a BLASTx output may be named BlastOut.blx while a BLASTN output
may be names BlastOut.bln. FASTA files must end with a fasta or fa extension.
For examples must have names like my_seq.fasta or my_seq.fa.

=back

=head1 SEE ALSO

The batch_blast.pl program is part of the DAWG-PAWS package of genome
annotation programs. See the DAWG-PAWS web page 
( http://dawgpaws.sourceforge.net/ )
or the Sourceforge project page 
( http://sourceforge.net/projects/dawgpaws ) 
for additional information about this package.

=head1 LICENSE

GNU GENERAL PUBLIC LICENSE, VERSION 3

http://www.gnu.org/licenses/gpl.html   

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 08/06/2007

UPDATED: 03/02/2010

VERSION: $Rev: 948 $

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 08/06/2007
# -Program started, all base subfunctions added
#
# 12/04/2007
# - POD Documentation updated
# - Version number changed to SVN Revision Number
# - Changed terminal output to STDERR
# - Moved POD documentation to end of the code
# - Trying to add POD select to print help and usage
#   message from the POD documentation
# - Addd an end statement
#
# 12/05/2007
# - Created the print_help as a subfunction that
#   takes type of help to print and POD source as 
#   variables.
# - Added 
#
# 01/04/2008
# - Fixed stupid typo where I had if instead of elseif
#
# 08/05/2008
# - Adding option to parse -m8 or -m9 BLAST output
#   The will be implemented as the blast-option
#   variable. The BLAST option will use the numbers 
#   as indicated in BLAST
#     0 = pairwise, 
#     1 = query-anchored showing identities,
#     2 = query-anchored no identities,
#     3 = flat query-anchored, show identities,
#     4 = flat query-anchored, no identities,
#     5 = query-anchored no identities and blunt ends,
#     6 = flat query-anchored, no identities and blunt ends,
#     7 = XML Blast output,
#     8 = tabular,
#     9 tabular with comment lines
#     10 ASN, text
#     11 ASN, binary [Integer]
#  Currently only 8,9 and 0 are implemented. 
#
# 03/02/2010
# - Adding support for GFF3 format output
# 03/03/2010
# - Finalized support for GFF3 format output
