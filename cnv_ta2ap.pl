#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_ta2ap.pl - Convert TriAnnot GFF3 to GFF output        |
#                                                           |
#-----------------------------------------------------------+ 
#                                                           |
#  AUTHOR: James C. Estill                                  |
# STARTED: 02/09/2007                                       |
# UPDATED: 12/11/2007                                       |
# DESCRIPTION:                                              |
#  Converts gff data tracks from gff format to the game     |
#  xml format for use in the Apollo Genome Annotation       |
#  Curation program. This is for use with output from the   |
#  TriAnnot pipeline.                                       |
#                                                           |
# DEPENDENCIES:                                             |
# *Apollo                                                   |
#    Requires that apollo be installed on the local machine |
#    since apollo is being used as the engine to do the     | 
#    conversion between formats.                            |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+ 

# TO DO:
# [ ] Grep the GFF directory to load the array of GFF files
#     into the Files2Convert array
# [ ] DON'T HARD CODE APOLLO PATH USE ENV VAR

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

#-----------------------------+
# SET VARIABLE SCOPE          |
#-----------------------------+
my $InFastaFile;               # FASTA fiel that GFF files will be mapped onto
my $InputRoot;                 # Base directory that contains all the GFF files
my $OutputRoot;                # Output directory
my $CatOut;                    # Concatenated output of all gff files,
                               #  not required
#my %Options;                   # Options hash to hold command line options
my $Valid = '0';               # Boolean [0,1] : A valid gff file was 
                               # found for conversion


# BOOLEANS
my $show_help = 0;             # Show program help
my $show_version = 0;          # Show program version
my $show_man = 0;              # Show program manual page using peldoc
my $show_usage = 0;            # Show program usage command             
my $quiet = 0;                 # Boolean for reduced output to STOUT
my $test = 0;                  # Run the program in test mode
my $verbose = 0;               # Run the program in verbose mode
my $create_game = 0;           # Create game.xml output file

# Command vars with default values
my $ap_path = "apollo";        # Apollo binary, assumes this is in the
                               # user's Path

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(
		    # Required Arguments
                    "g|gffdir=s"   => \$InputRoot,
		    "o|outdir=s"   => \$OutputRoot,
		    "catout=s"     => \$CatOut,
		    # Optional strings
		    "ap-path=s"    => \$ap_path,
		    "i|infile=s"   => \$InFastaFile,
		    # Booleans
		    "game"         => \$create_game,
		    "verbose"      => \$verbose,
		    "test"         => \$test,
		    "usage"        => \$show_usage,
		    "version"      => \$show_version,
		    "man"          => \$show_man,
		    "h|help"       => \$show_help,
		    "q|quiet"      => \$quiet,);

print STDERR "The ApConvert pogram has started.\n" if $verbose;

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
# CHECK FOR EXISTENCE OF      |
# INPUT FILES AND OUTPUT DIR  |
#-----------------------------+
# Check for  input FASTA file if it is passed at the 
# command line
if ($InFastaFile) { 
    unless (-e $InFastaFile) {
	print "The fasta file does not exist\n:".$InFastaFile."\n";
	exit;
    }
}

# Input GFF directory
unless (-e $InputRoot) {
    print "The input directory can not be found at:\n";
    print "$InputRoot\n";
    exit;
}

# MAKE OUTPUT DIR IF IT DOES NOT EXIST
unless (-e $OutputRoot) {
    mkdir $OutputRoot, 0777;
}

#-----------------------------+
# LOAD GFF FILE LIST TO ARRAY |
#-----------------------------+
opendir(INDIR, $InputRoot) ||
    die "Can not open input directory:\n$InputRoot: $!";
my @Files2Convert = grep /gff$/, readdir INDIR;
closedir INDIR;

#-----------------------------+
# OPEN CATOUTPUT              |
#-----------------------------+
if ($CatOut)
{
    open (CATOUT, ">".$CatOut) ||
	die "Can not open concatenated gff output file:\n$CatOut\n";
}

my $NumFiles = $#Files2Convert;

# Remove gff from Files2Convert
for (my $i=0; $i<=$NumFiles; $i++)
{
    if ($Files2Convert[$i] =~ /(.*)\.gff/)
    {
	print $1."\n";
	$Files2Convert[$i] = $1;
    }
} # End of for i loop

foreach my $FileRoot (@Files2Convert)
{   

    my $TriGffPath = $InputRoot.$FileRoot.".gff";      # TriAnnot Gff
    my $ApGffPath = $InputRoot.$FileRoot.".ap.gff";    # Apollo Gff
    my $OutPath = $OutputRoot.$FileRoot.".game.xml";   # Game XML output 

    # SHOW PROGRAM STATUS
    print "Processing: $TriGffPath\n";
    
    #-----------------------------+
    # CHECK FOR FILE EXISTENCE    |
    # BEFORE WE TRY TO USE IT     |
    #-----------------------------+
    if (-e $TriGffPath)
    {
	
	#-----------------------------+
	# CONVERT TRIANNOT GFF TO     |
	# APOLLO USABLE GFF           |
	#-----------------------------+
	# This will need to be a different convert for each
	# of the input files
	
	#---------------+
	# FGENESH       |
	#---------------+
	if ($FileRoot =~ "2fGh") {
	    &ModFGENESH ( $TriGffPath, $ApGffPath );
	    $Valid = '1';
	}
	
	#---------------+
	# GeneMark.hmm  |
	# Ta Matrix     |
	#---------------+
	elsif ($FileRoot =~ "2gmTa") {
	    &ModGenMark ( $TriGffPath, $ApGffPath );
	    $Valid = '1';
	}
	
	#---------------+
	# GeneMark.hmm  |
	# Os Matrix     |
	#---------------+
	elsif ($FileRoot =~ "2gmOs") {
	    &ModGenMark ( $TriGffPath, $ApGffPath );
	    $Valid = '1';
	}
	
	#---------------+
	# GeneMark.hmm  |
	# Hv Matrix     |
	#---------------+
	elsif ($FileRoot =~ "2gmHv") {
	    &ModGenMark ( $TriGffPath, $ApGffPath );
	    $Valid = '1';
	}
	
	#---------------+
	# GeneMark.hmm  |
	# Zm Matrix     |
	#---------------+
	elsif ($FileRoot =~ "2gmZm") {
	    &ModGenMark ( $TriGffPath, $ApGffPath );
	    $Valid = '1';
	}
	
	#---------------+
	# GeneId        |
	#---------------+
	elsif ($FileRoot =~ "2gID") {
	    &ModGeneId ( $TriGffPath, $ApGffPath );
	    $Valid = '1';
	}
	
	#---------------+
	# Eugene Os     |
	#---------------+
	elsif ($FileRoot =~ "2eugOs") {
	    &ModEugene ( $TriGffPath, $ApGffPath );
	    $Valid = '1';
	}

	#---------------+
	# TRF           |
	#---------------+
	# Tandem Repeat Finder
	elsif ($FileRoot =~ "1trf") {
	    &ModTRF ( $TriGffPath, $ApGffPath );
	    $Valid = '1';
	}

	#---------------+
	# UNSUPPORTED   |
	#---------------+
	else {
	    $Valid = '0';
	    print "\a\nERROR\n";
	    print "The file $FileRoot is not currently supported.\n";
	} # End of FileRoot selection of parser
	
	#-----------------------------+
	# Convert each gff file game  |
	# xml using the ApolloConvert |
	# subfunction                 |
	#-----------------------------+
	# Will comment this out while I am working
	# on the FGENESH and GenMark conversion subfunctions
	if ($Valid == '1') {

	    # Convert GFF file to game format
	    if ($create_game) {
		&ApolloConvert ($ApGffPath, "gff", $OutPath, 
				"game", $InFastaFile, "NULL", $ap_path);
	    }
	    
	    # Add apollo formatted GFF output to concatenated
            # GFF output file
	    if ($CatOut)
	    {
		open (APIN,$ApGffPath);
		while (<APIN>)
		{
		    print CATOUT $_;
		}
		close APIN;
	    } # End of if CatOut
	} # End of if Valid is true

    }
    else {
	print "\a\aERROR: The input Gff file does not exist at:".
	    "\n$TriGffPath\n";
    } # End of if input file exists

} # End of for each file in Files2Convert

#-----------------------------+
# CONVERT CONCATENATED GFF    |
# FILE TO THE APOLLO FORMAT   |
#-----------------------------+
if ($CatOut) {    
    if ($create_game) {
	print "Converting the concatenated GFF file\n";
	close CATOUT;
	
	my $CatGameOut = $CatOut.".game.xml";
	
	&ApolloConvert ($CatOut, "gff", $CatGameOut, 
			"game", $InFastaFile, "NULL", $ap_path);
    }
}

exit;

#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub ApolloConvert {
#-----------------------------+
# CONVERT AMONG FILE FORMATS  |
# USING THE APOLLO PROGRAM    |
#-----------------------------+
# Converts among the various data formats that can be used 
# from the command line in tbe Apollo program. For example
# can convert GFF format files into the game XML format.
# NOTES:
#  - Currently assumes that the input file is in the correct
#    coordinate system.
#  - GFF files will require a sequence file
#  - ChadoDB format will require a db password


    # ApPath - the path of dir with the Apollo binary
    #          Specifying the path will allow for cases
    #          where the program is not in the PATHS
    # ApCmd  - the apollo commands to run

    my $InFile = $_[0];        # Input file path
    my $InForm = $_[1];        # Output file format:
                               # game|gff|gb|chadoxml|backup
    my $OutFile = $_[2];       # Output file path
    my $OutForm = $_[3];       # Ouput file foramt
                               # chadoDB|game|chadoxml|genbank|gff|backup
    my $SeqFile = $_[4];       # The path of the sequence file
                               # This is only required for GFF foramt files
                               # When not required this can be passed as na
    my $DbPass = $_[5];        # Database password for logging on to the 
                               # chado database for reading or writing.
    my $ApPath = $_[6];        # Path to the Apollo binary
    my ( $ApCmd );


# The following path is for the jlb10 machine
#    $ApPath = "/home/jestill/Apps/Apollo/Apollo";
#    $ApPath = "/home/jestill/Apps/Apollo_1.6.5/apollo/bin/apollo";
#    $ApPath = "/Applications/Apollo/bin/apollo";

    # Set the base command line. More may need to be added for different
    # formats. For example, GFF in requires a sequence file and the CHADO
    # format will require a database password.
    $ApCmd = $ApPath." -i ".$InForm." -f ".$InFile.
	" -o ".$OutForm." -w ".$OutFile;

    # Make sure that that input output formats are in lowercase
    # may need to add something here to avoid converting chadoDB
    $InForm = lc($InForm);
    $OutForm = lc($OutForm);
    
    # Determine the proper command to use based on the input format
    # since GFF file also require a sequence file
    if ($InForm =~ "gff" )
    {
	$ApCmd = $ApCmd." -s ".$SeqFile;
    }
    
    if ($InForm =~ "chadodb")
    {
	$ApCmd = $ApCmd." -D ".$DbPass;
    }

    # Do the apollo command
    system ( $ApCmd );

}


sub ModFGENESH {
#-----------------------------+
# MODIFY THE FGENESH GFF FILE |
# TO WORK WITH APOLLO         |
#-----------------------------+
# InFile is the full path to a valid gff infile
# OutFile is the full path for the modified file

    my $InFile = $_[0];
    my $OutFile = $_[1];
    my $LineCount = 0;  # LineCount used for debu
    
    open (IN, "<$InFile") ||
	die "Can not open infile:\n$InFile\n";
    open (OUT, ">$OutFile") ||
	die "Can not open outfile:\n$OutFile\n";
    
    while (<IN>)
    {
	$LineCount++;
	chomp;                  # Remove the newline char
	my @Sp = split(/\t/);   # Split by the tab character

	if ($Sp[2] =~ "exon")
	{

	    # Split out the name to use for the
	    # exon
	    my @Id = split(/_/, $Sp[8]);

	    # Print to screen for debug
	    print $LineCount."\t".$Sp[2]."\n";
	    print "\t".$Id[1]."\n";

	    # PRINT OUTPUT TO GFF FILE
	    print OUT $Sp[0]."\t";
	    print OUT $Sp[1]."\t";
	    print OUT $Sp[2]."\t";
	    print OUT $Sp[3]."\t";
	    print OUT $Sp[4]."\t";
	    print OUT $Sp[5]."\t";
	    print OUT $Sp[6]."\t";
	    print OUT $Sp[7]."\t";
	    print OUT $Id[1]."\n";

	}
	

    }

    close IN;
    close OUT;

}



sub ModGenMark {
#-----------------------------+
# MODIFY THE FGENESH GFF FILE |
# TO WORK WITH APOLLO         |
#-----------------------------+
# InFile is the full path to a valid gff infile
# OutFile is the full path for the modified file

    my $InFile = $_[0];
    my $OutFile = $_[1];
    my $LineCount = 0;  # LineCount used for debu
    
    open (IN, "<$InFile") ||
	die "Can not open infile:\n$InFile\n";
    open (OUT, ">$OutFile") ||
	die "Can not open outfile:\n$OutFile\n";
    
    while (<IN>)
    {
	$LineCount++;
	chomp;                  # Remove the newline char
	my @Sp = split(/\t/);   # Split by the tab character

	if ($Sp[2] =~ "exon")
	{

	    # Split out the name to use for the
	    # exon
	    my @Id = split(/;/, $Sp[8]);
	    my @ModId = split(/_/, $Id[0]);
	    # Print to screen for debug
	    #print "\t".$Id[0]."\n";
	    print $LineCount."\t".$ModId[2]."\n";

	    # PRINT OUTPUT TO GFF FILE
	    print OUT $Sp[0]."\t";
	    print OUT $Sp[1]."\t";
	    print OUT $Sp[2]."\t";
	    print OUT $Sp[3]."\t";
	    print OUT $Sp[4]."\t";
	    print OUT $Sp[5]."\t";
	    print OUT $Sp[6]."\t";
	    print OUT $Sp[7]."\t";
	    print OUT $ModId[2]."\n";

	}
	

    }

    close IN;
    close OUT;

}

sub ModEugene {
#-----------------------------+
# MODIFY THE EUGENE GFF FILE  |
# TO WORK WITH APOLLO         |
#-----------------------------+
# InFile is the full path to a valid gff infile
# OutFile is the full path for the modified file

    my $InFile = $_[0];
    my $OutFile = $_[1];
    my $LineCount = 0;  # LineCount used for debu
    
    open (IN, "<$InFile") ||
	die "Can not open infile:\n$InFile\n";
    open (OUT, ">$OutFile") ||
	die "Can not open outfile:\n$OutFile\n";
    
    while (<IN>)
    {
	$LineCount++;
	chomp;                  # Remove the newline char
	my @Sp = split(/\t/);   # Split by the tab character

	if ($Sp[2] =~ "exon")
	{

	    # Split out the name to use for the
	    # exon
	    my @Id = split(/;/, $Sp[8]);
	    my @ModId = split(/_/, $Id[0]);
	    # Print to screen for debug
	    #print "\t".$Id[0]."\n";
	    print $LineCount."\t".$ModId[4]."\n";

	    # PRINT OUTPUT TO GFF FILE
	    print OUT $Sp[0]."\t";
	    print OUT $Sp[1]."\t";
	    print OUT $Sp[2]."\t";
	    print OUT $Sp[3]."\t";
	    print OUT $Sp[4]."\t";
	    print OUT $Sp[5]."\t";
	    print OUT $Sp[6]."\t";
	    print OUT $Sp[7]."\t";
	    print OUT "transcript_".$ModId[4]."\n";

	}
	

    }

    close IN;
    close OUT;

}

sub ModGeneId {
#-----------------------------+
# MODIFY THE FGENESH GFF FILE |
# TO WORK WITH APOLLO         |
#-----------------------------+
# InFile is the full path to a valid gff infile
# OutFile is the full path for the modified file

    my $InFile = $_[0];
    my $OutFile = $_[1];
    my $LineCount = 0;  # LineCount used for debu
    
    open (IN, "<$InFile") ||
	die "Can not open infile:\n$InFile\n";
    open (OUT, ">$OutFile") ||
	die "Can not open outfile:\n$OutFile\n";
    
    while (<IN>)
    {
	$LineCount++;
	chomp;                  # Remove the newline char
	my @Sp = split(/\t/);   # Split by the tab character

	if ($Sp[2] =~ "exon")
	{

	    # Split out the name to use for the
	    # exon
	    my @Id = split(/_/, $Sp[8]);

	    # Print to screen for debug
	    print $LineCount."\t".$Sp[2]."\n";
	    print "\t".$Id[1]."\n";

	    # PRINT OUTPUT TO GFF FILE
	    print OUT $Sp[0]."\t";
	    print OUT $Sp[1]."\t";
	    print OUT $Sp[2]."\t";
	    print OUT $Sp[3]."\t";
	    print OUT $Sp[4]."\t";
	    print OUT $Sp[5]."\t";
	    print OUT $Sp[6]."\t";
	    print OUT $Sp[7]."\t";
	    print OUT $Id[1]."\n";

	}
	

    }

    close IN;
    close OUT;

}



sub ModTRF {
#-----------------------------+
# MODIFY THE TRF GFF FILE     |
# TO WORK WITH APOLLO         |
#-----------------------------+
# TRF: Tandem Repeat Finder
# InFile is the full path to a valid gff infile
# OutFile is the full path for the modified file

    my $InFile = $_[0];
    my $OutFile = $_[1];
    my $LineCount = 0;  # LineCount used for debug
    
    open (IN, "<$InFile") ||
	die "Can not open infile:\n$InFile\n";
    open (OUT, ">$OutFile") ||
	die "Can not open outfile:\n$OutFile\n";
    
    while (<IN>)
    {
	$LineCount++;
	chomp;                  # Remove the newline char
	my @Sp = split(/\t/);   # Split by the tab character

	if ($Sp[2] =~ "tandem_repeat")
	{

	    # Split out the name to use for the
	    # exon
	    my @Id = split(/;/, $Sp[8]);
	    # SmID for Small ID, Additional information
	    # can be parsed from the TriAnnot pipeline GFF file
	    # including the length and repeat motif
	    my @SmId = split(/_/, $Id[0]);

	    # Print to screen for debug
	    print $LineCount."\t".$Sp[2]."\n";
	    print "\tTRF".$SmId[2]."\n";

	    # PRINT OUTPUT TO GFF FILE
	    print OUT $Sp[0]."\t";
	    print OUT $Sp[1]."\t";
	    print OUT $Sp[2]."\t";
	    print OUT $Sp[3]."\t";
	    print OUT $Sp[4]."\t";
	    print OUT $Sp[5]."\t";
	    print OUT $Sp[6]."\t";
	    print OUT $Sp[7]."\t";
	    # Small Id is just TRF####
	    print OUT "TRF".$SmId[2]."\n";
	}
	

    }

    close IN;
    close OUT;

} # End of ModTRF



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


# THE FOLLOWING IS DEPRECATED

sub PrintHelp  {
#-----------------------------+
# PRINT HELP STATEMENT        |
#-----------------------------+
    my $FullUsage = "\n".
	"+-----------------------------+\n".
	"| ApConvert                   |\n". 
	"| v. 04/12/2007               |\n".
	"+-----------------------------+\n".
	"Usage:\n".
	"ApConvert.pl -i InFastaFile -g GffDir -o OutputDir\n".
	"\n".
	"+-----------------------------+\n".
	"| ARGUMENTS                   |\n".
	"+-----------------------------+\n".
	" -i Fasta File [String]\n".
	"    Full path to FASTA file that GFF results refererence\n".
	" -g Gff Directory [String]\n".
	"    Full path to directory containing the multiple GFF files.\n".
	" -o Output Directory [String]\n".
	"    Full path to place the GAME xml files\n".
	"    The directory will be created if it does not exist\n".
	" -O Concatenated output file in gff format [String]\n".
	" -h Print help statement [Boolean]\n".
	" \n".
	"+-----------------------------+\n".
	"| SUPPORTED TriAnnot FILES    |\n".
	"+-----------------------------+\n".
	" 1trf ---> Tandem Repeat Finder\n".
	" 2eugOs -> Eugene - Oryza sativa\n".
	" 2fGh ---> Fgenesh - Gene Model Prediction\n".
	" 2gID ---> GeneID - Gene Model Prediction\n".
	" 2gmHv --> GeneMark Hordeum vulgare\n".
	" 2gmOs --> GeneMark Oryza sativa\n".
	" 2gmTa --> GeneMark Triticum aestevum\n".
	" 2gmZm --> GeneMakr Zea mays\n".  
	"\n";
    print $FullUsage;
}

=head1 NAME

cnv_ta2ap.pl - Convert TriAnnot GFF3 files to Apollo format

=head1 SYNOPSIS

=head2 Usage

    cnv_ta2ap.pl -i InFastaFile -g GffDir -o OutputDir

=head2 Required Arguments

    -g, --gffdir   # Directory containing TriAnnot gff files
    -o, --outdir   # Path to the base output directory

=head1 DESCRIPTION

Converts gff data tracks from gff format to the game
xml format for use in the Apollo Genome Annotation
Curation program. This is for use with output from the
TriAnnot pipeline.

=head1 REQUIRED ARGUMENTS

=over

=item -g,--gffdir

This is the directory that contains the gff output from the TriAnnot
program that you want to convert.

=item -o,--outdir

Path of the directory to place the program output.

=back

=head1 OPTIONS

=over 2

=item --catout

Optional path to a single gff output file that will contain all
of the converted files concatenated into a single gff file.

=item -i,--infile

Input FASTA file. The input fasta file will be used as the base file
to map gff sequence features onto when converting from gff to
game.xml output.

=item --ap-path

Full path to the binary file for the Apollo Genome Annotation
Curation program.

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

=back

=head1 DIAGNOSTICS

Error messages generated by this program and possible solutions are listed
below.

=over 2

=item ERROR: Could not create the output directory

The output directory could not be created at the path you specified. 
This could be do to the fact that the directory that you are trying
to place your base directory in does not exist, or because you do not
have write permission to the directory you want to place your file in.

=back

=head1 CONFIGURATION AND ENVIRONMENT

The program cnv_ta2ap.pl does not currently use an external configuration
file or make use of variables set the user's environment.

=head1 DEPENDENCIES

=head2 Required Software

=over

=item * Apollo Genome Annotation Tool

To convert output from gff output to the game.xml file format you
will need to have a version of Apollo Genome Annotation program
on your local machine AND have access to the Apollo GUI. Apollo 
can be downloaded at:
http://www.fruitfly.org/annot/apollo/

=back

=head2 Required Perl Modules

=over

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

=item * Limited support for TriAnnot output

Conversion subfunctions are specific to the data type
being converted. Since the DAWG-PAWS process makes limited use of
the TriAnnot output, support is limited to the following
output files.

=over

=item 1trf 

Tandem Repeat Finder

=item 2eugOs 

Eugene - Oryza sativa.

=item 2fGh 

Fgenesh - Gene Model Prediction

=item 2gID 

GeneID - Gene Model Prediction.

=item 2gmHv 

GeneMark Hordeum vulgare

=item 2gmOs 

GeneMark Oryza sativa.

=item 2gmTa 

GeneMark Triticum aestevum.

=item 2gmZm 

GeneMark Zea mays

=back

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

STARTED: 02/09/2007

UPDATED: 12/11/2007

VERSION: $Rev: 948 $

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
# 
# 02/09/2007
#  - Program started basic funnctionality of main body
#  - Added the ApolloConvert subfunction from pre-existing
#    source
#
# 02/12/2007
#  - Added ModFGENESH subfunction to convert the TriAnnot
#    GFF file to a format that Apollo can use.
#  - Added ModGenmark subfunction to convert the GenMark
#    files to the format Apollo can ues.
#  - Added code to the body of the script to use the
#    ModFGENESH and ModGenmark subfunctions
#
# 02/13/2007
#  - Added ModTRF subfunction. The names are being
#    lost when this GFF file is parsed, but leaving
#    this alone for now.
#
# 02/19/2007
# - Working on a modification to read fasta files
#   and load to a Files to convert. This will need
#   to also be aware of where the gff files are
#   for the fasta file of interest.
#
# 04/12/2007
# - Changing variables from hard coded to command line
# - Added the PrintHelp subfunction
#
# 04/16/2007
# - Get list of gff files from the input directory
# - Added Eugene format conversion
# - Added option for concatenated GFF output as $CatOut
# - Added option to convert concate GFF to GAME xml
#
# 06/22/2007
# - Adding POD documentation and cleaning up code
#
# 12/11/2007
# - Moved POD documentation to the end of the script
# - Updated POD documentatin
# - Added print_help subfunction
# - Deprecated PrintHelp subfunction
# - Added SVN Rev tracking
# - Converted options vars to long format
# - Added use strict and cleaned up problems
# - Commented out all refernces to the ApolloConvert
#   subfunction.
# - The fasta file is only used for conversion to game.xml
# - Added the create_game boolean
# - Added option to provide the path to apollo at the command line
# - Wow I hope nothing broke
