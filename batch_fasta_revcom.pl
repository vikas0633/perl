#!/usr/bin/perl -w

use Bio::SeqIO;                # Get seq objects
use Getopt::Long;              # Get command line options

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+

# Required variables
my $indir;
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

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|indir=s"     => \$indir,
                    "o|outdir=s"    => \$outdir,
		    # ADDITIONAL OPTIONS
		    "usage"         => \$show_usage,
		    "version"       => \$show_version,
		    "man"           => \$show_man,
		    "h|help"        => \$show_help,);




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
# CREATE THE OUT DIR          |
# IF IT DOES NOT EXIST        |
#-----------------------------+
unless (-e $outdir) {
    print "Creating output dir ...\n" if $verbose;
    mkdir $outdir ||
	die "Could not create the output directory:\n$outdir";
}

#-----------------------------+
# Get the FASTA files from the|
# directory provided by the   |
# var $indir                  |
#-----------------------------+
opendir( DIR, $indir ) 
    || die "Can't open directory:\n$indir"; 
my @fasta_files = grep /\.txt$|\.fasta$|\.fa$/, readdir DIR ;
closedir( DIR );
my $num_fasta_files = @fasta_files;

if ($num_fasta_files == 0) {
    print "\a";
    print "No fasta files were found in the input direcotry:\n";
    print "$indir\n";
    print "Fasta file MUST have the fasta or fa extension to be".
	" recognized as fasta files\n";
    exit;
}


my $fasta_file_num =0;
for my $ind_fasta_file (@fasta_files) {
    
    my $ind_report_num=0;
    
    $fasta_file_num++;
    
    print STDERR "Processing $fasta_file_num of $num_fasta_files\n";

    #-----------------------------+
    # GET ROOT FILE NAME          |
    #-----------------------------+
    if ($ind_fasta_file =~ m/(.*)\.masked\.fasta$/ ) {	    
	$name_root = "$1";
    }  
    elsif ($ind_fasta_file =~ m/(.*)\.fasta$/ ) {	    
	$name_root = "$1";
    }  
    elsif ($ind_fasta_file =~ m/(.*)\.fa$/ ) {	    
	$name_root = "$1";
    } 
    elsif ($ind_fasta_file =~ m/(.*)\.txt$/ ) {	    
	$name_root = "$1";
    } 
    else {
	$name_root = "UNDEFINED";
    }


    my $in_fasta_path = $indir.$ind_fasta_file;
    my $out_fasta_path = $outdir.$ind_fasta_file;

    print STDERR "IN: $in_fasta_path\n" if $verbose;
    print STDERR "OUT: $out_fasta_path\n" if $verbose;
    
    my $seq_in = Bio::SeqIO->new('-file'   => "<$in_fasta_path");
    my $seq_out = Bio::SeqIO->new('-file' => ">$out_fasta_path");

#    open (FASTAOUT, ">$out_fasta_path") ||
#	die "Can not open output file:\n$out_fasta_path\n";

    #-----------------------------+
    # PROCESS AND FLIP            |
    #-----------------------------+
    while (my $seq = $seq_in->next_seq) {
	
	my $revcom = $seq->revcom();
	
	my $id = $seq->display_id();
	$id = $id."_rev";
	
	$revcom->display_id($id);
	$seq_out->write_seq($revcom);

    }

}

exit;
