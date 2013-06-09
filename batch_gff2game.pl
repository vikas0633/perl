#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_convert.pl - Convert gff to game xml format         |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_at_gmail.com                         |
# STARTED: 08/07/2007                                       |
# UPDATED: 08/09/2007                                       |
#                                                           |
# DESCRIPTION:                                              | 
# Short script to convert and copy the wheat BACs           |
# Run this in the parent dir that the HEX* dirs exist       |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+
# This program needs a lot of cleaning up.
# Switch to a config file or fasta input list

use File::Copy;

print "Begin the batch conversion has\n";

# ap_path is the path to apollo on your local machine

my $ap_path = "/home/jestill/Apps/Apollo_1.6.5/apollo/bin/apollo";
my @root_names = (
		  #-----------------------------+
		  # CONVERTED
		  # 08/28/2007
		  #-----------------------------+
		  "HEX0014K09",
		  "HEX0238L06",
		  "HEX0298C05",
		  "HEX0325M12",
		  "HEX0350E24",
		  "HEX0358K19",
		  "HEX0411K13",
		  "HEX0412L20",
		  "HEX0454N17",
		  "HEX0593O07",
		  "HEX0651O10",
		  "HEX0693B09",
		  "HEX0720B04",
		  "HEX0769A02",
		  "HEX0769A15",
		  "HEX0769B12",
		  "HEX0769B23",
		  "HEX0769C08",
		  "HEX0769C13",
		  "HEX0769D05",
		  "HEX0769G19",
		  "HEX0769H06",
		  "HEX0769H11",
		  "HEX0769O12",
		  "HEX0769P14",
		  "HEX0769P20",
		  "HEX0788J20",
		  "HEX0865B02",
		  "HEX0865B23",
		  "HEX0865C12",
		  "HEX0865C20",
		  "HEX0865P15",
		  #-----------------------------+
		  # CONVERTED
		  # 08/22/2007
		  #-----------------------------+
		  #"HEX0865E02",
		  #"HEX0865E14",
		  #"HEX0865G02",
		  #"HEX0865G11",
		  #"HEX0865H06",
		  #"HEX0865I08",
		  #"HEX0865I12",
		  #"HEX0865J23",
		  #"HEX0865K11",
		  #"HEX0865K15",
		  #"HEX0865L17",
		  #"HEX0865M02",
		  #"HEX0865M13",
		  #"HEX0865N22",
		  #"HEX0865P01",
		  #"HEX0865P10",
		  #"HEX0905N22",
		  #"HEX0915K16",
		  #"HEX0922O11",
		  #"HEX0961B19",
		  #"HEX0961C24",
		  #"HEX0961D14",
		  #"HEX0961E20",
		  #"HEX0961F06",
		  #-----------------------------+
		  # CONVERTED
		  # 08/09/2007
		  #-----------------------------+
		  #"HEX0961F07",
		  #"HEX0961F22",
		  #"HEX0961G12",
		  #"HEX0961H02",
		  #"HEX0961H13",
		  #"HEX0961H19",
		  #"HEX0961I05",
		  #"HEX0961K08",
		  #"HEX0961L11",
		  #"HEX0961L14",
		  #"HEX0961L18",
		  #"HEX0961L22",
		  #"HEX0961N09",
		  #"HEX0961N16",
		  #"HEX0961N22",
		  #"HEX0961O15",
		  #"HEX0987A24",
		  #"HEX1011A11",
		  #"HEX1057B13",
		  #"HEX1057C05",
		  #"HEX1057D07",
		  #"HEX1057D08",
		  #"HEX1057D14",
		  #"HEX1057D20",
		  #"HEX1057E09",
		  #-----------------------------+
		  # CONVERTED
		  # 08/09/2007
		  #-----------------------------+
		  #"HEX1057E15",
		  #"HEX1057E21",
		  #"HEX1057F03",
		  #"HEX1057F11",
		  #"HEX1057F21",
		  #"HEX1057G04",
		  #"HEX1057G13",
		  #"HEX1057H15",
		  #"HEX1057H16",
		  #"HEX1057H22",
		  #"HEX1057I17",
		  #"HEX1115D24",
		  #"HEX1118O08",
		  #"HEX1184K09",
		  #"HEX1198N12",
		  #"HEX1224I17",
		  #"HEX1311B22",
		  #"HEX1332H23",
		  #"HEX1332N20",
		  #"HEX1549K12",
		  #"HEX1564B16",
		  #"HEX1745I24",
		  #"HEX1773P12",
		  #"HEX1804D04",
		  #"HEX1819O17",
		  #"HEX1866D14",
		  #"HEX2029F09",
		  #"HEX2358G08",
		  #"HEX2385J18",
		  #"HEX2436E11",
		  #"HEX2439N11",
		  #"HEX2444K23",
		  #"HEX2517B24",
		  #"HEX2520B04",
		  #"HEX2584A05",
		  #"HEX2874I04",
		  #"HEX2884C05",
		  #"HEX2903P03",
		  #"HEX2986I03",
		  #"HEX3045G05",
		);

my $num_files_to_convert = @root_names;

print "$num_files_to_convert files to convert\n";

for my $bat_root_name (@root_names) {

    my $game_created = $bat_root_name."/".$bat_root_name.".game.xml"; 
    my $game_new_copy = $bat_root_name.".game.xml";
    my $gff_src = "$bat_root_name/gff/$bat_root_name.gff";

    # Convert gff file to game
    if (-e $gff_src) {
	my $cnv_cmd = "cnv_gff2game.pl".
	    " -g $bat_root_name/gff/$bat_root_name.gff".
	    " --ap-path ".$ap_path.
	    " -i $bat_root_name/rm/".$bat_root_name."_TREP9.masked.fasta".
	    " -o $bat_root_name/$bat_root_name.game.xml";
	system ($cnv_cmd);
	#print "$cnv_cmd\n";
    }
    else {
	print "Source file not found:\n\t$gff_src";
    }

    # Copy game.xml file to main directory
    # This will facilitate the copy command to move these to jlb3
    #my $cp_cmd = ""

    if (-e $game_created) {
	copy ($game_created, $game_new_copy);
    } # End of if file exists
    
} # End of for each base ind_root_name in the root names list

exit;

#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+



