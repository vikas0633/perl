#!/usr/bin/perl
#----------------------------------------------------------#
#        Author: Douglas Senalik dsenalik@wisc.edu         #
#----------------------------------------------------------#
# "Black Box" program series
=bb
Find open reading frames in a DNA or RNA sequence
=cut bb
use strict;
use warnings;
use Getopt::Long;      # for getting command line parameters
# 1.0.0 - Mar 9, 2011
# 1.1.0 - June 8, 2011 = Add --trimheader, fix major bug in sorting alpha instead of numeric



############################################################
# configuration variables
############################################################
my @startcodons = ( "ATG" );
my @fullstartcodons = ( "ATG", "GTG", "CTG", "TTG" );
my @stopcodons = ( "[TU]AA", "[TU]AG", "[TU]GA" );
# thresholds for guessing orientation
my $guessthreshold = 0.10;  # "ambiguous" if both directions
my $guessminbp = 100;  # "tooshort" if below this
my %guessstats = ();



############################################################
# global variables
############################################################
my $ansiup    = "\033[1A";  # terminal control
(my $prognopath = $0) =~ s/^.*[\/\\]//;
my $currhdr = "";
my $currseq = "";
my $guessdata = "";
my $numseq = 0;
my %orfperid = ();
my $outfilename2 = "";  # --nonorffasta file name



############################################################
# command line parameters
############################################################
my $infilename  = "";   # input file name
my $outfilename = "";   # output file name
my $anystart    = 0;    # start also at beginning of sequence
my $fullstart   = 0;    # use all 4 start codons
my $minlen      = 100;  # minimum orf length in b.p.
my $guessorient = 0;    # different output format: guess orientation
my $fasta       = 0;    # different output format: FASTA
my $nonorffasta = 0;    # second file with all sequence not in the first one
my $fastacollapse = 0;  # no duplication in --fasta output
my $fastalargest = 0;   # only larger of any overlapping orfs is kept
my $trimheader  = 0;    # keep header only up to first white space
my $help        = 0;    # print help and exit
my $quiet       = 0;    # only show errors
my $debug       = 0;    # print extra debugging information
GetOptions (
            "infile=s"         => \$infilename,        # string
            "outfile=s"        => \$outfilename,       # string
            "fullstart"        => \$fullstart,         # flag
            "anystart"         => \$anystart,          # flag
            "minlen=i"         => \$minlen,            # integer
            "guessorientation" => \$guessorient,       # flag
            "fasta"            => \$fasta,             # flag
            "nonorffasta"      => \$nonorffasta,       # flag
            "fastacollapse"    => \$fastacollapse,     # flag
            "fastalargest"     => \$fastalargest,      # flag
            "trimheader"       => \$trimheader,        # flag
            "help"             => \$help,              # flag
            "quiet"            => \$quiet,             # flag
            "debug"            => \$debug);            # flag
# required parameters
unless ( $infilename ) { $help = 1; }
unless ( $outfilename ) { $help = 1; }
# debug implies not quiet
if ( $debug ) { $quiet = 0; }



############################################################
# print help screen
############################################################
if ( $help )
  {
    print "$prognopath

This program will detect open reading frames in FASTA
DNA or RNA sequences. This is similar to the NCBI program at 
http://www.ncbi.nlm.nih.gov/gorf/orfig.cgi

Required parameters:
  --infile=xxx       input file name
  --outfile=xxx      output file name, use \"-\" for stdout
Optional parameters:
  --fullstart        use full set of start codons: ATG GTG CTG TTG
                     the default is to only use ATG
  --anystart         start of sequence is also a valid orf start
  --minlen=xxx       minimum orf length in b.p., default=$minlen
  --guessorientation guess orientation based on strand with the
                     most total orfs, this data will be be saved
                     instead of the list of orfs
  --fasta            output file is in FASTA format, each orf is
                     a separate sequence, information is in header
  --nonorffasta      second FASTA file with all sequence not in
                     the first one. File name is --outfile name
                     with \"nonorf\" inserted
  --fastacollapse    combine overlapping sequence in the FASTA file
  --fastalargest     if two orfs overlap, keep only the larger one
  --trimheader       remove any text in the FASTA header after 
                     the first occurrence of white space
  --help             print this screen
  --quiet            only print error messages
  --debug            print extra debugging information
";
    exit 1;
  } # if ( $help )



############################################################
# perform orf detection
############################################################
my $OUTF = stdopen ( ">", $outfilename );
my $OUTF2;
if ( $nonorffasta )
  {
    ( $outfilename2 = $outfilename ) =~ s/(\.[^\.]*)$/.nonorf$1/;
    $OUTF2 = stdopen ( ">", $outfilename2 );
  }
unless ( ( $guessorient ) or ( $fasta ) )
  { print $OUTF "#ID\torf\tFrame\tStart\tCodon\tStop\tCodon\tLength\n"; }

my $INF = stdopen ( "<", $infilename );
while ( my $aline = <$INF> )
  {
    $aline =~ s/[\r\n]//g;
    if ( $aline =~ m/^>(.*)$/ )
      {
        my $seqhdr = $1;
        if ( $trimheader ) { $seqhdr =~ s/\s.*$//; }  # Trim anything after first white space
        unless ( $seqhdr ) { die "Error, FASTA header is missing, or starts with white space: \"$aline\"\n"; }
        process();
        $currhdr = $seqhdr;
        $currseq = "";
      }
    else
      {
        $currseq .= $aline;
      }
  } # while <$INF>
stdclose ( $INF );
process();

if ( $guessorient )
  {
    print $OUTF "#ID\tSeq Len\tbp Fwd\tbp Rev\tOrientation\n", $guessdata;
    print "Guess Orientation Summary:\n";
    while ( my ( $key, $value ) = each ( %guessstats ) )
      { print commify($value), "\t", $key, "\n"; }
  }

unless ( $quiet )
  { print commify($numseq), " sequences processed\n"; }

stdclose ( $OUTF );
if ( $nonorffasta )
  { stdclose ( $OUTF2 ); }



############################################################
sub process {
############################################################
  if ( $currseq )  # first call will have no sequence
    {
      $numseq++;
      # display progress for large input files
      if ( ( ( $numseq % 10 ) == 0 ) and ( ! $quiet ) )
        { print commify($numseq), "\n", $ansiup; }

      my @outdata = ();
      my $currseqrc = revcomp($currseq);
      my @start = ();
      my @startrc = ();
      my @stop = ();
      my @stoprc = ();
      my $currlen = length($currseq);
      # variables used for guessing the orientation
      my $bpplus = 0;
      my $bpminus = 0;

      # only use stop codon once, later hits will be for shorter orf
      my %stopused = ();
      my %stopusedrc = ();

      # to avoid column alignment problems, make sure header has no tabs in it
      $currhdr =~ s/\t/ /g;

      debugmsg ( "Processing \"$currhdr\", " . commify($currlen) . " b.p." );

      # if --anystart, add start of sequence positions in all reading frames
      if ( $anystart )
        {
          push ( @start, 0, 1, 2 );
          push ( @startrc, 0, 1, 2 );
        }

      foreach my $acodon ( @startcodons )
        {
          while ( $currseq =~ m/($acodon)/ig )
            {
              debugmsg ( "Fwd: Start codon $1 at ".(pos($currseq)-2) );
              push ( @start, pos($currseq)-2 );
            }
          while ( $currseqrc =~ m/($acodon)/ig )
            {
              debugmsg ( "Rev: Start codon $1 at ".(pos($currseqrc)-2) );
              push ( @startrc, pos($currseqrc)-2 );
            }
        }
      foreach my $acodon ( @stopcodons )
        {
          while ( $currseq =~ m/($acodon)/ig )
            {
              debugmsg ( "Fwd: Stop codon $1 at ".(pos($currseq)-2) );
              push ( @stop, pos($currseq)-2 );
            }
          while ( $currseqrc =~ m/($acodon)/ig )
            {
              debugmsg ( "Rev: Stop codon $1 at ".(pos($currseqrc)-2) );
              push ( @stoprc, pos($currseqrc)-2 );
            }
        }

      # sort arrays to put positions in increasing order
      @start = sort ( { $a <=> $b } @start );
      @stop = sort ( { $a <=> $b } @stop );
      @startrc = sort ( { $a <=> $b } @startrc );
      @stoprc = sort ( { $a <=> $b } @stoprc );

      # find start-stop pairs
      foreach my $astart ( @start )
        {
          my $frame = ( $astart % 3 );
          foreach my $astop ( @stop )
            {
              # frame must match
              if ( ( $astop % 3 ) != $frame )
                {
                  debugmsg ( "Fwd, Start=$astart, Stop=$astop, Frame mismatch $frame vs. ". ( $astop % 3 ) );
                  next;
                }

              # stop must be > $start
              if ( $astop <= $astart )
                {
                  debugmsg ( "Fwd, Start=$astart, Stop=$astop, Stop $astop < Start $astart" );
                  next;
                }

              # ignore this start codon if we hit a used stop codons ( meaning a start codon existed earlier upstream )
              if ( $stopused{$astop} )
                {
                  debugmsg ( "Fwd, Start=$astart, Stop=$astop, stop codon already used" );
                  last;
                }

              # length must be >= --minlen parameter
              my $orflen = ( $astop - $astart + 3 );
              if ( $orflen >= $minlen )
                {
                  # [0]ID  [1]orf  [2]Frame  [3]Start  [4]Codon  [5]Stop  [6]Codon  [7]Length
                  $orfperid{$currhdr}++;
                  my @outline = ( $currhdr,
                                  $orfperid{$currhdr},
                                  "+".($frame?$frame:3),
                                  $astart,
                                  substr($currseq,$astart-1,3),
                                  ($astop+2),
                                  substr($currseq,$astop-1,3),
                                  $orflen );
                  push ( @outdata, \@outline );
                  $bpplus += $orflen;
                  debugmsg ( "Fwd, Valid orf at $astart to $astop, length=$orflen" );
                }
              else
                { debugmsg ( "Fwd, Discarding too-short orf ( $orflen b.p. ) at $astart to $astop" ); }

              # if a stop codon is used, don't use it again as that will be for
              # a shorter orf inside a larger one. So, set flag to ignore it
              $stopused{$astop} = 1;

              # exit inner loop when first orf found, even if too short
              last;
            }
        } # foreach my $astart ( @start )

      foreach my $astart ( @startrc )
        {
          my $frame = ( $astart % 3 );
          foreach my $astop ( @stoprc )
            {
              # frame must match
              if ( ( $astop % 3 ) != $frame )
                {
                  debugmsg ( "Rev, Start=$astart, Stop=$astop, Frame mismatch $frame vs. ". ( $astop % 3 ) );
                  next;
                }

              # stop must be > $start
              if ( $astop <= $astart )
                {
                  debugmsg ( "Rev, Start=$astart, Stop=$astop, Stop $astop < Start $astart" );
                  next;
                }

              # ignore this start codon if we hit a used stop codons ( meaning a start codon existed earlier upstream )
              if ( $stopusedrc{$astop} )
                {
                  debugmsg ( "Rev, Start=$astart, Stop=$astop, stop codon already used" );
                  last;
                }

              # length must be >= --minlen parameter
              my $orflen = ( $astop - $astart + 3 );
              if ( $orflen >= $minlen )
                {
                  # [0]ID  [1]orf  [2]Frame  [3]Start  [4]Codon  [5]Stop  [6]Codon  [7]Length
                  $orfperid{$currhdr}++;
                  my @outline = ( $currhdr,
                                  $orfperid{$currhdr},
                                  "-".($frame?$frame:3),
                                  ($currlen-$astart+1),
                                  substr($currseqrc,$astart-1,3),
                                  ($currlen-($astop+2)+1),
                                  substr($currseqrc,$astop-1,3),
                                  $orflen );
                  push ( @outdata, \@outline );
                  $bpminus += $orflen;
                  debugmsg ( "Rev, Valid orf at $astart to $astop, length=$orflen" );
                }
              else
                { debugmsg ( "Rev, Discarding too-short orf ( $orflen b.p. ) at $astart to $astop" ); }

              # if a stop codon is used, don't use it again as that will be for
              # a shorter orf inside a larger one. So, set flag to ignore it
              $stopusedrc{$astop} = 1;

              # exit inner loop when first orf found, even if too short
              last;
            }
        } # foreach my $astart ( @startrc )

      # save output data to file
      # sort output data by orf length, largest first
      @outdata = sort { $b->[7] <=> $a->[7] } @outdata;

      # save guess orientation data
      if ( $guessorient )
        {
          my $txt;
          if ( ( $bpplus + $bpminus ) < 1 )
            { $txt = "no orfs"; }
          elsif ( ( $bpplus + $bpminus ) < $guessminbp )
            { $txt = "orfs too short"; }
          elsif ( ( abs( $bpplus - $bpminus ) / ( $bpplus + $bpminus ) ) < $guessthreshold )
            { $txt = "ambiguous"; }
          elsif ( $bpplus > $bpminus )
            { $txt = "forward"; }
          else
            { $txt = "reverse"; }

          $guessdata .= join ( "\t", $currhdr, $currlen, $bpplus, $bpminus, $txt ) . "\n";
          $guessstats{$txt}++;
        } # if ( $guessorient )

      elsif ( $fasta )
        {
          my %used = ();  # for --fastacollapse or --fastalargest only, key is start, value is end
          foreach my $rowref ( @outdata )
            {
              # [0]ID  [1]orf  [2]Frame  [3]Start  [4]Codon  [5]Stop  [6]Codon  [7]Length
              my $header = ">" . $rowref->[0].".".$rowref->[1] . " " .  join ( "; ", @{$rowref}[2..7] );
              my $len = $rowref->[7];
              my $seq = "error";
              if ( ( $fastacollapse ) or ( $fastalargest ) )
                {
                  # get range, start always the lower value
                  my $start = $rowref->[3];
                  my $end = $rowref->[5];
                  if ( $end < $start )
                    {
                      $start = $rowref->[5];
                      $end = $rowref->[3];
                    }

                  # see if this range overlaps an existing one
                  my $inserted = 0;
                  foreach my $key ( keys %used )
                    {
                      if ( ( $start >= $key ) and ( $end <= $used{$key} ) )
                        {
                          # full subset, can ignore it completely for both --fastacollapse and --fastalargest
                          debugmsg ( "Sequence \"$currhdr\": orf $start..$end is subset of orf $key..$used{$key}" );
                          $inserted = 1;
                          last;
                        } # if ( ( $start >= $key ) and ( $end <= $used{$key} ) )
                      elsif ( ( $end > $used{$key} ) and ( $start <= $used{$key} ) )
                        {
                          if ( $fastacollapse )
                            {
                              debugmsg ( "Sequence \"$currhdr\": orf $start..$end extends end of orf $key..$used{$key}" );
                              # extend end side
                              $used{$key} = $end;
                              # modify start to detect join events
                              if ( $start > $key ) { $start = $key }
                            } # if ( $fastacollapse )
                          if ( $fastalargest )
                            {
                              if ( ( $end - $start ) > ( $used{$key} - $key ) )  # if new orf larger
                                {
                                  debugmsg ( "Sequence \"$currhdr\": orf $start..$end larger than orf $key..$used{$key}" );
                                  $used{$start} = $end;
                                  delete ( $used{$key} );
                                }
                            } # if ( $fastalargest )
                          $inserted = 1;
                        } # elsif ( ( $end > $used{$key} ) and ( $start <= $used{$key} ) )
                      elsif ( ( $start < $key ) and ( $end >= $key ) )
                        {
                          if ( $fastacollapse )
                            {
                              debugmsg ( "Sequence \"$currhdr\": orf $start..$end extends start of orf $key..$used{$key}" );
                              # extend start side
                              $used{$start} = $used{$key};
                              delete ( $used{$key} );
                              # modify end to detect join events
                              if ( $end < $used{$start} ) { $end = $used{$start} }
                            } # if ( $fastacollapse )
                          if ( $fastalargest )
                            {
                              if ( ( $end - $start ) > ( $used{$key} - $key ) )  # if new orf larger
                                {
                                  debugmsg ( "Sequence \"$currhdr\": orf $start..$end larger than orf $key..$used{$key}" );
                                  $used{$start} = $end;
                                  delete ( $used{$key} );
                                }
                            } # if ( $fastalargest )
                          $inserted = 1;
                        } # elsif ( ( $start < $key ) and ( $end >= $key ) )
                    } # foreach

                  # no action taken above, so new non-overlapping range found
                  unless ( $inserted )
                    {
                      $used{$start} = $end; 
                      debugmsg ( "Sequence \"$currhdr\": orf $start..$end is new" );
                    } # unless ( $inserted )

                } # if ( ( $fastacollapse ) or ( $fastalargest ) )
              else  # do not collapse overlaps
                {
                  if ( $rowref->[2] =~ m/\-/ ) # if reverse complement
                    {
                      $seq = substr( $currseq, $rowref->[5]-1, $len );
                      $seq = revcomp ( $seq );
                      $header .= " [RC]";
                    }
                  else  # forward orientation
                    {
                      $seq = substr( $currseq, $rowref->[3]-1, $len );
                    }
                  print $OUTF $header, "\n";
                  print $OUTF $seq, "\n";
                } # else do not collapse overlaps
            } # foreach my $rowref ( @outdata )

          if ( ( $fastacollapse ) or ( $fastalargest ) )
            {
              debugmsg ( "Sequence \"$currhdr\" has ".scalar ( keys %used )." orf sequences to save" );
              my $numseq = 0;
              my $prevend = 0;
              foreach my $key ( sort { $a <=> $b } keys %used )
                {
                  $numseq++;
                  my $len = ( $used{$key} - $key + 1 );
                  my $header = ">$currhdr.$numseq $key..$used{$key} = $len b.p.";
                  my $seq = substr( $currseq, $key-1, $len );
                  debugmsg ( "Save FASTA of orf $key..$used{$key} length $len" );
                  print $OUTF $header, "\n";
                  print $OUTF $seq, "\n";
                  if ( $nonorffasta )
                    {
                      my $nolen = ( $key - $prevend - 1 );
                      my $noseq = substr( $currseq, $prevend, $nolen );
                      my $noheader = ">$currhdr.$numseq ".($prevend+1)."..".($key-1)." = $nolen b.p.";
                      debugmsg ( "Save FASTA of non-orf ".($prevend+1)."..".($key-1)." length $nolen" );
                      $prevend = $used{$key};
                      print $OUTF2 $noheader, "\n";
                      print $OUTF2 $noseq, "\n";
                    } # if ( $nonorffasta )
                } # foreach
              if ( $nonorffasta )
                {
                  $numseq++;
                  my $nolen = ( length($currseq) - $prevend );
                  my $noseq = substr( $currseq, $prevend, $nolen );
                  my $noheader = ">$currhdr.$numseq ".($prevend+1)."..".length($currseq)." = $nolen b.p.";
                  debugmsg ( "Save FASTA of last non-orf ".($prevend+1)."..".length($currseq)." length $nolen" );
                  print $OUTF2 $noheader, "\n";
                  print $OUTF2 $noseq, "\n";
                } # if ( $nonorffasta )
            } # if ( ( $fastacollapse ) or ( $fastalargest ) )

        } # if ( $fasta )

      else  # normal output
        {
          foreach my $rowref ( @outdata )
            { print $OUTF join ( "\t", @{$rowref} ), "\n"; }
        }

      @outdata = ();
    }
} # sub process



############################################################
sub debugmsg { my ( $text, $noreturn, $nolinenum ) = @_;
############################################################
  if ( $debug )
    {
      my ($package, $filename, $line, $sub) = caller(0);
      unless ( $nolinenum ) { $text = "Line $line: " . $text; }
      if ( ! ( $noreturn ) ) { $text .= "\n"; }
      print $text;
    } # if ( $debug )
} # sub debugmsg



###############################################################
sub timestr {
###############################################################
  @_ = localtime(shift || time);
  return(sprintf("%04d/%02d/%02d %02d:%02d", $_[5]+1900, $_[4]+1, $_[3], @_[2,1]));
} # sub timestr



###############################################################
sub commify {
###############################################################
# http://perldoc.perl.org/perlfaq5.html#How-can-I-output-my-numbers-with-commas
  local $_ = shift;
  1 while s/^([-+]?\d+)(\d{3})/$1,$2/;
  return $_;
} # commify



############################################################
sub revcomp { my ( $dna ) = @_;
############################################################
# standard DNA reverse complement, including degenerate bases
  my $revcomp = reverse ( $dna );
  $revcomp =~ tr/AaCcTtGgMmRrYyKkVvHhDdBb/TtGgAaCcKkYyRrMmBbDdHhVv/;
  return $revcomp;
} # sub revcomp



###############################################################
sub stdopen { my ( $mode, $filename, $extratext ) = @_;
###############################################################
# a replacement for the three-parameter open which also allows
# the use of "-" as the file name to mean STDIN or STDOUT
  my $fh;  # the file handle
  if ( $filename eq "-" )  # only exact match to "-" has special meaning
    {
      if ( $mode =~ m/>/ )
        { $fh = *STDOUT }
      else
	{ $fh = *STDIN }
    }
  else
    {
      # supplemental passed text for error messages, need one more space
      if ( defined $extratext )
        { $extratext .= " " }
      else
	{ $extratext = "" }

      my $text;  # this is only used for error message
      if ( $mode =~ m/^\+?>>/ )  # ">>" or "+>>"
        { $text = "append" }
      elsif ( $mode =~ m/^\+?>/ )  # ">" or "+>"
        { $text = "output" }
      elsif ( $mode =~ m/^\+?</ )  # "<" or "+<"
        { $text = "input" }
      else
	{ die "Error, unsupported file mode \"$mode\" specified to stdopen( $mode, $filename, $extratext )\n"; }
      open ( $fh, $mode, $filename ) or die ( "Error opening ${extratext}file \"$filename\" for $text: $!\n" );
    }
  # return the opened file handle to the caller
  return $fh;
} # stdopen



###############################################################
sub stdclose { my ( $fh ) = @_;
###############################################################
# same as built-in close, except in case of STDIN or STDOUT,
# and in this case the file handle is not closed

  unless ( fileno($fh) <= 2 )  # if file number is this low, is stdin or stdout or stderr
    { close ( $fh ) or die ( "Error closing file handle: $!\n" ); }

} # sub stdclose



# eof

=test_data
>A1Contig1
TCACACTTTTAGCTAACAGATTTACATATCTCCGTTGGAGATGCTCTAACATGACTAGAATGAGATGGCCCTCGTTAAAA
TTTTAAAATCTCGGCAGAGAAGTTGTTAAAAAGAGCACAATTCATGAGTTCGATAACTGATGAGAAACTGAGTGAACATA
AAGCAACACAGAGGGCAGTACAACTCTTATTTAAGAACGAGATTAGAAATGTCCAACATCTTCACAAATGAAGGACCGAA
GTTTGATATGGAATTCTGACTTGATAGAGACATTAGAGTTGGAAAACTTGTTGATCAATGCAGCTATTACCATGCATTCA
GCTGAAGCAAGAAAAGAGAGCAGAGGAGCTCATGCTCGTGAAGATTTTACGAAAAGAGATGATGAGAATTGGATGAAGCA
TTCATTGGGATACTGGGAGAATGAGAAGGTACGGCTGGACTACAGGCCTGTTCACATGAACACATTAGACGATGAGATTG
AAACGTTCCCACCTAAAGCACGTGTGTACTGAGGTATGGTATTGAAGGATACATGTGGTTGGGGAAAAATATAATATTTT
CTTCTAAGAAGTCCGCAATAAATTTACTGAGACCTAGAAAATTTCTAAATAATGACGTCATTTGTCAAACTGCAGTCGGG
GTGAGTAAGCTTGTCTCTGTCATGAGCAAGGGTGTGAGAATACCTCACTTGTATTGATCATGTTCCTGGCAACAAATATT
TTATACTACTCATAAAAATGACTCCAAACTTTCAATTACTAGTAACAACTGGATTGGAAGTTAAATCAAAGTGTTATCTT
ATTATTGGTAGCTCAGTTGAACTTTGTTTTAGAACAAAATGCTTACTATAATTTCGTGAGCATTATGACTTTGGATGACA
TTTTGGTCATTGTATTTGTGTATCAAATGTGAGATTAGTTGCTGTCTACT

NCBI results from http://www.ncbi.nlm.nih.gov/gorf/orfig.cgi
Frame		from		to	Length
 +3		228	..	512	285
 +2		533	..	697	165
 -1		493	..	621	129
 +3		738	..	851	114
 -3		203	..	313	111
 +1		124	..	231	108

Results from this program
#ID     orf     Frame   Start   Codon   Stop    Codon   Length
A1Contig1       2       +3      228     ATG     512     TGA     285
A1Contig1       3       +2      533     ATG     697     TGA     165
A1Contig1       5       -1      621     ATG     493     TAG     129
A1Contig1       4       +3      738     ATG     851     TAA     114
A1Contig1       6       -3      313     ATG     203     TAA     111
A1Contig1       1       +1      124     ATG     231     TGA     108

=cut test_data
