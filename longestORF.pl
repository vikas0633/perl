#!/usr/bin/perl

use Getopt::Long;
use Bio::SeqIO;

my ($sequencefile,$sequenceformat,$help) = (undef, 'fasta',undef);

&GetOptions('input|i=s'              => \$sequencefile,
            'format|f=s'             => \$sequenceformat,
            'help|h'                 => \$help,
            );

if ($help) {
   exec('perldoc', $0); # when I can be bothered
   die;
}

sub longestORF {
   my $best=0;
   my ($bests, $beste, $beststrand) = (-1,-1,0);

   my $dna=Bio::Seq->new(-seq => $_[0]);
   my %strand=('+'=>$dna->seq,
               '-'=>$dna->revcom->seq);

   foreach $direction (keys %strand) {
      my @starts;
      my @stops;
      while ($strand{$direction}=~m/(atg|^!taa|^!tga|^!tag)/gi) {
         push @starts,pos($strand{$direction})-2;
      }

      while ($strand{$direction}=~m/(taa|tga|tag|...$)/gi) {
         push @ends,pos($strand{$direction})-2;
      }

      for my $s (@starts) {
         for my $e (@ends) {
            if ($e%3==$s%3 and $e>$s) {
               if ($e-$s>$best) {
                  $best=$e-$s;
                  ($bests,$beste,$beststrand) = ($s, $e, $direction);
                  $bestorf=Bio::Seq->new(-seq=>$strand{$direction})->subseq($s,$e);
               }
            last
            } else {
               next
            }
         }
      }
      @starts=0;
      @ends=0;
   }
   return ($best, $bests, $beste, $beststrand, $bestorf);
}


my $seqio = new Bio::SeqIO('-format' => $sequenceformat,
                           '-file'   => $sequencefile );

while (my $dna = $seqio->next_seq) {

   print ">",$dna->display_id," ",$dna->desc,": ";
   $seq=$dna->seq;
   ($length,$start,$end,$direction,$sequence)=longestORF($seq);
   print $length,"\n",Bio::Seq->new(-seq=>$sequence)->translate->seq,"\n";
}
