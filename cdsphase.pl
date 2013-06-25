#!/usr/bin/perl

=head1 cdsphase.pl

read GFF, genome fasta, check CDS_exon features for proper phase; add if '.', and change if wrong?
d.gilbert, aug2006, gilbertd@indiana.edu

=head1 usage

 perl cdsphase.pl $sc/caf1a/dgil/dyak_caf1_DGIL_SNO.gff.gz $sc/dyak3/perchr/ \
  > dyak_caf1_DGIL_SNO.gffp

 options: -v verbose ; -tCDS restrict to this feature (CDS,exon default)
  could add option to output aa translations, cds dna

  Uses dna fasta file per scaffold/chromosome as $dnapath/$ref.fa

=head1 Dros. species gene predictions

-- add phase column 8 to SNAP predictions GFF
./cdsphase.pl $sc/caf1a/dgil/caf1a/dpse_caf1_DGIL_SNO.gff.gz $sc/dpse2/perchr/ > dpse_caf1_DGIL_SNO.gffp

-- check other predictions, dmoj case:
./cdsphase.pl -v -tCDS $sc/caf1a/bren/dmoj_caf1.gff3.gz $sc/dmoj3/perchr/
  -- mostly good, cdsphase changes some phases which may be due to different scoring of translations

./cdsphase.pl -v -tCDS $sc/caf1a/batz/contrast_na/dmoj.gff.gz $sc/dmoj3/perchr/
  -- partial genes have phase errors, needs corrected phases

./cdsphase.pl -v -tCDS $sc/caf1a/rgui/geneidv1.2/dmoj_caf1.gff3.gz $sc/dmoj3/perchr/ 
 -- this one has good phases, and has note on 5prime_partial, 3prime_partial 
   (but also has lots of (partial) predictions of a few aminos .. should trim all short predictions <~ 10-15 aa

./cdsphase.pl -v -tCDS $sc/caf1a/ncbi/dmoj_caf1_NCBI_GNO.gff.gz $sc/dmoj3/perchr/ 
 -- this one has wrong phase values, needs correction

=cut

use strict;
use lib("/bio/argos/common/perl/lib/");
use Bio::SeqIO;

my $verbose= 0;
my $outh = *STDOUT;
my $gffin= *STDIN;

# each predictor uses diff feature types; do we also need mRNA,gene parent?
my @cdsTypes= qw(CDS exon); # FIXME
my %cdsTypes= map { $_ => 1; } @cdsTypes;
my %seenType=();

# expect args:  cdsphase.pl  data.gff[.gz] dnafastapath/
my %args=();
my $gff= shift(@ARGV);
while($gff =~ s/^-(.)//){
  $args{$1}= $gff || 1;
  $gff= shift(@ARGV);
  }
my $genomefasta= shift(@ARGV);
unless( -f $gff && -d $genomefasta) {
  die "Usage:  cdsphase.pl [-v] [-tCDS] annot.gff[.gz]  dnafastapath/\n" ;
  }
my $ok;
if($gff =~ /\.gz$/) { $ok= open(GFF, "gunzip -c $gff|"); $gffin = *GFF; }
else { $ok= open(GFF,$gff); $gffin = *GFF; }

if($args{v}){ $verbose= $args{v}; }
if($args{t}){ @cdsTypes=split(/[,;.-]/,$args{t}); }

processGff($gffin,$genomefasta,$outh);
close($gffin);
close($outh);

#--------------------

sub processGff {
  my($gffHandle,$orggenomefasta,$outh)= @_;

  %seenType=();
  %cdsTypes= map { $_ => 1; } @cdsTypes;
  
  warn "# cdsphase: $orggenomefasta\n" if $verbose;
  
  my($faIO, $refdna)= (undef, undef);
  my($l_id, $l_ref, $skipref, @feats,@cds);

  # any patch for $ref to gff-ref ?  e.g. dpse 'Ch'; dper/dsec super/scaffold; ...
  ## assume cds grouped by ID in gff

  while(my $gffin= <$gffHandle>) {
    if($gffin =~ m/^\W/) { print $outh $gffin; next; }
    
    ## need to collect all CDS_exons/gene-mRNA first ...
    my @gff= split "\t",$gffin,9;
    my $ref = $gff[0];
    my $type= $gff[2];
    chomp($gff[8]);
    my $id= $gff[8];  $id =~ s/^(ID|Parent)=//; $id =~ s/[;,\s].*$//;
    $seenType{$type}++;
    # check here for exon AND CDS : drop exon if have CDS
    %cdsTypes=("CDS" => 1) if( $type eq "CDS" ); #$cdsTypes{"exon"}= 0 

    if($id ne $l_id && @cds) {
      my $cdsnew= cdsPhase( $refdna, \@cds); # sorted @cds / gene; not all @feats
      foreach my $gff (@$cdsnew) { print $outh join("\t",@$gff),"\n"; }
      @cds=(); # @feats=(); 
      }

    if($ref ne $l_ref) {
      undef $refdna; undef $faIO;
      $skipref=0;
      my $fafile= "$orggenomefasta/$ref.fa";
      # check here for $ref.fa variations: dpse: Ch2..ChX..
      unless(-e $fafile) {
        my $fref= $ref;
        if($fref =~ /^Ch/) { $fref =~ s/^Ch//; }
        elsif($fref =~ /^super/) { $fref =~ s/^super/scaffold/; }
        $fafile= "$orggenomefasta/$fref.fa"; 
      }
      unless(-e $fafile) {
        warn "*** Cant find dna.fa for $ref at $fafile\n#*** skipping cdsPhase edits for $ref\n";
        $skipref=1;  
      } else {
        $faIO = Bio::SeqIO->new(-file => $fafile,-format => 'Fasta');
        $refdna = $faIO->next_seq(); # wait till have gff ?
      }
      # warn "# cdsphase[$ref]: $fafile\n" if $verbose;
      $l_ref= $ref;
      }
      
    if( $cdsTypes{$type} && !$skipref) { push(@cds,\@gff); } 
    else { print $outh $gffin; }   #? dont need to keep
    
    $l_id= $id; $l_ref= $ref;
    }
    
  if(@cds) {
    my $cdsnew= cdsPhase( $refdna, \@cds); # sorted @cds / gene; not all @feats
    foreach my $gff (@$cdsnew) { print $outh join("\t",@$gff),"\n"; }
    @cds=();  
    }
    
  undef $refdna; undef $faIO;
}



sub cdsPhase {
  my($refdna, $cdsA)= @_;
  ## assume cdsA are all/only cds exon set for one gene/mrna
  return $cdsA unless(ref $refdna);
  
  my $cstrand= $cdsA->[0]->[6];
  my $isrev= ($cstrand eq '-' || $cstrand < 0);
  
  my @cds;   # sort by start
  if ($isrev) { @cds= sort{ $b->[3] <=> $a->[3] } @$cdsA; } # end 1st
  else { @cds= sort{ $a->[3] <=> $b->[3] } @$cdsA; } # start 1st
  my $nt_length= 0;
  my $ispartial= 0;
  
  foreach my $ix (0 .. $#cds)  {
    my($ref,$src,$type,$start,$stop,$score,$strand,$phase,$attr)= @{$cds[$ix]};
    my $id=$attr; $id =~ s/^(ID|Parent)=//; $id =~ s/[;,\s].*$//;
    # do we check exon ordering? use as given in gff?  need 1st .. last, differs for strands
    
    if($ix == 0) { 
      ## 1st exon; find start ATG; ** only need 3 bases at start, not all
      $ispartial= 0;
      my($bstart,$blen)= ($start - 1, $stop - $start + 1);
      #my($bstart,$blen)= ($isrev) ? ($stop-8,8) : ($start-1, 8);
      my $exondna  = substr( $refdna->seq(), $bstart, $blen);
      if($isrev) {  
        $exondna = reverse $exondna;
        $exondna =~ tr/gatcGATC/ctagCTAG/;
        }
      my $inc5= 0;
      for (; $inc5<=3; $inc5++) {
        my $atg= substr($exondna, $inc5, 3);
        last if($atg =~ /atg/i);
        }

      ## fixme, if $ispartial probably need check best aa translation frame
      ## yes; need full cds/all exons and translate() method      
      if ($inc5 > 2) { 
        $nt_length = 0;  $ispartial=1; $inc5 = 0; #start not found/incomplete prot ?
        my $cdsdna= $exondna;
        foreach my $ex (1 .. $#cds) {
          my($ref1,$src1,$type1,$start1,$stop1,$score1,$strand1,$phase1,$attr1)= @{$cds[$ex]};
          my($bstart,$blen)= ($start1 - 1, $stop1 - $start1 + 1);
          my $exon2dna  = substr( $refdna->seq(), $bstart, $blen);
          if($strand1 eq '-' || $strand1 < 0) {  
            $exon2dna = reverse $exon2dna;
            $exon2dna =~ tr/gatcGATC/ctagCTAG/;
            }
          $cdsdna.= $exon2dna;
          }
        $inc5 = getBestFrame( $cdsdna, $id);
        }
      
      if ($inc5 == 1) { $nt_length = 2; }
      elsif ($inc5 == 2) { $nt_length = 1; }
      else  { $nt_length = 0; }
    }
    
    my($inc5,$inc3,$elength,$frame);
    $elength = $stop - $start + 1;
		$nt_length  += $elength;
		$inc3        = $nt_length % 3;
		$inc5        = ($elength - $inc3) % 3; # only care about this one
		$frame       = ($start + $inc5) % 3;
		if ($inc5 == -1) { $inc5 = 2; }
    
    my $changed=0;
    if ($phase eq '.') {  $changed=1; }
    elsif ($phase ne $inc5 ) { 
      $changed=2; 
      warn "# phase change exon[$ix]: $phase => $inc5; $ref:$start-$stop/$strand,$type:$src,$id\n" if $verbose;
      } 
    if($changed) { $cds[$ix]->[7]= $inc5; }  
    if($ispartial && $ix == 0) { $cds[$ix]->[8] .= ";partial_gene=true"; } # 5prime_partial=true; 3prime..
    }
    
  return \@cds;
}



my @s5CodonTable = ();
BEGIN{
 @s5CodonTable = (
	 [
		 ['K','N','K','N','X',],
		 ['T','T','T','T','T',],
		 ['R','S','R','S','X',],
		 ['I','I','M','I','X',],
		 ['X','X','X','X','X',],
	],
	 [
		 ['Q','H','Q','H','X',],
		 ['P','P','P','P','P',],
		 ['R','R','R','R','R',],
		 ['L','L','L','L','L',],
		 ['X','X','X','X','X',],
	],
	 [
		 ['E','D','E','D','X',],
		 ['A','A','A','A','A',],
		 ['G','G','G','G','G',],
		 ['V','V','V','V','V',],
		 ['X','X','X','X','X',],
	],
	 [
		 ['*','Y','*','Y','X',],
		 ['S','S','S','S','S',],
		 ['*','C','W','C','X',],
		 ['L','F','L','F','X',],
		 ['X','X','X','X','X',],
	],
	 [
		 ['X','X','X','X','X',],
		 ['X','X','X','X','X',],
		 ['X','X','X','X','X',],
		 ['X','X','X','X','X',],
		 ['X','X','X','X','X',],
	],

);
}

sub ibase {
  my $c= substr($_[0],$_[1],1);
  return 0 if ($c eq 'A');
  return 1 if ($c eq 'C');
  return 2 if ($c eq 'G');
  return 3 if ($c eq 'T');
  return 4;
}  
  
sub translate {
  my($cds, $offset)= @_;
  $cds = uc($cds); ## fix chars ??
  my $aa="";
  my $aa_length = int((length($cds) - $offset) / 3);
	for (my $i = 0; $i < $aa_length; $i++) {
		my $idx = 3 * $i + $offset;
		$aa .= $s5CodonTable[ ibase($cds,$idx)][ ibase($cds,$idx+1) ][ ibase($cds,$idx+2) ];
	}
  return $aa; 
}

sub getBestFrame {
  my($cds, $id)= @_;
  my ($bestscore,$besti)= (-999,0);
  for (my $i= 0; $i<3; $i++) {
    my $pro= translate( $cds,$i );
    my $score = $pro =~ tr/*/*/; # has_internal_stops($pro);
    #is inner M bad?# $score += $pro =~ tr/M/M/;   # has_internal_starts($pro);
    $score *= -3; 
    if (substr($pro,length($pro)-1,1) eq '*') { $score += 4; } # adj internal == end
    if (substr($pro,0,1) eq 'M') { $score += 1; }
    if ($score > $bestscore) { $besti= $i; $bestscore=$score; }
 warn("# bestFrame[$i,$id]: $score ; $pro \n") if($verbose); # debug  
    }
  return $besti;
}

__END__

#item  fasta db

  # dbm based index; not platform independent
  my $dna_db = Bio::DB::Fasta->new('/path/to/fasta/files');

  my $seq     = $db->seq('CHROMOSOME_I',4_000_000 => 4_100_000);
  my $revseq  = $db->seq('CHROMOSOME_I',4_100_000 => 4_000_000);

  my $obj     = $db->get_Seq_by_id('CHROMOSOME_I');
  my $seq     = $obj->seq;
  my $subseq  = $obj->subseq(4_000_000 => 4_100_000);

  #-- or plat-independent (inherits from above) --
  my $lucene = new Bio::DB::GFF::Adaptor::lucene(xxx);
  my $dna_db = new Bio::DB::GFF::Adaptor::LuceneFasta( $fafile, _adaptor => $lucene ); 

#cut



## ATG is universal start codon for euk. nuclear genes

## for start= 0 .. 2, find ATG 
## ? and check aatrans, internal *-stops, aa[0] == 'M', aa[-1]= '*'
## parts from zoeCDS.c of Snap/I.Korf 

	if (exon->strand == '+') {
		c1 = dna->s5[exon->start];
		c2 = dna->s5[exon->start +1];
		c3 = dna->s5[exon->start +2];
		if (c1 == 0 && c2 == 3 && c3 == 2) ef.start = 1; //# ATG
	} else if (exon->strand == '-') {
		c1 = dna->s5[exon->end];
		c2 = dna->s5[exon->end -1];
		c3 = dna->s5[exon->end -2];
		if (c1 == 3 && c2 == 0 && c3 == 1) ef.start = 1; //# TAC (~ATG)
	}
#..........

	best_score = -1000000;
	best_idx   = -1;
	for (i = 0; i <= 2; i++) {
		pro[i] = zoeTranslateDNA(tx->def, tx, i);
		score = - 3 * has_internal_stops(pro[i]);
		if (pro[i]->seq[pro[i]->length -1] == '*') score++;
		if (pro[i]->seq[0] == 'M') score++;
		if (score > best_score) {
			best_score = score;
			best_idx = i;
		}
	}
	
	bt.inc5 = best_idx;
	bt.inc3 = (tx->length - best_idx) % 3;
	bt.aa   = pro[best_idx];
#..........

## find phase == inc5 for all exons:

##	/* label with correct phase and frame */
	if      (cds->inc5 == 1) nt_length = 2;
	else if (cds->inc5 == 2) nt_length = 1;
	else                     nt_length = 0;
		
	elength = 0;
  /*orig* for (i = 0; i < cds->exons->size; i++)  */
  { int i, exstart, exend, exinc;
    if(cds->strand == '-') { exstart= cds->exons->size-1; exend=-1; exinc=-1; }
    else { exstart= 0; exend= cds->exons->size; exinc=1; }
  for (i = exstart; i != exend; i += exinc )
  {
		exon        = cds->exons->elem[i];
		elength     = exon->end - exon->start + 1;
		nt_length  += elength;
		inc3        = nt_length % 3;
		inc5        = (elength - inc3) % 3;
		frame       = (exon->start + inc5) % 3;
		exon->frame = frame;
		exon->inc5  = inc5;
		if (exon->inc5 == -1) exon->inc5 = 2;
		exon->inc3  = inc3;
		if (!zoeVerifyFeature(exon)) {
			zoeWarn("exon does not validate after correcting phase & frame");
			zoeExit("%s", cds->name);
		}
   }
	}
	