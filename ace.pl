#!/usr/local/bin/perl -- w #-*-Perl-*-

#   USAGE:    
#
#   AUTHOR:  William Hayes  -  //97
#
#   INPUT:   
#            
#            
#   OUTPUT:  
#            
#
#   Subroutine:
#        
#          
#          
#            
#   Notes:   
#            
#            
#            
#   Changes:  WSH - (//97)
#	      WSH - (//97)	       
#             
#     

# Get option arguments
$genesdir="./genes";
$testsetdir="./testsets";

use Getopt::Long;
&options;

# Location of alignment program Myers and Miller, CABIOS (1989)
$align="align"; 

#$tseqset="Test";
#@sets=qw(anno genscan);

########## Global Variables ######### 
# Hold results of the accuracy evaluations
%results=();

#Sequence Object - contains all sequence related information
#%seqobj{Seq|SymSeq|Set|Unannotated}
#$seqobj{Seq}=$seq;
#@{$seqobj{SymSeq}}=@symseq;
#$seqobj{Set}{$set}{$group}{Dir|Exon|ProtSeq|ProtWarn|Strand}=$strand|{$lend}=$rend|$protseq|$protwarn;
#$seqobj{Unannotated}{$lend}=$rend;

########## Start Main ######### 
@test=&testset($tseqset);


# Check to see that each result file exists
$notexistflag=0;
print OUT "$sepline\n"; 
foreach $sequence (@test) {
    foreach $set (@sets) {
	if (!(-e "./$sequence.$set.gff") && !(-e "$genesdir/$sequence.$set.gff")) {
	    print STDERR "GFF feature file does not exist in: ./$sequence.$set.gff or $genesdir/$sequence.$set.gff\n";
	    print OUT "GFF feature file does not exist in: ./$sequence.$set.gff or $genesdir/$sequence.$set.gff\n";
	    $notexistflag=1;
	    $skipseq{$sequence}=1;
	}
    }
}

#if ($notexistflag) {
#    print STDERR "Cannot run ace without the preceding files\n";
#    exit;
#}

foreach $sequence (@test) {
    next if defined $skipseq{$sequence};
    print "Evaluating Sequence: $sequence\n";
    undef %seqobj;
    print "  Collecting\n";
    %seqobj=&collect($sequence,@sets);
    %results=&evaluate(\%seqobj,\@sets,\%results);
}

# Present results
&display($tseqset,%results);


#########################################################################
#	Subroutines
######################################################################### 

#########################################################################
#	Evaluate sequence and prediction results
######################################################################### 
sub evaluate {
    my($seqobj,$sets,$results)=@_;
    my %seqobj=%{$seqobj};
    my %results=%{$results};
    my @sets=@{$sets};

    print "  Evaluating NT\n";
    # Nucleotide level evaluation
    %results=&nteval(\%seqobj,\%results,\@sets);

    print "  Evaluating Exons\n";
    # Exon level evaluation
    %results=&exoneval(\%seqobj,\%results);

    
    # Protein level evaluation
    if (!$nogeneacc) {print "  Evaluating Prot\n";%results=&proteval(\%seqobj,\%results);}

    %results;
}

#########################################################################
#	Display results
######################################################################### 
sub display {
    my($tseqset,%results)=@_;
    
    my $sepline = "#" x 80;
    my $annoset=$main::sets[0];
    my $TN,$FN,$TP,$FP;



    for $i (1..$#main::sets) {

	$pred=$sets[$i];
	open(OUT,">>$tseqset.$pred.results") || die "Couldn't open $tseqset.$pred.results: $!\n";

	if ($notexistflag) {
	    print OUT "$sepline\n";
	}
	else {
	    print OUT "$sepline\n$sepline\n";
	}
	print OUT "Parameters:  TestSet=$tseqset  AnnoSet=$annoset  Restrict=$restrict  FirstIsoform=$firstisoform\n";
	print OUT "Date: ",`date`;
	print OUT "$sepline\n";
	# NT results
	printf OUT "Tested Sequences: $tseqset   Annotation Set: $sets[0]\n";
	printf OUT "Number of Seqs: %5d  Total length of Sequences: %10d\n",$results{NumSeq},$results{TotalSeqLen};
	print OUT "\nNucleotide results\n";

	$set=$main::sets[$i];
	$TP=$results{NT}{$set}{TP};$TN=$results{NT}{$set}{TN};
	$FP=$results{NT}{$set}{FP};$FN=$results{NT}{$set}{FN};

	$Sn= $TP+$FN==0 ? 0 : $TP/($TP+$FN);
	$Sp= $TP+$FP==0 ? 0 : $TP/($TP+$FP);
	
	# Calculate Correlation Coefficient
	$ccdenom=($TP+$FN)*($TN+$FP)*($TP+$FP)*($TN+$FN);
	$cc= $ccdenom>0 ? (($TP*$TN)-($FP*$FN))/sqrt($ccdenom) : 0;
	
	# Calculate ACP
	$acp1= $TP+$FN==0 ? 0 : $TP/($TP+$FN);
	$acp2= $TP+$FP==0 ? 0 : $TP/($TP+$FP);
	$acp3= $TN+$FP==0 ? 0 : $TN/($TN+$FP);
	$acp4= $TN+$FN==0 ? 0 : $TN/($TN+$FN);
	$acp=($acp1 + $acp2 + $acp3 + $acp4)/4;
	$ac=($acp-0.5)*2;

	printf OUT "   Set: %-10s  Sn: %5.2f  Sp: %5.2f  SnSp: %5.2f  CC: %5.2f  AC: %5.2f  TP: %5d  TN: %5d  FP: %5d  FN: %5d\n",$set,$Sn,$Sp,($Sn+$Sp)/2,$cc,$ac,$TP,$TN,$FP,$FN;
    }

    # Exon results
    print OUT "\nExon results\n";
    for $i (1..$#main::sets) {
	$set=$main::sets[$i];
	if ($results{Exon}{$annoset}{Cnt}==0) {$Sn=$Missing=0;}
	else {
	    $Sn=$results{Exon}{$set}{Exact}/$results{Exon}{$annoset}{Cnt};
	    $Missing=$results{Exon}{$set}{Missing}/$results{Exon}{$annoset}{Cnt};
	}
	if ($results{Exon}{$set}{Cnt}==0) {$Sp=$Wrong=0;}
	else {
	    $Sp=$results{Exon}{$set}{Exact}/$results{Exon}{$set}{Cnt};
	    $Wrong=$results{Exon}{$set}{Wrong}/$results{Exon}{$set}{Cnt};
	}

	printf OUT "   Set: %-10s  Sn: %5.2f  Sp: %5.2f  SnSp: %5.2f  ME: %5.2f  WE: %5.2f  Counts -> Exact: %4d  Anno: %4d  %10s: %4d  MEcnt: %4d  WEcnt: %4d\n",$set,$Sn,$Sp,($Sn+$Sp)/2,$Missing,$Wrong,$results{Exon}{$set}{Exact},$results{Exon}{$annoset}{Cnt},$set,$results{Exon}{$set}{Cnt},$results{Exon}{$set}{Missing},$results{Exon}{$set}{Wrong};
    }

    if (!$nogeneacc) {
	# Gene results
	print OUT "\nGene Results\n";
	for $i (1..$#main::sets) {
	    $set=$main::sets[$i];
	    if ($results{Prot}{$set}{LongestSeqLen}==0) {$fracsimilaritylongest=0;}
	    else {$fracsimilaritylongest=$results{Prot}{$set}{SimResidues}/$results{Prot}{$set}{LongestSeqLen};}
	    if ($results{Prot}{$set}{AnnoSeqLen}==0) {$fracsimilaritybyanno=0;}
	    else {$fracsimilaritybyanno=$results{Prot}{$set}{SimResidues}/$results{Prot}{$set}{AnnoSeqLen};}
#	    print "R: $results{Prot}{$set}{SimResidues}  A: $results{Prot}{$set}{AnnoSeqLen} L: $results{Prot}{$set}{LongestSeqLen}\n";

	    $poveralla= $results{Prot}{$annoset}{Cnt}==0 ? -1 : ($results{Prot}{$set}{Cnt}-$results{Prot}{$set}{Wrong})/$results{Prot}{$annoset}{Cnt};
	    $aoverallp= $results{Prot}{$set}{Cnt}==0 ? -1 : ($results{Prot}{$annoset}{Cnt}-$results{Prot}{$set}{Missing})/$results{Prot}{$set}{Cnt};
	    $povera= ($results{Prot}{$annoset}{Cnt}-$results{Prot}{$set}{Missing})==0 ? -1 : ($results{Prot}{$set}{Cnt}-$results{Prot}{$set}{Wrong})/($results{Prot}{$annoset}{Cnt}-$results{Prot}{$set}{Missing});
	    $MG=$results{Prot}{$annoset}{Cnt}==0 ? 0 : $results{Prot}{$set}{Missing}/$results{Prot}{$annoset}{Cnt};
	    $WG= $results{Prot}{$set}{Cnt}==0 ? 0 : $results{Prot}{$set}{Wrong}/$results{Prot}{$set}{Cnt};


	    printf OUT "   Set: %-10s  FSL: %5.3f FSA: %5.3f  P/allA: %5.3f  A/allP: %5.3f  P/A: %5.3f  MG: %5.3f  WG: %5.3f  Counts -> Anno: %4d  %10s: %4d  MGcnt: %4d  WGcnt: %4d\n",$set,$fracsimilaritylongest,$fracsimilaritybyanno,$poveralla,$aoverallp,$povera,$MG,$WG,$results{Prot}{$annoset}{Cnt},$set,$results{Prot}{$set}{Cnt},$results{Prot}{$set}{Missing},$results{Prot}{$set}{Wrong};
	}
	print OUT "\n\n";
    }
}
 

#########################################################################
#	Nucleotide level evaluation
######################################################################### 
sub nteval {
    my($seqobj,$results)=@_;
    my %seqobj=%{$seqobj};
    my %results=%{$results};
    my ($seqlen,$i,$j,$symnt,$set);
    my ($annobits,$compbits,$annocoding);

    my @symseq=@{$seqobj{SymSeq}};

    my %cnt=();
    $results{NumSeq}++;
    # Seqlen is the part that is annotated, i.e. not X'd out on the symbolic sequence
    $seqlen=0;
    foreach $i (0..$#symseq) {if ($symseq[$i] ne "X") {$seqlen++;}}
    $results{TotalSeqLen}+=$seqlen;

    foreach $i (0..$#symseq) {
	$symnt=$symseq[$i];
	next if $symnt eq "X";

	# Compare all predicted sets against Anno or base comparison set 
	for $j (1..$#main::sets) {
	    $set=$main::sets[$j];
	    @results=&bittest($symnt,$j);
	    foreach $result (@results) {
#		print "I: $i  J: $j  Set: $set  Result: $result\n";
		$cnt{$set}{$result}++;
	    }
	}
    }

    foreach $set (keys %cnt) {
	foreach $cnt (keys %{$cnt{$set}}) {
	    next if $cnt eq "TN";
	    $results{NT}{$set}{$cnt}+=$cnt{$set}{$cnt};
	}
	$results{NT}{$set}{TN}+=($seqlen*2)-$cnt{$set}{TP}-$cnt{$set}{FP}-$cnt{$set}{FN};
#	print "SeqLen: $seqlen  Symseqlen: $#symseq  TN: $results{NT}{$set}{TN} ",($seqlen*2),"-$cnt{$set}{TP}-$cnt{$set}{FP}-$cnt{$set}{FN}\n";
    }
    %results;
}
 
######################################################################### 
#	Exon level evaluation
######################################################################### 
sub exoneval {
    my($seqobj,$results)=@_;
    my %seqobj=%{$seqobj};
    my %results=%{$results};
    my ($annoset,$bitstrand);
    $annoset=$main::sets[0];
    # Loop through prediction result sets to count Exact Exons
    foreach $annogroup (keys %{$seqobj{Set}{$annoset}}) {
	$annostrand=$seqobj{Set}{$annoset}{$annogroup}{Strand};
	foreach $annolend (keys %{$seqobj{Set}{$annoset}{$annogroup}{Exon}}) {
	    $annorend=$seqobj{Set}{$annoset}{$annogroup}{Exon}{$annolend};
	    for $j (1..$#main::sets) {
		$predset=$main::sets[$j];
		foreach $predgroup (keys %{$seqobj{Set}{$predset}}) {
		    if (defined $seqobj{Set}{$predset}{$predgroup}{Exon}{$annolend} && 
			$seqobj{Set}{$predset}{$predgroup}{Exon}{$annolend} == $annorend &&
			$seqobj{Set}{$predset}{$predgroup}{Strand} eq $annostrand) {
			$results{Exon}{$predset}{Exact}++;
		    }
		}
	    }
	}
    }

    # Missing or Wrong Exons
    @symseq=@{$seqobj{SymSeq}};
    for $i (0..$#main::sets) {
	$set=$main::sets[$i];
	foreach $group (keys %{$seqobj{Set}{$set}}) {
	    foreach $lend (keys %{$seqobj{Set}{$set}{$group}{Exon}}) {
		$rend=$seqobj{Set}{$set}{$group}{Exon}{$lend};
		$results{Exon}{$set}{Cnt}++;
		undef %missedflag;
		$overlapflag=0;
		for $j ($lend-1..$rend-1) {
		    $symnt=$symseq[$j];
		    next if $symnt eq "X";
		    # Is annotation base considered coding?
		    $abitcoding = $symnt & 2**0 ? 1 : 0;
		    $abitstrand = $symnt & 2**1 ? 1 : 0;
		    if ($i>0) {
			# Determine if base is predicted coding
			$pbitcoding = $symnt & 2**($k*2) ? 1 : 0;
			$pbitstrand = $symnt & 2**($k*2+1) ? 1 : 0;
		    }
		    if ($i==0) {
			# Missing exons test
			for $k (1..$#main::sets) {
			    $missedset=$main::sets[$k];
			    if (!defined $missedflag{$missedset}) {$missedflag{$missedset}=0;}
			    # Set predicted set strand=1 if plus strand coding
			    $pbitcoding = $symnt & 2**($k*2) ? 1 : 0;
			    $pbitstrand = $symnt & 2**($k*2+1) ? 1 : 0;
			    if ($pbitcoding && $abitstrand==$pbitstrand) {
				$missedflag{$missedset}=1;
			    }
			}
		    }
		    # Wrong exons test
		    elsif ($abitcoding && $abitstrand==$pbitstrand) {$overlapflag++;}
		}
		# Wrong exons count
		if ($verbose==3) {
		    print "I: $i  Group: $group Lend: $lend  Overlap: $overlapflag Symnt: $symnt ";
		    print "AC:$abitcoding AS:$abitstrand PS:$pbitstrand\n";
		}
		
		if ($set ne $annoset && $overlapflag==0) {
		    $results{Exon}{$set}{Wrong}++;
		    if ($verbose==3) {
			print "I: $i  $set ne $annoset Group: $group  LEnd: $lend Wrong: $lend $rend\n\n";
		    }
		}
		# Missing exons count
		elsif ($i == 0) {
		    foreach $set (keys %missedflag) {
			if ($missedflag{$set}==0) {
			    $results{Exon}{$set}{Missing}++;
			    if ($verbose==3) {
				print "I: $i  Group: $group  LEnd: $lend Missing: $lend $rend\n\n";
			    }
			}
		    }
		}
	    }
	}
    }
	
    %results;
}

#########################################################################
#	Protein level evaluation
#
#       Preliminary Proposal:
#
#	1.  List %similarity from the best (or all overlapping)
#	    predictions - using align program 
#           (based upon the longest aligned seq length)
#	2.  List predicted genes overlapping annotated genes vs all annotated for P/A
#	3.  List Missing genes and Wrong genes
#           only genes not overlapping with either annotated (for Wrong)
#           or predicted (for Missing) genes 
######################################################################### 
sub proteval {
    my($seqobj,$results)=@_;
    my %seqobj=%{$seqobj};
    my %results=%{$results};

    my %overlap=();

    $annoset=$main::sets[0];

    # Determine which genes overlap with which gene predictions
    for $i (1..$#main::sets) {
	$set=$main::sets[$i];
	foreach $group (keys %{$seqobj{Set}{$set}}) {
	    $strand=$seqobj{Set}{$set}{$group}{Strand};
	    foreach $lend (keys %{$seqobj{Set}{$set}{$group}{Exon}}) {
		$rend=$seqobj{Set}{$set}{$group}{Exon}{$lend};
		foreach $annogroup (keys %{$seqobj{Set}{$annoset}}) {
		    $annostrand=$seqobj{Set}{$annoset}{$annogroup}{Strand};
		    foreach $annolend (keys %{$seqobj{Set}{$annoset}{$annogroup}{Exon}}) {
			$annorend=$seqobj{Set}{$annoset}{$annogroup}{Exon}{$annolend};
			$overlapresult=&overlap($lend,$rend,$annolend,$annorend);
			if ($strand eq $annostrand && $overlapresult=~/PL|E|S|B|PG/) {
			    $overlap{$set}{$group}{$annoset}{$annogroup}=1;
			    $overlap{$annoset}{$annogroup}{$set}{$group}=1;
			}
		    }
		}
	    }
	}
    }

    # Count Missing genes
    foreach $annogroup (keys %{$seqobj{Set}{$annoset}}) {
	$results{Prot}{$annoset}{Cnt}++;
	foreach $compset (keys %{$seqobj{Set}}) {
	    next if $compset eq $annoset;
	    if (!defined $overlap{$annoset}{$annogroup}{$compset}) {
		$results{Prot}{$compset}{Missing}++;
	    }
	}
    }

    # Count Wrong genes
    for $i (1..$#main::sets) {
	$set=$main::sets[$i];
	foreach $group (keys %{$seqobj{Set}{$set}}) {
	    $results{Prot}{$set}{Cnt}++;
	    if (!defined $overlap{$set}{$group}{$annoset}) {
		$results{Prot}{$set}{Wrong}++;		
	    }
	}
    }

    # get %similarity scores
    foreach $group (keys %{$seqobj{Set}{$annoset}}) {
	foreach $overlapset (keys %{$overlap{$annoset}{$group}}) {
	    $maxsimres=-1;$cnt=0;
	    foreach $overlapgroup (keys %{$overlap{$annoset}{$group}{$overlapset}}) {
		$cnt++;
		next if $seqobj{Set}{$annoset}{$group}{ProtSeq}!~/\S/;
		next if $seqobj{Set}{$overlapset}{$overlapgroup}{ProtSeq}!~/\S/;
		($simres,$annoseqlen,$longestseqlen)=&align($seqobj{Set}{$annoset}{$group}{ProtSeq},
					     $seqobj{Set}{$overlapset}{$overlapgroup}{ProtSeq});
		if ($simres>$maxsimres) {
		    $maxsimres=$simres;
		    $maxlongestseqlen=$longestseqlen;
		    $maxannoseqlen=$annoseqlen;
		}
	    }
	    if ($maxsimres>-1) {
		$results{Prot}{$overlapset}{SimResidues}+=$maxsimres;
		$results{Prot}{$overlapset}{AnnoSeqLen}+=$maxannoseqlen;
		$results{Prot}{$overlapset}{LongestSeqLen}+=$maxlongestseqlen;
	    }
	    else {
		print "Couldn't find overlapping prediction for $group in $overlapset\n";
	    }
	}
    }

    %results;
}

#########################################################################
#	Get number of aligned residues using the 'align' program
######################################################################### 
sub align {
    my($stdseq,$compseq)=@_;  # Annotated sequence, Predicted sequence order
    my($seqlen,$simres,$identity);
    open(STD,">/tmp/$$.stdseq") || die "Couldn't open /tmp/$$.stdseq: $!\n"; 
    print STD ">Standard sequence for align\n$stdseq\n";
    close STD;

    open(COMP,">/tmp/$$.compseq") || die "Couldn't open /tmp/$$.compseq: $!\n"; 
    print COMP ">Compare sequence for align\n$compseq\n";
    close COMP;

    system("$main::align -O /tmp/$$.align.out  /tmp/$$.stdseq /tmp/$$.compseq  >/dev/null 2>/dev/null");
    open(IN,"/tmp/$$.align.out") || die "Couldn't open /tmp/$$.align.out: $!\n"; 
    while (<IN>) {
	if (/^\s*(\S+)\% identity;/) {$identity=$1;last;}
    }
    close IN;

    system("/bin/rm /tmp/$$.*");

    $stdseq=~s/\*$//;
    chomp $stdseq;
    chomp $compseq;
    $stdseqlen=length $stdseq;
    $longestseqlen= length $stdseq > length $compseq ? length $stdseq : length $compseq;
    $simres=($identity/100)*$longestseqlen;
    $simres=int($simres+0.5);
    ($simres,$stdseqlen,$longestseqlen);
}

#########################################################################
#	Collect sequence object information
######################################################################### 
sub collect {
    my($sequence,@sets)=@_;
    my ($seq);
    # Read DNA sequence
    open(SEQ,"./$sequence.seq") || open(SEQ,"$genesdir/$sequence.seq") || die "Couldn't open ./$sequence.seq or $genesdir/$sequence.seq: $!\n";
    while (<SEQ>) {
	next if /^>/;
	next if /^\s*\#/;
	$seq.=$_;
    }
    $seq=~s/\d//g;$seq=~s/\s//g;
    $seqobj{Seq}=$seq;

    # Read in each GFF file
    foreach $sets (@sets) {
	($flag,%seqobj)=&collectGFF($sequence,\%seqobj,$sets);
	$altsplflag=$altsplflag | $flag;
    }
    # Select best Isoform between each Annotated Isoform set and set of Isoform predictions
    if ($altsplflag && $main::firstisoform) {
	%seqobj=&selectFirstIsoform(\%seqobj,\@sets);
    }
    elsif ($altsplflag) {
	%seqobj=&selectIsoform(\%seqobj,\@sets);
    }

    # Collect protein sequence related to each CDS if gene level accuracy is to be determined
    if (!$nogeneacc) {%seqobj=&cds2prot(%seqobj);}

    # Mask sequence against annotation/predictions and unannotated sequence
    %seqobj=&maskseq(\%seqobj,\@sets);

    %seqobj;
}

#########################################################################
#	Select First occuring Isoform
######################################################################### 
sub selectFirstIsoform {
    my($seqobj,$sets)=@_;
    my %seqobj=%{$seqobj};
    my @sets=@{$sets};
    my %best;
    my ($prefix,$isonumber,$set,$group);

    foreach $set (@sets) {
	%best=();
	foreach $group (keys %{$seqobj{Set}{$set}}) {
	    if ($group=~/^(\S+)\.AltSpl(\d+)$/) {
		$prefix=$1;$isonumber=$2;
		if ($isonumber < $best{$prefix}{Number} || !defined $best{$prefix}{Number}) {
		    $best{$prefix}{Group}=$group;
		    $best{$prefix}{Number}=$isonumber;
		}
	    }
	}
	foreach $group (keys %{$seqobj{Set}{$set}}) {
	    if ($group=~/^(\S+)\.AltSpl(\d+)$/) {
		$prefix=$1;$isonumber=$2;
		if (!defined $best{$prefix}{Group}) {
		    print "Something wrong with $group in &selectFirstIsoform subroutine\n";
		}
		elsif ($group ne $best{$prefix}{Group}) {
		    delete $seqobj{Set}{$set}{$group};
		}
	    }
	}	
    }
    %seqobj;
}

#########################################################################
#	Select best Isoform between annotated and predicted Isoform sets
######################################################################### 
sub selectIsoform {
    my($seqobj,$sets)=@_;
    my %seqobj=%{$seqobj};
    my @sets=@{$sets};

    # Select via number of matching Exons
    $anno=$sets[0]; $pred=$sets[1];
    foreach $set ($anno,$pred) {
	%max=();
	if ($verbose==3) {print "Set: $set\n";}
	foreach $group (sort keys %{$seqobj{Set}{$set}}) {
	    if ($verbose==3) {print "Group: $group\n";}
	    $cnt=0;
	    next if $group!~/AltSpl/i;
	    $prefix=$group;$prefix=~s/\.AltSpl\d+//;
	    $strand=$seqobj{Set}{$set}{$group}{Strand};
	    foreach $lend (sort {$a<=>$b} keys %{$seqobj{Set}{$set}{$group}{Exon}}) {
		$rend=$seqobj{Set}{$set}{$group}{Exon}{$lend};
		if ($verbose==3) {print "Lend: $lend Rend: $rend  Strand: $strand\n";}
		# Test each of the predicted CDS's including different Isoforms
		foreach $compset ($anno,$pred) {
		    next if $compset eq $set;
		    if ($verbose==3) {print "Compset: $compset\n";}
		    foreach $compgroup (sort keys %{$seqobj{Set}{$compset}}) {
			$compstrand=$seqobj{Set}{$compset}{$compgroup}{Strand};
			if ($verbose==3) {print "Compgroup: $compgroup\n";}
			foreach $complend (sort {$a<=>$b} keys %{$seqobj{Set}{$compset}{$compgroup}{Exon}}) {
			    $comprend=$seqobj{Set}{$compset}{$compgroup}{Exon}{$complend};
			    if ($verbose==3) {print "CompLend: $complend  CompRend: $comprend  CompStrand: $compstrand\n";}
			    if ($lend==$complend && $rend==$comprend && $strand eq $compstrand) {
				$cnt++;
			    }
			    if ($verbose==3) {print "Cnt: $cnt\n";}
			}
			# Collect best scoring Isoforms
			if ($max{$prefix}{Cnt}<$cnt) {
			    $max{$prefix}{Group}=$group;
			    $max{$prefix}{Cnt}=$cnt;
			}
		    }
		}
	    }
	}
	foreach $prefix (keys %max) {
	    $bestgroup=$max{$prefix}{Group};
	    foreach $group (keys %{$seqobj{Set}{$set}}) {
		next if $bestgroup eq $group;
		if ($group=~/$prefix\.AltSpl\d+/i) {
		    if ($verbose==3) {print "Test: $sequence Set: $set Group: $group Prefix: $prefix Best: $max{$prefix}{Group}\n";}
		    delete $seqobj{Set}{$set}{$group};
		}
	    }
	}
    }
	
    %seqobj;
}

#########################################################################
#	Collect GFF information for each comparison set
######################################################################### 
sub collectGFF {
    my($sequence,$seqobj,$set)=@_;
    my %seqobj=%{$seqobj};
    my %annotated=();
    my %phase;
    my $altsplflag=0;

    open(GFF,"./$sequence.$set.gff") || open(GFF,"$genesdir/$sequence.$set.gff") || die "Couldn't open ./$sequence.$set.gff or $genesdir/$sequence.$set.gff: $!\n";
    while (<GFF>) {
	next if /^\s*\#/;
	next if /^\s*$/;
	chomp $_;
	@gff=split(/\t/,$_);
	next if $gff[8]=~/^S\.\d+/i;

	# Isoform checks - cannot handle more than 2 test sets when trying to
	#   find best match between annotation and prediction Isoforms 
	#   unless taking the first Isoform out of any set of Isoforms
	if (!$altsplflag && $gff[8]=~/AltSpl/i && $#sets>=2 && !$firstisoform) {
	    print STDERR "Cannot have more than two <CompareSets> when Isoforms exist for the annotation or prediction CDS's\n";
	    &usage;
	}
	# Set flag is Isoforms exist
	elsif (!$altsplflag && $gff[8]=~/AltSpl/i) {$altsplflag=1;}

	# Only collect GFF features that have CDS or annotated in the feature field
	if ($gff[2]=~/CDS/) {
	    $seqobj{Set}{$set}{$gff[8]}{Strand}=$gff[6];
	    $seqobj{Set}{$set}{$gff[8]}{Exon}{$gff[3]}=$gff[4];
	    $seqobj{Set}{$set}{$gff[8]}{Phase}{$gff[3]}=$gff[7];
	}
	elsif ($set eq $main::restrict && $gff[2]=~/annotated/i) {$annotated{$gff[3]}=$gff[4];}
    }
    close GFF;
    
    # Define Unannotated sequence from the GFF 'annotated' features
    if (defined %annotated) {
	$start=1;
	foreach $lend (sort {$a<=>$b} keys %annotated) {
	    $rend=$annotated{$lend};
	    if ($lend-$start>0) {$seqobj{Unannotated}{$start}=$lend-1;}
	    $start=$rend+1;
	}
	$seqlen=length $seqobj{Seq};
	if ($seqlen-$start>0) {$seqobj{Unannotated}{$start}=$seqlen;}
    }
    undef %annotated;

    ($altsplflag,%seqobj);
}

#########################################################################
#	Mask sequence with CDS locations
######################################################################### 
sub maskseq {
    my($seqobj,$sets)=@_;
    my %seqobj=%{$seqobj};
    my @sets=@{$sets};
    my $seqlen=length $seqobj{Seq};
    my @symseq=();
    my ($group,$lend,$rend,$i,$j,$mask);

    # Initialize @symseq
    foreach $i (0..$seqlen-1) {$symseq[$i]=0;}

    # Mask each Compare group against the SYMbolic SEQuence (symseq)
    for ($i=0;$i<=2*$#sets;$i+=2) {
	$set=$sets[$i/2];
	foreach $group (keys %{$seqobj{Set}{$set}}) {
	    foreach $lend (keys %{$seqobj{Set}{$set}{$group}{Exon}}) {
		$rend=$seqobj{Set}{$set}{$group}{Exon}{$lend};
		foreach $j ($lend-1..$rend-1) {
		    #Use two bit mask for determining Coding/Noncoding and +/- strand
		    if ($seqobj{Set}{$set}{$group}{Strand} eq "+") {$mask=2**$i + 2**($i+1);}
		    else {$mask=2**$i;}
		    $symseq[$j]=$symseq[$j] | $mask;
		}
	    }
	}
    }

    # Mask unannotated sequence with X's
    foreach $lend (keys %{$seqobj{Unannotated}}) {
	$rend=$seqobj{Unannotated}{$lend};
	for $i ($lend-1..$rend-1) {
	    $symseq[$i]="X";
	}
	# Remove Exons that are fully overlapped by unannotated sequence
	foreach $set (keys %{$seqobj{Set}}) {
	    foreach $group (keys %{$seqobj{Set}{$set}}) {
		foreach $exonlend (sort {$a<=>$b} keys %{$seqobj{Set}{$set}{$group}{Exon}}) {
		    $exonrend=$seqobj{Set}{$set}{$group}{Exon}{$exonlend};
		    $overlap=&overlap($lend,$rend,$exonlend,$exonrend);
		    if ($overlap=~/[ES]/) {
			delete $seqobj{Set}{$set}{$group}{Exon}{$exonlend};
		    }
		}
	    }
	}
    }

    # Remove Group entries that do not have any Exon members 
    foreach $set (keys %{$seqobj{Set}}) {
	foreach $group (keys %{$seqobj{Set}{$set}}) {
	    $exoncnt=0;
	    foreach $exonlend (sort {$a<=>$b} keys %{$seqobj{Set}{$set}{$group}{Exon}}) {
		$exoncnt++;
	    }
	    if ($exoncnt==0) {delete $seqobj{Set}{$set}{$group};}
	}
    }
    
    @{$seqobj{SymSeq}}=@symseq;
    
    %seqobj;
}

#########################################################################
#	Read in list of test sequences
######################################################################### 
sub testset {
    my($testset)=@_;
    open(IN,"./$testset.set") || open(IN,"$testsetdir/$testset.set") || die "Couldn't open ./$testset.set or $testsetdir/$testset.set: $!\n"; 
    while (<IN>) {
	next if /^\s*\#/;
	next if /^\s*$/;
	if (/^\s*(\S+)\s*/) {push @set,$1;}
    }
    @set;
}

#########################################################################
#	Collect protein sequence for each CDS
######################################################################### 
sub cds2prot {
    my(%seqobj)=@_;

    my($seq,$dir,$lend,$rend,$cds,$warning);

    $seq=$seqobj{Seq};
    foreach $set (keys %{$seqobj{Set}}) {
	foreach $group (keys %{$seqobj{Set}{$set}}) {
	    $cds="";$phase=-1;
	    $strand=$seqobj{Set}{$set}{$group}{Strand};
	    foreach $lend (sort {$a<=>$b} keys %{$seqobj{Set}{$set}{$group}{Exon}}) {
		$rend=$seqobj{Set}{$set}{$group}{Exon}{$lend};
		if (($phase==-1 && $strand eq "+") || $strand eq "-") {
		    $phase=$seqobj{Set}{$set}{$group}{Phase}{$lend};
		} #SRIDHAR fix
		$cds.=substr($seq,$lend-1,$rend-$lend+1);
	    }

	    # If complement strand take reverse complement of $cds
	    if ($strand eq "-") {$cds=&revcomp($cds);}
	    
#	    print "Set: $set   Group: $group  Phase: $phase  Remainder: $remainder\n";
	    if ($set!~/genie/i) {substr($cds,0,$phase)="";}
	    $remainder=(length $cds)%3;
	    substr($cds,-$remainder,$remainder)="";

	    # Check CDS for ATG start, mult of three and a stop codon
	    $problem=&checkcds($cds);

	    if($set eq $annoset && $problem =~ /Internal-Stop-Codon/) { #SRIDHAR addition
		print OUT "Seq: $sequence Set: $set  CDS=$group has Internal Stop Codon\n";
	    }
	    # If CDS is not Incomplete, determine protein translation of CDS
	    elsif ($problem ne "") {
		print STDERR "Seq: $sequence Set: $set  Group: $group  CDS translation to protein Problem: $problem\n";
	    }

	    $prot=&protran($cds);
	    $seqobj{Set}{$set}{$group}{ProtSeq}=$prot;

	}
    }
    %seqobj;
}

#########################################################################
#	Check Coding/Noncoding, +/- strand bits for symbolic sequence
#	  each set including the annotation set has bits set in the 
#	  @symseq array corresponding to each base of the sequence
######################################################################### 
sub bittest {
    my($symnt,$set)=@_;
    my($annobits,$compbits,$annocoding,@results);

    # Check if Annotation is coding or not
    $annocoding = $symnt & 2**($k*2) ? 1 : 0;
    $annostrand = $symnt & 2**($k*2+1) ? 1 : 0;
    # Determine pred sequence bits
    $predcoding = $symnt & 2**($set*2) ? 1 : 0;
    $predstrand = $symnt & 2**($set*2+1) ? 1 : 0;


    # Determine result
    if ($annocoding==1 &&  $predcoding==1 && $annostrand==$predstrand)            {push @results,"TP";}
    if ($annocoding==0 && ($annostrand!=$predstrand || $annocoding==$predcoding)) {push @results,"TN";}
    if ($predcoding==1 && ($annostrand!=$predstrand || $annocoding==0))           {push @results,"FP";}
    if ($annocoding==1 && ($annostrand!=$predstrand || $predcoding==0))           {push @results,"FN";}
    
    @results;
}

#########################################################################
#	Check CDS for ATG start, mult of three and a stop codon
######################################################################### 
sub checkcds {
    my($cds)=@_;
    my $problem="";
    my($i, $codon); #SRIDHAR
    my $flag=0;
    my $stop=substr($cds,-3,3);
    $codon = substr($cds,0,3);
    if ($cds!~/^ATG/) {$problem.="Non-ATG start $codon ";}
    if ($stop!~/TGA|TAG|TAA/) {$problem.="No stop codon ";}
    if ((length $cds)%3!=0) {$problem.="Non-triplet ";}

    for($i=0;$i<length($cds)-3;$i+=3) {
	$codon = substr($cds,$i,3);
	if($codon =~ /TGA|TAG|TAA/) {
	    $problem.="Internal-Stop-Codon $codon at position $i ";
	    last;
	}
    }
    $problem;
}


#########################################################################
# Overlap subroutine and key for result - checks two subsequences for overlaps
#        -----         CDS Annotation - 2nd subseq is L|PL|E|S|B|PG|G than 1st subseq
#  -----               L = less than
#      -----           PL = partially less than
#        -----         E = equal to
#         ---          S = smaller than  (fully contained)
#      ---------       B = bigger than (fully overlaps)
#           -----      PG = partially greater than
#               -----  G = greater than
########################################################################## 
sub overlap {
    my($lend1,$rend1,$lend2,$rend2)=@_;
    my @ends=($lend2,$rend2,$lend1,$rend1);
    my $result="";
    if    ($ends[1]<$ends[2])                    {$result="L";}
    elsif ($ends[0]<$ends[2]  && $ends[1]<$ends[3])  {$result="PL";}
    elsif ($ends[0]==$ends[2] && $ends[1]==$ends[3]) {$result="E";}
    elsif ($ends[0]>=$ends[2]  && $ends[1]<=$ends[3])  {$result="S";}
    elsif ($ends[0]<=$ends[2]  && $ends[1]>=$ends[3])  {$result="B";}
    elsif ($ends[0]>$ends[3])                    {$result="G";}
    elsif ($ends[0]>$ends[2]  && $ends[1]>$ends[3])  {$result="PG";}
    $result;
}


##########################################################################################
#  Translate nt's from region or orf into protein sequence
#       what cannot be translated is turned into X's
#       stop codons are *'s
##########################################################################################
sub protran {
    my($seq)=@_;
    my($prot,$i,%tranaa);
    $seq=~ tr/a-z/A-Z/;
    $seq=~tr/a-z/A-Z/;

    %tranaa=qw(TTT F TTC F TTA L TTG L CTT L CTC L CTA L CTG L ATT I ATC I ATA I 
	       ATG M GTT V GTC V GTA V GTG V TCT S TCC S TCA S TCG S CCT P CCC P 
	       CCA P CCG P ACT T ACC T ACA T ACG T GCT A GCC A GCA A GCG A TAT Y 
	       TAC Y TAA * TAG * CAT H CAC H CAA Q CAG Q AAT N AAC N AAA K AAG K 
	       GAT D GAC D GAA E GAG E TGT C TGC C TGA * TGG W CGT R CGC R CGA R 
	       CGG R AGT S AGC S AGA R AGG R GGT G GGC G GGA G GGG G);

    $prot="";
    for ($i=0;$i<length($seq);$i+=3) {
	$triplet=substr($seq,$i,3);
	if ($tranaa{$triplet} eq "") {$prot.="X";}
	else {$prot.=$tranaa{$triplet};}
    }
    $prot;
}


##########################################################################################
#  Get reversed complement of DNA strand
##########################################################################################
sub revcomp {
# $revcompseq = &revcomp($seq); # Get reverse complement of sequence 
    my($cseq) = @_;
    $cseq=~ tr/a-z/A-Z/;
    $cseq = reverse($cseq);
    $cseq =~ tr/A/1/;$cseq =~ tr/T/A/;$cseq =~ tr/1/T/;		
    $cseq =~ tr/G/2/;$cseq =~ tr/C/G/;$cseq =~ tr/2/C/;	
    $cseq;
}

#########################################################################
#	Print Sequence in FASTA format with 60bp lines SRIDHAR
######################################################################### 

sub printSeq{
    my($seq, $name) = @_;
    my($i);
    my($len) = length($seq);
    print ">$name\n";
    for($i=0;$i < $len;$i+=60){
	print substr($seq,$i,60)."\n";
    }
}
#########################################################################
#	Get options and print usage if necessary
#########################################################################
sub options {
    my $help = 0;		# handled locally
    my $ident = 0;		# handled locally

    # Process options.
    if ( @ARGV > 0 && $ARGV[0] =~ /^[-+]/ ) {
	&usage 
	  unless &GetOptions (	'help'  => \$help,
				'tseqset=s' => \$tseqset,
				'annoset=s' => \$annoset,
				'genesdir=s' => \$genesdir,
				'testsetdir=s' => \$testsetdir,
				'printsets' => \$printseqsets,
				'firstisoform' => \$firstisoform,
				'restrict=s' => \$restrict,
				'nogeneacc'=> \$nogeneacc,
				'verbose' => \$verbose,
			     )
		&& !$help;
    }
    if (defined $printseqsets) {
	opendir(DIR,"$testsetdir");
	@files = grep(/\.set/,readdir(DIR));
	closedir(DIR);
	opendir(DIR,".");
	@localfiles = grep(/\.set/,readdir(DIR));
	closedir(DIR);
	print "Test Sets: ";
	foreach $fn (@files,@localfiles) {
	    $fn=~s/\.set//;
	    print "$fn ";
	}
	print "\n";
	exit 1;
    }
    if (!defined $tseqset) {print STDERR "Please give a Test sequence set name, i.e. $0 -set Tigger\n";&usage;}
    if (!defined $annoset) {print STDERR "Please give an annotation set name or base comparison result set\n";
    }
    if ($#ARGV<0) {
	print STDERR "Please give at least one <PredictionSet>\n";
	&usage;
    }
    else {
	@pset=@ARGV;
    }
    if (!defined $verbose) {$verbose=0;}
    if (!defined $restrict) {$restrict=$annoset;}

    # First value of @sets is the base comparison set or Annotated set against which the other
    #   sets will be compared.  One can only have one prediction set if Isoforms are involved.
    @sets=($annoset,@pset);
    
}

sub usage {
    print STDERR <<EndOfUsage;

Usage: $0 [ options ] <PredictionSets>
    -help		                    This message
    -tseqset <SeqSET>                       Test set of sequences to evaluate                  
    -annoset <AnnotatedSet>                 Set of Annotations to compare predictions against
    -genesdir <GeneResultsDir>              Directory where all Annotations, Prediction results and Sequences are stored
    -testsetdir <TestSetDir>                Directory where Test Sets are stored
    -restrict <AnnotatedSet|PredictionSet>  See below for description (Default is <AnnotatedSet>)  
    -firstisoform                           Select first Isoform for accuracy evaluation
    -nogeneacc                              Don't compute gene (protein) level accuracy
    -printsets                              Used alone to print out what test sets are available
    -verbose                                Prints messages to STDOUT to keep one abreast of 
                                            what\'s going on


This program compares GFF files, annotation vs prediction, to determine the accuracy of the gene prediction method.  A test set of sequences, <SeqSET>, are used where the names of the sequences are <NAME> with the associated Sequence and <AnnotatedSet|PredictionSet> files <NAME>.seq, <NAME>.anno.gff, <NAME>.genscan.gff.  The <NAME>.seq file is a FASTA file.

The <CompareSets> GFF files need to have the standard GFF format.  Blank lines and lines starting with a \# are ignored.  All fields are separated via the TAB character.  The feature field must have CDS in it to become part of a predicted gene.  The group field for a set of CDS exons must be the same for all CDS exons from the same gene.  Isoforms (Alternative Splice forms) must have a ".AltSpl0"..".AltSpl<N>" post-fixed to the group identifier for a CDS.  One can also indicate "Annotated" sequence using Annotated in the Feature field with the appropriate LEnd and REnd.

The GFF Annotated feature allows one to predict genes in a genomic sequence where the annotation is only known for a fraction of the sequence.  All predictions, annotations will be ignored for any sequence that is not in the GFF Annotated range of sequence.  GFF Annotated features will be taken only from one set, <AnnotatedSet> or a <PredictionSet>.  If one has generated a new set of prediction GFF's called "perfect", but these predictions only apply to a small portion of the sequence, then one can put GFF Annotation features in the GFF file of the <SeqSet>.perfect.gff file indicating the range of sequence to which the predictions apply and also set '-restrict perfect'.  This way one does not get hit for missed annotations that fall outside of the sequence that was predicted for genes.

If Isoforms exist for the annotation or prediction CDS's, one can only give one <AnnotationCompareSet> and one <PredictionCompareSet>.

If the 'align' program is not in your path, please explicitly set the location of this program in ace.pl by setting the $align variable at the beginning of ace.pl

Results of the comparisons will be given in the form of: 

NT Results

Exon Results

Gene Results

FSA: 0.000  (Similar residues/annotated gene length) using 'align' program (Myers & Miller, 1988) 
FSL: 0.000  (Similar residues/longest {predicted or annotated} gene length) using 'align' program (Myers & Miller, 1988) 

P/allA: 1.000  Number of predicted genes that overlap annotated genes divided by total number of annotated genes 

A/allP: 0.250  Number of annotated genes that overlap predicted genes divided by total number of predicted genes 

P/A: Number of predicted genes that overlap annotated genes divided by number of annotated genes that overlap predicted genes

MG: 0.000  Number of annotated genes with no predicted gene overlap divided by the total number of annotated genes

WG: 0.750  Number of wrong gene predictions with no annotated gene overlap divided by the total number of gene predictions


EndOfUsage

    if (!$ENV{'POSIXLY_CORRECT'}) {
          print STDERR <<EndOfAbbrev
       BTW, you can abbreviate option names!

EndOfAbbrev
    }
    exit 1;
}



## POD Documentation:
=head1 NAME

ace - Perl script for generating Accuracy Evaluations of GFF\'d Gene Predictions

=head1 SYNOPSIS

=head2 Parsing Blast reports


=head1 INSTALLATION
  
=head1 DESCRIPTION

B<FEATURES:>


B<Results:>

NT Results

Exon Results

Gene Results

FSA: 0.000  (Similar residues/annotated gene length) using 'align' program (Myers & Miller, 1988) 
FSL: 0.000  (Similar residues/longest {predicted or annotated} gene length) using 'align' program (Myers & Miller, 1988) 

P/allA: 1.000  Number of predicted genes that overlap annotated genes divided by total number of annotated genes 

A/allP: 0.250  Number of annotated genes that overlap predicted genes divided by total number of predicted genes 

P/A: Number of predicted genes that overlap annotated genes divided by number of annotated genes that overlap predicted genes

MG: 0.000  Number of annotated genes with no predicted gene overlap divided by the total number of annotated genes

WG: 0.750  Number of wrong gene predictions with no annotated gene overlap divided by the total number of gene predictions

=head1 USAGE

=head1 AUTHOR

William S. Hayes
Bioinformatics, SmithKline-Beecham R&D
William_S_Hayes@sbphrd.com

= head1 Notes


=cut
## End of POD Documentation
