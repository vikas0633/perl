#!/usr/local/bin/perl

#use strict;
use warnings;
use Getopt::Long;
use constant USAGE =><<END;

#http://www.daimi.au.dk/~kmt/scriptlist/gffintersect.txt
SYNOPSIS:

 gffintersect.pl [OPTIONS] file1 file2

DESCRIPTION:

Finds intersection of two GFF files. Prints lines from file2 that
intersect (or do not intersect) with file1

OPTIONS:

       --not
             To print lines from file2 that don't match the criteria rather than lines that do
       --all
            To print everything from file2
       --self
            for comparisons vs self (omit file2); also adds \"self&\" prefix to name1
       --minfrac1 xxx
            To specify minimum fractional overlap for file1 features
       --minfrac2 xxx
            To specify minimum fractional overlap for file2 features
       --near xxx
             To extend definition of \"overlap\" to nearby GFFs
       --minhits xxx
            To specify minimum number of hits required to
            print lines from file2 (default is 1)
       --maxhits xxx
            To specify maximum number of hits required to
            print lines from file2 (default is 1)
       --name1
            If file1 is a name that won't look good in \"intersect(..)=(....)\" messages
       --quiet
            To suppress \"intersect(..)=(....)\" messages altogether
       --print1
            To copy all lines from file1 to output
       --regex
            To specify a (perl) regular expression determining how
            much of the reference id in field one that is to be taken
            as reference id and how much that is just an exon, intron
            or alt.spl (Eg. --regex '^[A-Za-z0-9]+_[A-Za-z0-9]+''
            specifying that only NC_003279 of NC_003279_3 is to
            determine the reference sequence.)
       --shutit
            To suppress GFF format warnings
       --presorted
            The GFF files are already sorted.

EXAMPLES:



AUTHOR:

Ian Holmes / (Kasper Munch)

COPYRIGHT:

This program is free software. You may copy and redistribute it under
the same terms as Perl itself.

END

my $not = 0;
my $all = 0;
my $presorted = 0;
my $self = 0;
my $quiet = 0;
my $byparent = 0;
my $print1 = '';
my $regex = '';
my $minfrac1 = 0;
my $minfrac2 = 0;
my $minhits = 1;
my $maxhits = 1;
my $near = 0;
my $name1 = '';
my $shutit = 0;

GetOptions("not" => \$not,
           "all" => \$all,
           "presorted" => \$presorted,
           "self" => \$self,
           "quiet" => \$quiet,
           "byparent" => \$byparent,
           "print1" => \$print1,
           "regex=s" => \$regex,
           "minfrac1=f" => \$minfrac1,
           "minfrac2=f" => \$minfrac2,
           "minhits=i" => \$minhits,
           "maxhits=i" => \$maxhits,
           "near=i" => \$near,
           "name1=s" => \$name1,
           "shutit" => \$shutit) or die USAGE;

if (@ARGV==0 && $self) {
  die "Sorry - perl doesn't buffer STDIN, so you'll have to use a temporary file for that. Died";
}
if (@ARGV==1) {
  push @ARGV, $self ? $ARGV[0] : '-';
}
@ARGV==2 or die USAGE;
($file1,$file2) = @ARGV;

$name1 = $file1 unless $name1;
if ($self) {
  $name1 = "self&$name1";
}

# Put the files into arrays of the lines:
undef $/;
open my $fh1, "$file1";
$_ = <$fh1>;
$file1 = [split /^/];
open my $fh2, "$file2";
$_ = <$fh2>;
$file2 = [split /^/];
$/ = "\n";

# Make mappings from line to linenr:
for (my $i=0; $i < @$file1; $i++) {
  $file1{$$file1[$i]} = $i+1;
}
for (my $i=0; $i < @$file2; $i++) {
  $file2{$$file2[$i]} = $i+1;
}

unless ($presorted) {
  # Sort the arrays:
  $file1 = [sort generic @$file1];
  $file2 = [sort generic @$file2];
}

# $line1 holds next line from file1; @f1 holds fields of $line1
($line1,@f1) = getline($file1,\$n1);
# $line2 holds next line from file2; @f2 holds fields of $line2
($line2,@f2) = getline($file2,\$n2);

$aref = \@file1;



# Loop is over all lines in file2:
while ($line2) {

  $seqname2 = $f2[0];

  # printbuffer1 holds lines to be printed from file1. output these now, if it's time.
  while (@printbuffer1) {
    my @f = split /\t/, $printbuffer1[0];
    if ($f[0] lt $seqname2 or ($f[0] eq $seqname2 and $f[3] <= $f2[3])) {
      print shift @printbuffer1;
    } else {
      last;
    }
  }

  if ($seqname2 ne $lastseqname2) {
    @file1 = ();
    $lastseqname2 = $seqname2;
  }

  # Keep sucking lines from file1 and sticking GFF coords into @file1,
  # until past $line2 endpoint. Entries are "start end linenumber"
  #
  while ($line1 and $seqname1 = $f1[0], ($seqname1 lt $seqname2 or ($seqname1 eq $seqname2 and $f1[3] - $near <= $f2[4]))) {

    if ($seqname1 eq $seqname2) {
      push @file1, "@f1[3,4] $n1";
    }

    # put lines from file1 into file1 printing buffer if -print1 switch is set
    if ($print1) {
      if ($seqname1 lt $seqname2 or ($seqname1 eq $seqname2 and $f1[3] <= $f2[3])) {
        print $line1;
      } else {
        push @printbuffer1, $line1;
      }
    }
    ($line1,@f1) = getline($file1,\$n1);
  }



# use Data::Dumper;
# print Data::Dumper->Dump([$aref]);
# die;



  # loop over every entry in @$aref until past $line2 endpoint $skip
  # and $hitstart2yet are used to keep track of entries in @$aref that
  # may be discarded.
  #
  @lines = ();
  $hitstart2yet = $skip = 0;
  ($start2,$end2) = @f2[3,4];
  for ($i=0;$i<@$aref;$i++) {
    $entry1 = $aref->[$i];
    last if $entry1 - $near > $end2;
    ($start1,$end1,$n) = split /\s+/, $entry1;

    # check for entries that may be discarded
    #
    if ($end1 + $near < $start2) {
      if (!$hitstart2yet) {
        $skip = $i + 1;
      }
    } else {

      # if line from file1 overlaps line from file2, store linenumber in @lines array
      #
      $hitstart2yet = 1;
      ($effstart1,$effend1) = ($start1-$near/2,$end1+$near/2);
      ($effstart2,$effend2) = ($start2-$near/2,$end2+$near/2);
      $maxstart = $effstart1 > $effstart2 ? $effstart1 : $effstart2;
      $minend = $effend1 < $effend2 ? $effend1 : $effend2;
      $overlaplen = $minend + 1 - $maxstart;
      $len1 = $effend1 + 1 - $effstart1;
      $len2 = $effend2 + 1 - $effstart2;
      if (($len1 && $overlaplen/$len1 >= $minfrac1) && ($len2 && $overlaplen/$len2 >= $minfrac2)) {
        push @lines, $n;
      }
    }
  }
  while ($skip-- > 0) {
    shift @$aref;
  }

  # print if eligible
  #
  if ($maxhits != 1) {
    $test = (@lines >= $minhits && @lines <= $maxhits);
  } else {
    $test = (@lines >= $minhits);
  }
  if ($test && ($all || !$not)) {
    if ($quiet) {
      print $line2;
    } else {
      chomp $line2;
      print "$line2 intersect($name1)=(@{[sort {$a<=>$b} @lines]})\n";
    }
  } elsif (!$test && ($all || $not)) {
    print $line2;
  }

  ($line2,@f2) = getline($file2,\$n2);
}

print $line1 if $print1 && $line1;


sub getline {
  my ($file, $nref) = @_;
  my $line = '';

  while (@$file and !$line || $line !~ /\S/) {
    $line = shift @$file;
  }

  $$nref = $file1{$line};

  my @f = split(/\t/,$line,9);

  # Hack to make the script work with the parent info instead of the info in the main fields:
  if ($byparent) {
    my @p = split " ", $f[8];
    ($f[0], $f[1], $f[2], $f[3], $f[4], $f[5], $f[6], $f[7], $f[8]) = 
      ($p[1], '.', '.', $p[2], $p[3], '.', '+', '.', '.');
  }
  if (@f == 0) {
    return (undef);
  }
  if ($regex) {
    $f[0] =~ /($regex)/;
    $f[0] = $1;
  }
  if (!$shutit) {
    if (@f < 9) {
      warn "Warning: fewer than 9 fields at $file line $$nref\n";
    }
    if (join("",@f[0..7]) =~ / /) {
      warn "Warning: space in tab-delimited field at $file line $$nref\n";
    }
  }
  ($line,@f);
}

sub generic {
  my $sorted = 0;
  my @a = split /\t/, $a;
  my @b = split /\t/, $b;
  for (my $i = 0; $i<@a && !$sorted ; $i++) {
    if ($a[$i] =~ /^\d+$/ && $b[$i] =~ /^\d+$/) {
      $sorted = $a[$i] <=> $b[$i];
    } else {
      $sorted = $a[$i] cmp $b[$i];
    }
  }
  $sorted;
}

=head1 SYNOPSIS:

 gffintersect.pl [OPTIONS] file1 file2

=head1 DESCRIPTION:

Finds intersection of two GFF files. Prints lines from file2 that
intersect (or do not intersect) with file1

=head1 OPTIONS:

=over 4

=item --not

To print lines from file2 that don't match the criteria rather than lines that do

=item --all

To print everything from file2

=item --self

for comparisons vs self (omit file2); also adds \"self&\" prefix to name1

=item --minfrac1 xxx

To specify minimum fractional overlap for file1 features

=item --minfrac2 xxx

To specify minimum fractional overlap for file2 features

=item --near xxx

To extend definition of \"overlap\" to nearby GFFs

=item --minhits xxx

To specify minimum number of hits required to
print lines from file2 (default is 1)

=item --maxhits xxx

To specify maximum number of hits required to
print lines from file2 (default is 1)

=item --name1

If file1 is a name that won't look good in \"intersect(..)=(....)\" messages

=item --quiet

To suppress \"intersect(..)=(....)\" messages altogether

=item --print1

To copy all lines from file1 to output

=item --regex

To specify a (perl) regular expression determining how
much of the reference id in field one that is to be taken
as reference id and how much that is just an exon, intron
or alt.spl (Eg. --regex '^[A-Za-z0-9]+_[A-Za-z0-9]+''
specifying that only NC_003279 of NC_003279_3 is to
determine the reference sequence.)

=item --shutit

To suppress GFF format warnings

=item --presorted

The GFF files are already sorted.

=back

=head1 EXAMPLES:



=head1 AUTHOR:

Kasper Munch

=head1 COPYRIGHT:

This program is free software. You may copy and redistribute it under
the same terms as Perl itself.


=cut
