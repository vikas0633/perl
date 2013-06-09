#
#  FASTA splitter  -  a script for partitioning a FASTA file into pieces
#
#  Version 0.1.1 (March 2, 2012)
#
#  Copyright (c) 2012 Kirill Kryukov
#
#  This software is provided 'as-is', without any express or implied
#  warranty. In no event will the authors be held liable for any damages
#  arising from the use of this software.
#
#  Permission is granted to anyone to use this software for any purpose,
#  including commercial applications, and to alter it and redistribute it
#  freely, subject to the following restrictions:
#
#  1. The origin of this software must not be misrepresented; you must not
#     claim that you wrote the original software. If you use this software
#     in a product, an acknowledgment in the product documentation would be
#     appreciated but is not required.
#  2. Altered source versions must be plainly marked as such, and must not be
#     misrepresented as being the original software.
#  3. This notice may not be removed or altered from any source distribution.
#

use strict;

my $n_parts_total = 0;
my $n_parts_sequence = 0;
my $part_total_size = 0;
my $part_sequence_size = 0;
my $eol = "\x0A";
my $line_length = 60;
my $infile = '';

my $usage = q~Usage: perl fasta_splitter.pl method [options] fastafile
Where "method" is one of:
    -n-parts-total N       - split into N parts of similar total size
    -n-parts-sequence N    - split into N parts of similar sequence size
    -part-total-size N     - split into parts of at most N bytes each
    -part-sequence-size N  - split into parts containing at most N bp each
Options:
    -line-length N         - output line lenght, 60 by default
    -eol [dos|mac|unix]    - end-of-line of the output, unix by default
~;

my $argc = scalar(@ARGV);
if (!$argc) { print STDERR $usage; exit; }
my $n_methods = 0;
for (my $i=0; $i<$argc; $i++)
{
    if ($i < $argc-1)
    {
        if ($ARGV[$i] eq '-n-parts-total'     ) { $i++; $n_parts_total      = int($ARGV[$i]); if ($n_parts_total      < 0) { die; } $n_methods++; next; }
        if ($ARGV[$i] eq '-n-parts-sequence'  ) { $i++; $n_parts_sequence   = int($ARGV[$i]); if ($n_parts_sequence   < 0) { die; } $n_methods++; next; }
        if ($ARGV[$i] eq '-part-total-size'   ) { $i++; $part_total_size    = int($ARGV[$i]); if ($part_total_size    < 0) { die; } $n_methods++; next; }
        if ($ARGV[$i] eq '-part-sequence-size') { $i++; $part_sequence_size = int($ARGV[$i]); if ($part_sequence_size < 0) { die; } $n_methods++; next; }
        if ($ARGV[$i] eq '-line-length')        { $i++; $line_length        = int($ARGV[$i]); if ($line_length        < 1) { die; }               next; }
        if ($ARGV[$i] eq '-eol')
        {
            $i++;
            if (lc($ARGV[$i]) eq 'dos') { $eol = "\x0D\x0A"; }
            elsif (lc($ARGV[$i]) eq 'mac') { $eol = "\x0D"; }
            elsif (lc($ARGV[$i]) eq 'unix') { $eol = "\x0A"; }
            else { die "Error: Unknown -eol parameter\n"; }
            next;
        }
    }
    if (!$infile) { $infile = $ARGV[$i]; }
}
if (!$infile) { die "Error: Input file not specified\n"; }
if (!$n_methods) { die "Error: Method not specified\n"; }
if ($n_methods != 1) { die "Error: Multiple methods specified\n"; }
if (!-e $infile) { die "Error: Can't find the input file\n"; }

my $eol_length = length($eol);
my $base = $infile;
my $ext = '';
if ($base =~ /^(.+?)(\.(fasta|faa|fna|fa))$/) { ($base,$ext) = ($1,$2); }


if ($n_parts_total or $n_parts_sequence)
{
    my @names = ();
    my %seq = ();
    my $name = '';

    open (my $IN,'<',$infile) or die "Error: Can't open the input file\n";
    while (<$IN>)
    {
        chomp;
        if (/^>(.*)$/) { $name = $1; push @names, $name; next; }
        $seq{$name} .= $_;
    }
    close $IN;

    my $n_seq = scalar keys %seq;
    my ($total_sequence_size,$total_data_size) = (0,0);
    foreach my $name (keys %seq)
    {
        $total_sequence_size += length $seq{$name};
        $total_data_size += length($name) + $eol_length + 1;
        while ($seq{$name} =~ /(.{1,$line_length})/g) { $total_data_size += length($1) + $eol_length; }
    }
    print STDERR "$n_seq sequences, $total_sequence_size bp, $total_data_size bytes\n";

    my $i = 0;
    my $written = 0;
    my $max_part = $n_parts_total ? $n_parts_total : $n_parts_sequence;
    my $len = length($max_part);
    for (my $part = 1; $part <= $max_part; $part++)
    {
        my $should_write = sprintf("%.0f",$n_parts_total ? ($total_data_size * $part / $n_parts_total) : ($total_sequence_size * $part / $n_parts_sequence) );
        my $part_file = $base;
        if ($part_file !~ /\.part-\d+$/) { $part_file .= '.part'; }
        $part_file .= sprintf("-%0*d%s",$len,$part,$ext);
        open(my $OUT,'>',$part_file) or die "Error: Can't create file \"$part_file\"\n";
        binmode $OUT;
        while ($written < $should_write)
        {
            print $OUT '>', $names[$i], $eol;
            if ($n_parts_total) { $written += length($names[$i]) + $eol_length + 1; }
            while ($seq{$names[$i]} =~ /(.{1,$line_length})/g)
            {
                print $OUT $1, $eol;
                if ($n_parts_total) { $written += length($1) + $eol_length; }
            }
            if ($n_parts_sequence) { $written += length($seq{$names[$i]}); }
            $i++;
        }
        close $OUT;
    }
}
else
{
    my $in_size = -s $infile;
    open (my $IN,'<',$infile) or die "Error: Can't open the input file\n";
    my $name = '';
    my $data = '';
    my $part = 1;
    my $written_this_part = 0;
    my $max_part_size = $part_total_size ? $part_total_size : $part_sequence_size;
    my $len = length(int($in_size/$max_part_size));  # This is just a guess.
    my $out_file = sprintf("%s.part-%0*d%s",$base,$len,$part,$ext);
    open (my $OUT,'>',$out_file) or die "Can't create output file \"$out_file\"\n";
    binmode $OUT;
    while (<$IN>)
    {
        chomp;
        if (/^>(.*)$/)
        {
            my $new_name = $1;
            if ($name)
            {
                my $seq_size = $part_total_size ? (length($name) + $eol_length + 1) : length($data);
                if ($part_total_size) { while ($data =~ /(.{1,$line_length})/g) { $seq_size += length($1) + $eol_length; } }
                if ( ($written_this_part + $seq_size > $max_part_size) and ($written_this_part != 0) )
                {
                    close $OUT;
                    $part++;
                    $out_file = sprintf("%s.part-%0*d%s",$base,$len,$part,$ext);
                    open ($OUT,'>',$out_file) or die "Can't create output file \"$out_file\"\n";
                    binmode $OUT;
                    $written_this_part = 0;
                }
                dump_seq($OUT,$name,$data);
                $written_this_part += $seq_size;
            }
            $name = $new_name;
            $data = '';
            next;
        }
        $data .= $_;
    }
    if ($name)
    {
        my $seq_size = $part_total_size ? (length($name) + $eol_length + 1) : length($data);
        if ($part_total_size) { while ($data =~ /(.{1,$line_length})/g) { $seq_size += length($1) + $eol_length; } }
        if ( ($written_this_part + $seq_size > $max_part_size) and ($written_this_part != 0) )
        {
            close $OUT;
            $part++;
            $out_file = sprintf("%s.part-%0*d%s",$base,$len,$part,$ext);
            open ($OUT,'>',$out_file) or die "Can't create output file \"$out_file\"\n";
            binmode $OUT;
        }
        dump_seq($OUT,$name,$data);
    }
    close $IN;
    close $OUT;
}


sub dump_seq
{
    my ($file,$name,$data) = @_;
    if (!$name) { return; }
    print $file '>', $name, $eol;
    while ($data =~ /(.{1,$line_length})/g)
    {
        print $file $1, $eol;
    }
}
