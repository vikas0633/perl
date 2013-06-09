#!/usr/bin/perl -w

# WORKING COPY
# 11/12/2008

# THIS COPIED DIRECTLY FROM
# http://gmod.org/wiki/Subtrack_HOWTO

use strict;

my $name = shift;
my $desc = shift;

@ARGV or die "I need three args: name, desc, filename";

print qq(track type=wiggle_0 name="$name" description="$desc"\n);
while (<>) {
  chomp;
  my ($ref,$start,$end,$score) = (split)[0,3,4,5];
  $ref =~ s/CHROMOSOME_|chr//;
  $start--; # zero-based, half-open                 
  $score = 255 if $score !~ /^[-.Ee0-9]$/;
                                                                                                                     
  print join("\t",$ref,$start,$end,$score), "\n";
}
