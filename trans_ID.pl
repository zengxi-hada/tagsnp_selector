#!/usr/bin/perl

use strict;
use warnings;


my $frq = shift;
my $tagsnp = shift;
my $out = shift;

die "perl $0 <FRQ file> <tagsnp file> <out put file>\n" unless ($frq && $tagsnp && $out);

my %frq;
open FRQ, $frq or die $!;
while(<FRQ>){
    chomp;
    my @a = split;
    $frq{$a[1]} = $a[0];
}
close FRQ;

open SNP, $tagsnp or die $!;
open OUT, ">$out" or die $!;

while(<SNP>){
	chomp;
	my @a = split;
	print OUT $frq{$a[0]}, "\n" if exists $frq{$a[0]};
}

close SNP;  
close OUT;
