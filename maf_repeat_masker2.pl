#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;


=head1 function

   This application is designed to filter the frq file  by maf and cut the sequence 

=head1 usage

  perl $0
		-frq	<str>		input file
		-maf1	<num>		the threhold of maf [0.02]
		-maf2	<num>		the max threhold of maf [0.45]
		-l	<int>		length of each probe [91]
		-f	<str>		reference fasta file
		-o	<str>		output file (sequence)
		-re	<str>		repeat masker file
		-chr	<str>		name of chromosome
		-h|-?			help
=head1 version

  v1.0  2011-11-01

=head1 author

 zengxi@genomics.org.cn

=cut

my ($frq, $out, $help, $maf1, $length, $fa, $maf2, $chr, $repeat);
$maf1 = 0.02;
$maf2 = 0.45;
$length = 91;


GetOptions(
    'frq=s' => \$frq,
    'o=s' => \$out,
    'maf1=f' => \$maf1,
	'maf2=f' => \$maf2,
	'l=i' => \$length,
	'f=s' => \$fa,
	're=s' => \$repeat,
	'chr=s' => \$chr,
    'h|?' => \$help,
);

die `pod2text $0` unless ($frq && $out && $fa && $repeat && $chr);
die `pod2text $0` if $help;

my %h;

my %re; my $chr1;

if($repeat =~ /gz$/){
    open RE, "<:gzip", $repeat or die $!;
}else{
    open RE, $repeat or die $!;
}

while(<RE>){
    chomp;
    my @a = split;
    $chr1 = $a[5];
    next if $chr1 ne $chr;
    $re{$chr1}{"re\t$a[7]"} = $a[6];
}
close RE;


my @right_offset = keys %{$re{$chr}};

my ($tmp_right_offset);

my %ppos;

open FRQ, $frq or die $!;
while(<FRQ>){
	chomp;
	my @a = split;
	next if ($a[-1] < $maf1 || $a[-1] > $maf2);
	my $pos = $a[1];
	$ppos{"frq\t$pos"} = 1;
#	($zero, $tmp_pos) = &BinarySearch(\@right_offset, $pos);
#	$tmp_right_offset = $right_offset[$tmp_pos];
#    if($re{$chr}{$tmp_right_offset} > $pos){
# 	   next;
#	}	
#	$h{$a[1]} = $a[2];
}
close FRQ;

my @new_array = sort {(split /\s+/, $a)[1] <=> (split /\s+/, $b)[1]} (@right_offset, (keys %ppos));

for my $i(0..$#new_array){
	my $pos;
	if(exists $ppos{$new_array[$i]}){
		my $count = 0;
		while(1){
			$count++;
			if(exists $re{$chr}{$new_array[$i+$count]}){	
				$tmp_right_offset = $new_array[$i + $count];
				last;
			}
		}	
		$pos = (split /\s+/, $new_array[$i])[1];
		delete $ppos{$new_array[$i]} if($re{$chr}{$tmp_right_offset} <= $pos);
	}
}	

my $chr2;

open FA, $fa or die $!;
$/ = '>';
<FA>;
while(<FA>){
	chomp;
	($chr2)= $_=~ /^(chr\w+)\n/;
	s/\n//g;
	s/$chr2//;	
	$fa = $_;
}
close FA;

my $read;
my $pos;
my $count = 0;

open OUT, ">$out" or die $!;
for my $i(sort {(split /\s+/, $a)[1]<=>(split /\s+/, $b)[1]} keys %ppos){
	$pos = (split /\s+/, $i)[1];
	$count++;
	my $length1 = $length - 1;
	$read = substr($fa,$pos - $length1/2 - 1,$length1+1);
	my $start = $pos - $length1/2;
	my $end = $pos + $length1/2;
	print OUT ">$chr:$start..$end\n";
	print OUT "$read\n";
}
close OUT;	



################################# Subroutine ###################################


sub BinarySearch  ########### retrun right offset ##########
{
        my ($array1, $value) = @_;
        my @array = @$array1;
        my $middle = 0;
        my $left = 0;
        my $right = @array - 1;
        while($left <= $right)
        {
                $middle = ($left + $right) / 2;
                $middle = int($middle);
                if($array[$middle] == $value)
                {
                        return 1, $middle;
                }
                elsif($array[$middle] < $value)
                {
                        $left = $middle + 1;
                }
                else
                {
                        $right = $middle - 1;
                }
        }
        return 0, $left;
}

