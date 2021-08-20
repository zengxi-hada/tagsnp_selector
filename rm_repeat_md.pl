#!/usr/bin/perl

use strict;
use warnings;
use PerlIO::gzip;
use Getopt::Long;
use File::Basename;
use Cwd qw/abs_path/;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;

=head1 function
  
	Filter the set of Polymorphism points according to the result of Blast(only select unique mapping reads)

=head1 usage
  	perl $0	
			-b <str>	the result of blast
			-o <str>	the out put file
			-c <str>	the chr
			-m <int>	the max length of mismatch in the process of blast [5]		
			-h|-?		help

=head1 author
  zengxi@genomics.org.cn

=head1 version

  v1.0	2011-10-15

=cut

my ($blast, $chr, $repeat, $out, $help,$length, $max_mismatch);

$max_mismatch = 5;

GetOptions(
	'b=s' => \$blast,
	'o=s' => \$out,
	'c=s' => \$chr,
	'm=i' => \$max_mismatch,
	'h|?' => \$help,
);

die `pod2text $0` unless ($blast && $out && $chr);
die `pod2text $0` if $help;

my %re; my $chr1;

#if($repeat =~ /gz$/){
#	open RE, "<:gzip", $repeat or die $!;
#}else{
#	open RE, $repeat or die $!;
#}

#while(<RE>){
#	chomp;
#	my @a = split;
#	$chr1 = $a[5];
#	next if $chr1 ne $chr; 
#	$re{$chr1}{$a[7]} = $a[6];
#}
#close RE;

#my @right_offset = sort {$a<=>$b} keys %{$re{$chr}};

my @tagsnp;
my ($start, $pos, $end, $chr2, $zero, $tmp_pos, $tmp_right_offset);
my @a2; my $flag = 0;

open TAG, ">$out" or die $!;

open BLAST, $blast or die $!;
$/ = "# BLASTN 2.2.21 [Jun-14-2009]\n";
<BLAST>;
while(<BLAST>){
	chomp;
	s/# Query.*\n//;
	s/# Database.*\n//;
	s/# Fields.*\n//;
	my @a0 = split /\n/, $_;
	my @tmp = ();
	for my $i(@a0){
		my @a = split /\s+/, $i;		
		@a2 = @a;
		($start) = ($a[0] =~ /^chr\d+:(\d+).*\d+$/);
		($end) = ($a[0] =~ /\d+\.\.(\d+)/);
        ($chr2) = ($a[0] =~ /(chr\d+)/);
		$length =  $end - $start + 1;
		my $mismatch = $a[4] + $a[6] - 1 + $length - $a[7]; 
		my $identity = sprintf "%.0f", $a[2]*0.01*$a[3];
		my $gap_length = $a[3] - $a[4] - $identity; 
#		next if ($gap_length > 5);
#		next if ($mismatch > 5);
		next if ($gap_length+$mismatch > $max_mismatch);
		next if ($a[1] ne $chr2);
		$pos = $start + ($length-1)/2;
#        ($zero, $tmp_pos) = &BinarySearch(\@right_offset, $pos);
#		push @tmp, $i;
#		$tmp_right_offset = $right_offset[$tmp_pos];
#        if($re{$chr2}{$tmp_right_offset} > $pos){
#            push @tagsnp, $pos;
#        }
		push @tagsnp, $pos;
#		$flag = 0;
#		if (@tmp > 1){$flag = 1; last;}
		if (@tagsnp > 1){$flag = 1; last;}
	}
	if($flag == 0 and @tagsnp == 1){
		print TAG "$tagsnp[0]\n";
	}
	$flag = 0;		
#	}
	@tagsnp = ();
		
}
close BLAST;
close TAG;



sub BinarySearch
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

