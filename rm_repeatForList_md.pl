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
			-b <str>	the result file list of blast
			-o <str>	the out put file
			-m <int>    the max length of mismatch in the process of blast [5]				
			-c <str>	the chr
			-h|-?		help

=head1 author
  zengxi@genomics.org.cn

=head1 version

  v1.0	2011-10-15

=cut

my ($list, $chr, $repeat, $out, $help,$length, $blast, $max_mismatch);

$max_mismatch = 5;

GetOptions(
	'b=s' => \$list,
	'o=s' => \$out,
	'm=i' => \$max_mismatch,
	'c=s' => \$chr,
	'h|?' => \$help,
);

die `pod2text $0` unless ($list && $out && $chr);
die `pod2text $0` if $help;

my %re; my $chr1;

my @tagsnp;
my ($start, $pos, $end, $chr2, $zero, $tmp_pos, $tmp_right_offset);
my @a2; my $flag = 0;

open TAG, ">$out" or die $!;

open LIST, $list or die $!;
while(<LIST>){ 
	chomp;
	$blast = $_;
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
#			next if ($gap_length > 5);
#			next if ($mismatch > 5);
			next if ($gap_length+$mismatch > $max_mismatch);
			next if ($a[1] ne $chr2);
			$pos = $start + ($length-1)/2;
			push @tagsnp, $pos;
			if (@tagsnp > 1){$flag = 1; last;}
		}
		if($flag == 0 and @tagsnp == 1){
			print TAG "$tagsnp[0]\n";
		}
		$flag = 0;		
		@tagsnp = ();
			
	}
	close BLAST;
	$/="\n";
}

close TAG;
close LIST;


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

