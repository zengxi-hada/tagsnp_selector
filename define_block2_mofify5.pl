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

   This application is designed to define blocks that containing tag snps

=head1 usage

  perl $0
		-ld	<str>			input LD file
		-frq	<str>			input frq file
		-tagsnp	<str>			The alternative collection of tagsnp
		-o	<str>			output file
		-minr	<float>			min value of R squre [0.8]
		-rep	<float>			min value of represent rate [0.85]
		-stp	<float>			step in the process of defining block [0.01]
		-h | -?				help

=head1 author

  zengxi@genomics.org.cn

=head1 version

  2011-10-15: v1.00

=cut

my ($ld, $frq, $out, $help, $tagsnp);
my $min_rsqure=0.8; my $rep_rate=0.85; my $step=0.01;

GetOptions(
    'ld=s'=>\$ld,
	'frq=s'=>\$frq,
    'o=s'=>\$out,
    'tagsnp=s'=>\$tagsnp,
	'minr=f'=>\$min_rsqure,
	'rep=f'=>\$rep_rate,
	'stp=f'=>\$step,
    'h|?' =>\$help,
);

die `pod2text $0` if(!$ld || !$frq || !$out || !$tagsnp);
die `pod2text $0` if($help);

my %tag;
my %itag;
my %ld;
my $chr;
my %block;

if($ld=~/chr\d+/){
	$chr = $&;
}

my $tmp1; my $tmp2; my $tmp3; my $tmp_id3;

open TAG, "$tagsnp" or die $!;
while(<TAG>){
	chomp;
	my @tag = split;
	$itag{$tag[0]} = 1;  
}
close TAG;

open FRQ, $frq or die $!;

my $tmp_frq; my $tmp_maf; my $tmp; my $tmp_id;
my @frq = ();
my $count_frq = 0;
my $counta = 0;
my ($start, $end, $start1, $start2, $start3, $end1, $end2, $end3);

while(<FRQ>){
	chomp;
	my @a_a = split;
	
	$count_frq++;
	if($count_frq == 1){
		$tmp_frq = $a_a[1];
        $tmp_maf = $a_a[2];
		$tmp = $_;
		next;
	}
	if($a_a[1]-$tmp_frq > 60){
		push @frq, $tmp;
		$tmp_frq = $a_a[1];
        $tmp_maf = $a_a[2];
        $tmp = $_;
		
		$counta++;
		if(scalar @frq == 3){
			my (@a, @a2, @a3);			
			my $line1 = shift @frq;
			@a = split /\s+/, $line1;				
			$tag{$a[0]}{"info"} = "$chr\_$a[1]\_$a[2]";
			$block{$a[1]} = $a[0] if exists $itag{$a[0]};
##			$tag{$a[0]}{"block"} = $a[1];
										
			my $line2 = shift @frq;
			@a2 = split /\s+/, $line2;
			$tag{$a2[0]}{"info"} = "$chr\_$a2[1]\_$a2[2]";
			$block{$a2[1]} = $a2[0] if exists $itag{$a2[0]};
##     	    $tag{$a2[0]}{"block"} = $a2[1];	

			my $line3 = shift @frq;
			@a3 = split /\s+/, $line3;
	        $tag{$a3[0]}{"info"} = "$chr\_$a3[1]\_$a3[2]";
			$block{$a3[1]} = $a3[0] if exists $itag{$a3[0]};
##       	$tag{$a3[0]}{"block"} = $a3[1];
			if($counta == 3){
		        $start = int(($a[1] + $a2[1])/2);
   			    $end = int(($a2[1] + $a3[1])/2);
 	   		    $tag{$a2[0]}{"region"} = "$start\-$end";
	        	my $tmp_tmp = 0;
    	    	$tag{$a[0]}{"region"} = "$tmp_tmp"."-"."$start";
				
			}elsif($counta > 3){
				$start1 =int( ($tmp2 + $tmp3) / 2);
	            $end1 = int(($tmp3 + $a[1]) / 2);
	            $tag{$tmp_id3}{"region"} = "$start1\-$end1";
				$start2 = $end1;
				$end2 = int(($a[1] + $a2[1]) / 2);
				$tag{$a[0]}{"region"} = "$start2\-$end2";
				$start3 = $end2;
				$end3 = int(($a2[1] + $a3[1]) / 2);
				$tag{$a2[0]}{"region"} = "$start3\-$end3";
			}
			@frq = ();
			$tmp_id3 = $a3[0];
            $tmp1 = $a[1];
            $tmp2 = $a2[1];
            $tmp3 = $a3[1];
		}
	}elsif($a_a[1]-$tmp_frq <= 60){
		if($a_a[2] < $tmp_maf){
			$tmp_frq = $a_a[1];
			$tmp_maf = $a_a[2];
			$tmp = $_;
		}
	}
	$tmp_id = $a_a[0];
}
close FRQ;

$tag{$tmp_id}{"info"} = "$chr\_$tmp_frq\_$tmp_maf";
$tag{$tmp_id}{"region"} = "$tmp_frq\-end of the choromosome";

if(@frq == 1){
	my $line1 = shift @frq;
	my @a = split /\s+/, $line1;	
	my $startp1 = int( ($tmp2 + $tmp3) / 2);
	my $endp1 = int(($tmp3 + $a[1]) / 2);
	$tag{$tmp_id3}{"region"} = "$startp1\-$endp1";
	$tag{$a[0]}{"region"} = "$endp1\-end of the choromosome";
	$tag{$a[0]}{"info"} = "$chr\_$a[1]\_$a[2]";	
}elsif(@frq == 2){
	my $line1 = shift @frq;
	my @a = split /\s+/, $line1;
	$tag{$a[0]}{"info"} = "$chr\_$a[1]\_$a[2]";
	my $line2 = shift @frq;
	my @a2 = split /\s+/, $line1;
	$tag{$a2[0]}{"info"} = "$chr\_$a2[1]\_$a2[2]";
	my $startp1 = int( ($tmp2 + $tmp3) / 2);
	my $endp1 = int(($tmp3 + $a[1]) / 2);
	$tag{$tmp_id3}{"region"} = "$startp1\-$endp1";
	my $startp1_1 = $endp1;
	my $endp1_1 = int(($a[1]+$a2[1])/2);
	$tag{$a[0]}{"region"} = "$startp1_1\-$endp1_1";	
	$tag{$a2[0]}{"region"} = "$endp1_1\-end of the choromosome";
}elsif(@frq == 0){
	$tag{$tmp_id3}{"region"} = "$end3\-end of the choromosome";
} 

@frq = ();


open LD, $ld or die $!;
while(<LD>){
	chomp;
	next if (/L1/);
	my @a = split;
	next if ($a[2] < 0.8);
	next if (not exists $tag{$a[0]}); 
	next if (not exists $tag{$a[1]}); 
	$tag{$a[0]}{"R_squre"}{$a[1]} = $a[2] if(exists $itag{$a[0]});
	$tag{$a[1]}{"R_squre"}{$a[0]} = $a[2] if(exists $itag{$a[1]});

}
close LD;



########################## initiate blocks when the R_squre is 1 ###########################

my @tmp_array = ();
my $count_loop = 0;
my ($flag1_1, $flag2_1, $flag3_1) = (0, 0, 0);
my ($flag1, $flag2) = (0, 0);
my ($tmp_id1_1, $tmp_id2_1, $tmp_id3_1, $tmp_pos1, $tmp_pos2, $tmp_pos3);

for my $i(sort {$a<=>$b} keys %block){
	$count_loop++;
	push @tmp_array, $i;
	if(@tmp_array == 3){
		my $id1 = $block{$tmp_array[0]};
		my $id2 = $block{$tmp_array[1]};
		my $id3 = $block{$tmp_array[2]};
		if($count_loop == 3){	
			if(exists $tag{$id2}{"R_squre"}{$id1} && exists $tag{$id2}{"R_squre"}{$id3} && $tag{$id2}{"R_squre"}{$id1} == 1 && $tag{$id2}{"R_squre"}{$id3} == 1 && exists $itag{$id2}){
				$block{$tmp_array[0]} = join ' ', ($id1, $id2, $id3);
				delete $block{$tmp_array[1]};
				delete $block{$tmp_array[2]};
				$flag3_1 = 1;
			}elsif(exists $tag{$id2}{"R_squre"}{$id1} && $tag{$id2}{"R_squre"}{$id1} == 1 && ((exists $tag{$id2}{"R_squre"}{$id3} && $tag{$id2}{"R_squre"}{$id3} != 1) || not exists $tag{$id2}{"R_squre"}{$id3}) && (exists $itag{$id2} || exists $itag{$id1})){
				$block{$tmp_array[0]} = join ' ', ($id1, $id2);
				delete $block{$tmp_array[1]};
				$tmp_id1_1 = $id3;
				$tmp_pos1 = $tmp_array[2];
				$flag1_1 = 1;
			}elsif((exists $tag{$id2}{"R_squre"}{$id1} && $tag{$id2}{"R_squre"}{$id1} != 1) || not exists $tag{$id2}{"R_squre"}{$id1}){
				$tmp_id2_1 = $id2;
				$tmp_id3_1 = $id3;
				$tmp_pos2 = $tmp_array[1];
				$tmp_pos3 = $tmp_array[2];
				$flag2_1 = 1;
			}
		}else{
			if($flag1_1 == 1){
				$flag1_1 = 0;
				if(exists $tag{$id1}{"R_squre"}{$id2} && $tag{$id1}{"R_squre"}{$id2} == 1 && exists $tag{$id1}{"R_squre"}{$tmp_id1_1} && $tag{$id1}{"R_squre"}{$tmp_id1_1} == 1 && exists $itag{$id1}){	
					$block{$tmp_pos1} = join ' ', ($tmp_id1_1, $id1, $id2);
					delete $block{$tmp_array[0]};
					delete $block{$tmp_array[1]};
					$flag1_1 = 1;
					$tmp_id1_1 = $id3;
					$tmp_pos1 = $tmp_array[2];
				}elsif(((exists $tag{$id1}{"R_squre"}{$id2} && $tag{$id1}{"R_squre"}{$id2} != 1) || not exists $tag{$id1}{"R_squre"}{$id2}) && exists $tag{$id1}{"R_squre"}{$tmp_id1_1} && $tag{$id1}{"R_squre"}{$tmp_id1_1} == 1 && (exists $itag{$tmp_id1_1} || exists $itag{$id1})){
					$block{$tmp_pos1} = join ' ', ($tmp_id1_1, $id1);
					delete $block{$tmp_pos1};
					$tmp_id2_1 = $id2;
					$tmp_id3_1 = $id3;
					$tmp_pos2 = $tmp_array[1];
					$tmp_pos3 = $tmp_array[2];
					$flag2_1 = 1;
				}elsif((exists $tag{$id1}{"R_squre"}{$tmp_id1_1} && $tag{$id1}{"R_squre"}{$tmp_id1_1} != 1) || not exists $tag{$id1}{"R_squre"}{$tmp_id1_1}){
					if(exists $tag{$id2}{"R_squre"}{$id1} && $tag{$id2}{"R_squre"}{$id1} == 1 && exists $tag{$id2}{"R_squre"}{$id3} && $tag{$id2}{"R_squre"}{$id3} == 1 && exists $itag{$id2}){
						$block{$tmp_array[0]} = join ' ', ($id1, $id2, $id3);
						delete $block{$tmp_array[1]};
		                delete $block{$tmp_array[2]};
						$flag3_1 = 1;
					}elsif(exists $tag{$id2}{"R_squre"}{$id1} && $tag{$id2}{"R_squre"}{$id1} == 1 && ((exists $tag{$id2}{"R_squre"}{$id3} && $tag{$id2}{"R_squre"}{$id3} != 1) || not exists $tag{$id2}{"R_squre"}{$id3}) && (exists $itag{$id2} || exists $itag{$id1})){
						$block{$tmp_array[0]} = join ' ', ($id1, $id2);
						$flag1_1 = 1;
						$tmp_pos1 = $tmp_array[2];
						$tmp_id1_1 = $id3;
						delete $block{$tmp_array[1]};
					}elsif((exists $tag{$id2}{"R_squre"}{$id1} && $tag{$id2}{"R_squre"}{$id1} != 1) || not exists $tag{$id2}{"R_squre"}{$id1}){ 
						$tmp_id2_1 = $id2;
						$tmp_id3_1 = $id3;
						$tmp_pos2 = $tmp_array[1];
						$tmp_pos3 = $tmp_array[2];
						$flag2_1 = 1;
					}
				}
										
			}elsif($flag2_1 == 1){
				$flag2_1 = 0;
				if(exists $tag{$tmp_id3_1}{"R_squre"}{$id1} && $tag{$tmp_id3_1}{"R_squre"}{$id1} == 1){
					if(exists $tag{$tmp_id3_1}{"R_squre"}{$tmp_id2_1} && $tag{$tmp_id3_1}{"R_squre"}{$tmp_id2_1} == 1 && exists $itag{$tmp_id3_1}){
						$block{$tmp_pos2} = join ' ', ($tmp_id3_1, $tmp_id2_1, $id1);	
						delete $block{$tmp_pos3};
						delete $block{$tmp_array[0]};		
						$tmp_id2_1 = $id2;
						$tmp_id3_1 = $id3;
						$tmp_pos2 = $tmp_array[1];
						$tmp_pos3 = $tmp_array[2];
						$flag2_1 = 1;
					}elsif((exists $tag{$tmp_id3_1}{"R_squre"}{$tmp_id2_1} && $tag{$tmp_id3_1}{"R_squre"}{$tmp_id2_1} != 1) || not exists $tag{$tmp_id3_1}{"R_squre"}{$tmp_id2_1}){
						$tmp_id1_1 = $tmp_id3_1;
						$tmp_pos1 = $tmp_pos3;
						if(exists $tag{$id1}{"R_squre"}{$id2} && $tag{$id1}{"R_squre"}{$id2} == 1 && exists $tag{$id1}{"R_squre"}{$tmp_id1_1} && $tag{$id1}{"R_squre"}{$tmp_id1_1} == 1 && exists $itag{$id1}){
		                    $block{$tmp_pos1} = join ' ', ($tmp_id1_1, $id1, $id2);
		                    delete $block{$tmp_array[0]};
		                    delete $block{$tmp_array[1]};
							$flag1_1 = 1;
							$tmp_id1_1 = $id3;
							$tmp_pos1 = $tmp_array[2];	
 		               }elsif(((exists $tag{$id1}{"R_squre"}{$id2} && $tag{$id1}{"R_squre"}{$id2} != 1) || not exists $tag{$id1}{"R_squre"}{$id2}) && exists $tag{$id1}{"R_squre"}{$tmp_id1_1} && $tag{$id1}{"R_squre"}{$tmp_id1_1} == 1 && (exists $itag{$tmp_id1_1} || exists $itag{$id1})){
		                    $block{$tmp_pos1} = join ' ', ($tmp_id1_1, $id1);
		                    delete $block{$tmp_array[0]};
		                    $tmp_id2_1 = $id2;
		                    $tmp_id3_1 = $id3;
		                    $tmp_pos2 = $tmp_array[1];
		                    $tmp_pos3 = $tmp_array[2];
		                    $flag2_1 = 1;
		               }elsif((exists $tag{$id1}{"R_squre"}{$tmp_id1_1} && $tag{$id1}{"R_squre"}{$tmp_id1_1} != 1) || not exists $tag{$id1}{"R_squre"}{$tmp_id1_1}){
		                    if(exists $tag{$id2}{"R_squre"}{$id1} && $tag{$id2}{"R_squre"}{$id1} == 1 && exists $tag{$id2}{"R_squre"}{$id3} && $tag{$id2}{"R_squre"}{$id3} == 1 && exists $itag{$id2}){
		                        $block{$tmp_array[0]} = join ' ', ($id1, $id2, $id3);
		                        delete $block{$tmp_array[1]};
		                        delete $block{$tmp_array[2]};
								$flag3_1 = 1;
		                    }elsif(exists $tag{$id2}{"R_squre"}{$id1} && $tag{$id2}{"R_squre"}{$id1} == 1 && ((exists $tag{$id2}{"R_squre"}{$id3} && $tag{$id2}{"R_squre"}{$id3} != 1) || not exists $tag{$id2}{"R_squre"}{$id3}) && (exists $itag{$id2} || exists $itag{$id1})){
		                        $block{$tmp_array[0]} = join ' ', ($id1, $id2);
		                        $flag1_1 = 1;
								$tmp_id1_1 = $id3;
								$tmp_pos1 = $tmp_array[2];
		                        delete $block{$tmp_array[1]};
		                    }elsif((exists $tag{$id2}{"R_squre"}{$id1} && $tag{$id2}{"R_squre"}{$id1} != 1) || not exists $tag{$id2}{"R_squre"}{$id1}){
		                        $tmp_id2_1 = $id2;
		                        $tmp_id3_1 = $id3;
		                        $tmp_pos2 = $tmp_array[1];
		                        $tmp_pos3 = $tmp_array[2];
		                        $flag2_1 = 1;
		                    }
						}
					}

				}elsif((exists $tag{$tmp_id3_1}{"R_squre"}{$id1} && $tag{$tmp_id3_1}{"R_squre"}{$id1} != 1) || not exists $tag{$tmp_id3_1}{"R_squre"}{$id1}){
					if(exists $tag{$tmp_id3_1}{"R_squre"}{$tmp_id2_1} && $tag{$tmp_id3_1}{"R_squre"}{$tmp_id2_1} == 1){
						$block{$tmp_pos2} = join ' ', ($tmp_id2_1, $tmp_id3_1);
						delete $block{$tmp_pos3};
					}
#################  3 points #################################################
					if(exists $tag{$id2}{"R_squre"}{$id1} && $tag{$id2}{"R_squre"}{$id1} == 1 && exists $tag{$id2}{"R_squre"}{$id3} && $tag{$id2}{"R_squre"}{$id3} == 1 && exists $itag{$id2}){
                        $block{$tmp_array[0]} = join ' ', ($id1, $id2, $id3);
                        delete $block{$tmp_array[1]}; 
                        delete $block{$tmp_array[2]};
						$flag3_1 = 1;
                    }elsif(exists $tag{$id2}{"R_squre"}{$id1} && $tag{$id2}{"R_squre"}{$id1} == 1 && ((exists $tag{$id2}{"R_squre"}{$id3} && $tag{$id2}{"R_squre"}{$id3} != 1) || not exists $tag{$id2}{"R_squre"}{$id3}) && (exists $itag{$id2} || exists $itag{$id1})){
                        $block{$tmp_array[0]} = join ' ', ($id1, $id2);
						$tmp_id1_1 = $id3;
						$tmp_pos1 = $tmp_array[2];
                        $flag1_1 = 1;
                        delete $block{$tmp_array[1]};
                    }elsif((exists $tag{$id2}{"R_squre"}{$id1} && $tag{$id2}{"R_squre"}{$id1} != 1) || not exists $tag{$id2}{"R_squre"}{$id1}){
                        $tmp_id2_1 = $id2;
                        $tmp_id3_1 = $id3; 
                        $tmp_pos2 = $tmp_array[1];
                        $tmp_pos3 = $tmp_array[2];
                        $flag2_1 = 1;
                    }
				}
			}elsif($flag3_1 == 1){
				$flag3_1 = 0;
######################################################################################################
				if(exists $tag{$id2}{"R_squre"}{$id1} && $tag{$id2}{"R_squre"}{$id1} == 1 && exists $tag{$id2}{"R_squre"}{$id3} && $tag{$id2}{"R_squre"}{$id3} == 1 && exists $itag{$id2}){
                	$block{$tmp_array[0]} = join ' ', ($id1, $id2, $id3);
                 	delete $block{$tmp_array[1]};
                	delete $block{$tmp_array[2]};
        	        $flag3_1 = 1;
                }elsif(exists $tag{$id2}{"R_squre"}{$id1} && $tag{$id2}{"R_squre"}{$id1} == 1 && ((exists $tag{$id2}{"R_squre"}{$id3} && $tag{$id2}{"R_squre"}{$id3} != 1) || not exists $tag{$id2}{"R_squre"}{$id3}) && (exists $itag{$id2} || exists $itag{$id1})){
                	$block{$tmp_array[0]} = join ' ', ($id1, $id2);
                	$tmp_id1_1 = $id3;
              		$tmp_pos1 = $tmp_array[2];
               		$flag1_1 = 1;
        	        delete $block{$tmp_array[1]};
                }elsif((exists $tag{$id2}{"R_squre"}{$id1} && $tag{$id2}{"R_squre"}{$id1} != 1) || not exists $tag{$id2}{"R_squre"}{$id1}){
                	$tmp_id2_1 = $id2;
                	$tmp_id3_1 = $id3;
               		$tmp_pos2 = $tmp_array[1];
                	$tmp_pos3 = $tmp_array[2];
               		$flag2_1 = 1;
				}
			}
#########################################################################################################	
		}
			
		@tmp_array = ();
	}
}
							
####################################################################################					
######					cluster with taged point                     ###############
####################################################################################			
@tmp_array = ();

my $bul = 1;

$| = 1;

$min_rsqure = $min_rsqure - 0.000001;
my $count_while = 0;


for (my $t=1; $t>=$min_rsqure; $t=$t-$step){
	my $countipp = 0;
	$count_while = 0;
	my $loop_countp = 0;
	while($bul){
		$count_while++;
		$bul = 0;
		my @keys = sort {$a<=>$b} keys %block;
		my %block_record;
#################### the 1st core loop ##########################
		foreach my $block_posId(0..$#keys){
			if((keys %block_record) > 5){delete $block_record{$block_posId-5};} 
			my %candidate_block;
			my $startip = $block_posId-5;
			my $endip = $block_posId+5;
			if($block_posId-5 < 0){
				$startip = 0;
			}elsif($block_posId+5 >= @keys){
				$endip = @keys - 1;
			}elsif($block_posId-5 < 0 && $block_posId+5 >= @keys){
				$startip = 0;
				$endip = @keys - 1;
			} 
			my @bkArrray_keys = @keys[$startip..$endip];
			my $block_pos = $keys[$block_posId];
			next if not exists $block{$block_pos};
			my @tmp_block_pos2;
			if($block{$block_pos} =~ /\s+/){
				@tmp_block_pos2 = split /\s+/, $block{$block_pos};
			}else{
				@tmp_block_pos2 = ($block{$block_pos});
			}
####################### the 2nd core loop  #########################
			for my $i(@bkArrray_keys){
				my $flag9 = 0;
				my $bk_id = $i;
				my $tmp_combination = join ' ', sort($bk_id, $block_pos);
				for my $ti(1..5){
					if (exists $block_record{$block_posId - $ti}{$tmp_combination}){$flag9 = 1;last;}
				}
				next if  $flag9 == 1;
				$block_record{$block_posId}{$tmp_combination} = 1;
				next if not exists $block{$bk_id};
				my @tmp;
				if($block{$bk_id} =~ /\s+/){
					@tmp = split /\s+/, $block{$bk_id} if exists $block{$bk_id};	
				}else{
					@tmp = ($block{$bk_id}) if exists $block{$bk_id};
				}
				my %rep_num;
				next if not exists $block{$block_pos};
				$loop_countp++; 
				my @tmp_block_pos = sort (@tmp, @tmp_block_pos2);
				my $xh = join ' ', sort @tmp;
				my $zx = join ' ', sort @tmp_block_pos2;
				if ($xh eq $zx){next;}
				my %count_uniq2;
				for my $i2(@tmp_block_pos){
					if(not exists $itag{$i2}){next;}
					$rep_num{$i2} = 0;
					for my $i2_2(@tmp_block_pos){
						if(exists $tag{$i2}{"R_squre"}{$i2_2} && $tag{$i2}{"R_squre"}{$i2_2}>=$t){
							$rep_num{$i2}++;
						}		
					}
				}
				if((keys %rep_num)==0){next;}
				my $max_num = &maxi(values %rep_num); 
				my $iirep_rate = ($max_num+1)/@tmp_block_pos;
				if($iirep_rate > 0.85){
					my $tmp_rate_bkid = "$iirep_rate\-$bk_id";
					$candidate_block{$tmp_rate_bkid} = [@tmp];
#					$block_record{$bk_id} = 1;
				}				
			}
#################################################################################################################
			if((keys %candidate_block) != 0){
				$bul = 1;
				$countipp++;
				for my $candidate_key(sort{(split /-/,$b)[0]<=>(split /-/,$a)[0]} keys %candidate_block){				
					my @tmp_1 = (@tmp_block_pos2, @{$candidate_block{$candidate_key}});
					my %count_uniq;
					my $oo_block_id = (split /-/, $candidate_key)[1];
					my @blockId_array = ($block_pos, $oo_block_id);
					my $min_id = &min(@blockId_array);
					$block{$min_id} = join ' ', @tmp_1;
					for my $blockId(@blockId_array){
						if($min_id != $blockId){
							my @tmp_blockId;
							if($block{$block_pos} =~ /\s+/){
								@tmp_blockId = split /\s+/, $block{$blockId};
							}else{
								@tmp_blockId = ($block{$blockId});
							}
							delete $block{$blockId};
						}
					}	
					last;
				}
			}
#################################
		}
	}
	my $block_count = scalar (keys %block);
	print "$bul\t$t\t$count_while\t$countipp\t$block_count\t$loop_countp\n";
	$bul = 1; 
}

######################################################################################
#########################     Cluster with the untagsnp point   ######################
######################################################################################

my %char; my $no = 0;
my %record_pos1;my @one;
my %record_pos2;

my @apos; my %untag;


for my $i(keys %tag){
	my $pos = (split /_/, $tag{$i}{"info"})[1] if exists $tag{$i}{"info"};
	if(not exists $itag{$i})  {$untag{$pos} = $i;}  
}	


@apos = ((keys %block), (keys %untag));
@apos = sort {$a<=>$b} @apos;
my %reco_bkItem; 

for my $i(keys %block){
    my @block_item;
    if($block{$i} =~ /\s+/){
        @block_item = split /\s+/, $block{$i};
    }else{
        @block_item = ($block{$i});
    }
    my @map = sort @block_item;
    my %rep_num2; my @rep_num2_1;
    for my $rp1(@map){
        $rep_num2{$rp1} = 0;
        for my $rp2(@map){
            if(exists $tag{$rp1}{"R_squre"}{$rp2} && $tag{$rp1}{"R_squre"}{$rp2}>=0.8){
                $rep_num2{$rp1}++;
            }
        }
    }
    my $rpt_i;

    for my $rp3(sort keys %rep_num2){
        if (not exists $itag{$rp3}){
            $rpt_i = 0;
        }else{
            $rpt_i = ($rep_num2{$rp3}+1)/@map;
        }
        if($rpt_i>0.85){
            push @rep_num2_1, "$rpt_i";
        }else{
            push @rep_num2_1, "0";
        }
    }
    my $max_id2 = &max_id(@rep_num2_1);
    my $reco_item = $block_item[$max_id2];
	$reco_bkItem{$i} = $reco_item;
}

my %block2; my $flag11;

for my $order(1..$#apos){
	if(not exists $untag{$apos[$order]}){next;}
	my $flag1 = 0; my $flag2 = 0;
	my $ro = $order; my $lo = $order;
	my @array_ro; my @array_lo;
	while(1){
		$ro = $ro - 1;
		$lo = $lo + 1;
		if($ro >= 0 && $flag1 == 0){		
			if(exists $block{$apos[$ro]}){push @array_ro, $ro;}
			if(@array_ro >= 4){$flag1 = 1;}
		}else{
			$flag1 = 1;
		}

		if($lo < @apos && $flag2 == 0){
			if(exists $block{$apos[$lo]}){push @array_lo, $lo;} 
			if(@array_lo >= 4){$flag2 = 1;}
		}else{
			$flag2 = 1;
		}   	
		last if ($flag1==1 && $flag2==1);
	}
	
 	my @bok_order = (@array_ro, @array_lo);
	my $id = $untag{$apos[$order]}; 	
	my %item_rq; my %rep_num; my %rep_num_n;
	for my $i(@bok_order){
		$rep_num_n{$i} = 0;
		my @block_item;
		if($block{$apos[$i]} =~ /\s+/){
			@block_item = split /\s+/, $block{$apos[$i]};
		}else{
			@block_item = ($block{$apos[$i]});
		}
		my @map = sort @block_item;
		for my $rp1(@map){
			if(exists $tag{$id}{"R_squre"}{$rp1} && $tag{$id}{"R_squre"}{$rp1} >= 0.8){$rep_num_n{$i}++;}
		}
		$rep_num{"$apos[$i]\t$rep_num_n{$i}"} = 1;
		
		my $reco_item = $reco_bkItem{$apos[$i]};
		if(exists $tag{$reco_item}{"R_squre"}{$id}){
			my $rq = $tag{$reco_item}{"R_squre"}{$id};
			$item_rq{"$apos[$i]\t$rq"} = 1;
		}else{
			$item_rq{"$apos[$i]\t0"} = 1;
		}
	}
	my @key_item_rq = keys %item_rq;
	my ($flag10, $max_pos) = &max_space(@key_item_rq);
	if($flag10 == 0){
		$block{$max_pos} = join ' ', ($block{$max_pos}, $id);
	}elsif($flag10 == 2){
		my @array_r_count;
        for my $i5(@bok_order){
            my $reco_item = $reco_bkItem{$apos[$i5]};
            my @block_item;
            if($block{$apos[$i5]} =~ /\s+/){
                @block_item = split /\s+/, $block{$apos[$i5]};
            }else{
                @block_item = ($block{$apos[$i5]});
            }
            my $r_count = 0;
            for my $i6(@block_item){
                if(exists $tag{$reco_item}{"R_squre"}{$i6} && $tag{$reco_item}{"R_squre"}{$i6} >= 0.8){   ## modify at 10-22, 23:48
                    $r_count++;
                }
            }
            my $dealer = $r_count/(@block_item+1);
            push @array_r_count, "$apos[$i5]\t$dealer";
        }
        ($flag11,$max_pos) = &max_space2(@array_r_count);
        if($flag11==0){
            $block{$max_pos} = join ' ',($block{$max_pos},$id);
        }elsif($flag11==1){
            $block2{$apos[$order]} = $id;
        }
	}else{
		my @key_rep_num = keys %rep_num;
		$max_pos = &max_space1(@key_rep_num);
		$block{$max_pos} = join ' ', ($block{$max_pos}, $id);
	}
}				

undef %char;

########################## out put #########################################

my $count_id = 0; my $iregion_start; my $iregion_end; 
my $iregion_length; my $region;
@block{keys %block2} = values %block2;

open OUT, ">$out" or die $!;

for my $i(sort {$a<=>$b} keys %block){
	$count_id++;
	print OUT "# $count_id %% $i\n";
	my $item; my $maf;
	my ($item_id, $item_pos);
	my @print_block;
	if($block{$i} =~ /\s+/){
		@print_block = split /\s+/, $block{$i};
	}else{
		@print_block = ($block{$i});
	} 
	for my $zh(@print_block){
##		if(not exists $tag{$zh}{"info"}){next;}
		$item_id .= "$zh; ";
		my $xh_pos =  (split /_/, $tag{$zh}{"info"})[1];
		$item_pos .= "$xh_pos; ";
		my $imaf = (split /_/, $tag{$zh}{'info'})[-1];	
		$imaf = &Sprintf($imaf);
		$maf .= "$imaf; ";
		my $iregion = $tag{$zh}{'region'};
		eval {$region .= "$iregion; ";};
	}
	print OUT "\tid: $item_id\n";
	print OUT "\tposition: $item_pos\n";
######## out put "MAF" and rep_rate and rep_region of each point in every block #####
	my $rpt; my $propt;
	my @map = @print_block;
	my %rep_num2; my @rep_num2_1;
	for my $rp1(@map){
##		if(not exists $tag{$rp1}{"R_squre"}){next;}
		$rep_num2{$rp1} = 0;
		for my $rp2(@map){
			if(exists $tag{$rp1}{"R_squre"}{$rp2} && $tag{$rp1}{"R_squre"}{$rp2}>=0.8){
				$rep_num2{$rp1}++;			
			}
		}
	}
	my $rpt_i;
	for my $rp3(keys %rep_num2){
		if (not exists $itag{$rp3}){
			$rpt_i = 0;
		}else{
			$rpt_i = ($rep_num2{$rp3}+1)/@map;
		}
		$rpt_i = &Sprintf($rpt_i);
		$rpt .= "$rpt_i; ";
		if($rpt_i>0.85 && exists $itag{$rp3}){
			$propt .= "1; ";
			push @rep_num2_1, "$rpt_i";
		}else{
			$propt .= "0; ";
			push @rep_num2_1, "0";
		}
	}
	my $max_id2 = &max_id(@rep_num2_1);
	my $item_count = @map;
	print OUT "\tmaf:$maf\n";
	print OUT "\trep_rate: $rpt\n";
	print OUT "\tregion: $region\n";
	print OUT "\tproperty: $propt\n";
	print OUT "\trecomend tagsnp: $map[$max_id2]\n";
	print OUT "\tnumber of points: $item_count\n";
	$region = '';
}	

####################  Subrouting programs  ######################### 

sub max{
    my @a = @_;
    my $max = (split /-/,$a[0])[0];
	my $tmp = (split /-/,$a[0])[1];
	my $imax = "$max\-$tmp";
    for my $i(@a){
		my @b = split /-/, $i;
        if($b[0]>$max){
            $max = $b[0];
			$imax = "$max\-$b[1]";
        }
    }
	return $imax;
}

sub max_space{
    my @a = @_;
    @a = sort {(split /\s+/, $a)[1]<=>(split /\s+/, $b)[1]} @a;
    my $max = (split /\s+/,$a[0])[1];
    my $imax = (split /\s+/,$a[0])[0];
    my $flag = 0;
    my $start = $max;
    my $end = (split /\s+/, $a[$#a])[1];
    if($start==$end && $end!=0){$flag = 1; return $flag, $imax;}
    if($end == 0){$flag = 2; return $flag, $imax;}
    $imax = (split /\s+/, $a[$#a])[0];
    return $flag, $imax;
}


sub max_space1{
    my @a = @_;
    my $max = (split /\s+/,$a[0])[1];
    my $imax = (split /\s+/,$a[0])[0];
    for my $i(@a){
        my @b = split /\s+/, $i;
        if($b[1]>$max){
            $max = $b[1];
            $imax = $b[0];
        }
    }
    return $imax;
}

sub max_space2{
    my @a = @_;
    @a = sort {(split /\s+/, $a)[1]<=>(split /\s+/, $b)[1]} @a;
    my $max = (split /\s+/,$a[0])[1];
    my $imax = (split /\s+/,$a[0])[0];
    my $flag = 0;
    my $end = (split /\s+/, $a[$#a])[1];
    my $start = $max;
    if($end<0.85){$flag = 1; return $flag, $imax;}
    $imax = (split /\s+/, $a[$#a])[0];
    return $flag, $imax;
}


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
                        return $middle, $middle;
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
        return $right, $left;
}

sub Sprintf
{
	my ($a) = @_;
	my $b = sprintf ("%.2f", $a);
    return $b;
}

sub max_id
{
	my @a = @_;
	my $max_id = 0;
	for my $i(0..$#a){
		if($a[$i] > $a[$max_id]){
			$max_id = $i;
		}
	}
	return $max_id;
} 

sub min
{
	my @a = @_;
	my $min = $a[0];
	for my $i(@a){
		if($i < $min){
			$min = $i;
		}
	}
	return $min;
}

sub maxi 
{
	my @a = @_;
    my $max = $a[0];
    for my $i(@a){
        if($i > $max){
            $max = $i;
        }
    }
    return $max;
}

sub index_binarySerach
{
	my ($hash1, $value) = @_; 
	my $middle = 0;
	my $pos = 0;
	my $right = 0;
	my $left = 0;	
	my $length = length $value;
	my $head_chr = "x" x (8-$length);
	my $value_1 = $head_chr.$value; 
	my $key1 = substr($value_1, 0, 3);	
	my $key2 = substr($value_1, 3, 3);
	my $key3 = substr($value_1, 6);
	my %hash = %$hash1;
	if(exists $hash{$key1}){
		if(exists $hash{$key1}{$key2}){
			if(exists $hash{$key1}{$key2}{$key3}){
				$pos = $hash{$key1}{$key2}{$key3};
				return $pos, $pos;
			}else{
				my $a1; my $b1;
				my @key_1_1 = @{$record_pos2{$key1}{$key2}};
				my @key_1 = map {$_=~s/x//g; if($_ eq ""){$_=0;} $_} @{$record_pos2{$key1}{$key2}}; 
				($left, $right) = &BinarySearch (\@key_1, $value);
				($left, $right) = ($hash{$key1}{$key2}{$key_1_1[$left]}, $hash{$key1}{$key2}{$key_1_1[$left]}+1);
			}
		}else{
			my $b1; my $a1;
			my @key_2_1 = @{$record_pos1{$key1}};
			my @key_2 = map {$_=~s/x//g; if($_ eq ""){$_=0;} $_} @{$record_pos1{$key1}};
			($left, $right) = &BinarySearch (\@key_2, $value);
			my @array2 = @{$record_pos2{$key1}{$key_2_1[$left]}};
			($left, $right) = ($hash{$key1}{$key_2_1[$left]}{$array2[$#array2]}, $hash{$key1}{$key_2_1[$left]}{$array2[$#array2]}+1); 	
		}												
	}else{
		my $a1; my $b1;
		my @key_3_1 = @one;
		my @key_3 = map {$_=~s/x//g; if($_ eq ""){$_=0;} $_} @one;
		($left, $right) = &BinarySearch (\@key_3, $value);
		my @array1 = @{$record_pos1{$one[$left]}};
		my @array2 = @{$record_pos2{$one[$left]}{$array1[$#array1]}};
		($left, $right) = ($hash{$key_3_1[$left]}{$array1[$#array1]}{$array2[$#array2]}, $hash{$key_3_1[$left]}{$array1[$#array1]}{$array2[$#array2]}+1); 
	}
	return $left, $right;	
}
