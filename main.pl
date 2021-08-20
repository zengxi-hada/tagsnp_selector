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
use File::Path qw(make_path remove_tree);

=head1 function
  
  This program is designed to run the pipeline of selecting tagsnps automatically.
  It will deploy several programs including 'maf_repeat_masker2.pl', 'BLAST', 'define_block2_mofify5.pl', 'rm_repeatForList_md.pl',
  'rm_repeat_md.pl', 'trans_ID.pl',  and output the defined blocks using to select tagsnps. 

=head1 usage
  
  perl $0
		-frq	<str>		the input frq file
		-ld	<str>		the input ld file
		-m	<int>		the max number of mutation in the process of blast [5]
		-r	<num>		the min value of R_squre [0.8]
		-s	<num>		length of each step in the process of defining block [0.01]
		-maf1	<num>		the min threhold of maf [0.02]
		-maf2	<num>		the max threhold of maf [0.45]
		-l	<int>		length of each probe [91]	 
		-f	<str>		reference file  
		-rp	<str>		the repeat masker file of whole genome [default database]
		-e	<num>		the parameter e of BLAST, Expectation value [2.96e-43]
		-rep	<num>		min value of representive rate in the process difine block [0.85]
		-od	<str>		output directory 
		-cut			cut frq file or not, it must used with "-size" [default not]
		-size	<int>		the size that the querry file is going to cut [200000]
		-qsub			auto qsub or not [you must set "mem" when using this switch] 
		-mem	<num>		the memory for submitting jobs onto sge system
		-h|-?			help
		  

=head1 version
  
  v1.0  2011-11-01

=head1 author
  
  zengxi@genomics.org.cn

=head1 example

  perl ./main.pl -frq ./chr1.frq  ./chr1.LD -f ./chr1.fa -od . -cut -size 200000 -qsub -mem 5g

=cut

my ($frq, $ld, $outdir, $outshell, $help, $qsub, $mem);
my $maf1 = 0.02; my $maf2 = 0.45;
my ($mrsqure, $step, $max_mismatch, $length, $e) = (0.8, 0.01, 5, 91, 2.96e-43);
my $ref = "/ifs1/ST_REHEAL/USER/PUBLIC_database/database/Homo_sapiens/HG19_noRandom/blast_index/human.fa";
my $repeat_region = "/ifs1/ST_REHEAL/USER/PUBLIC_database/database/Homo_sapiens/HG19_noRandom/repeat_region/rmsk.txt.gz";
my $rep_rate = 0.85; my $size = 200000; my $cut;

GetOptions(
	'm=i' => \$max_mismatch,
	'r=f' => \$mrsqure,
	's=f' => \$step,
	'od=s' => \$outdir,
	'frq=s' => \$frq,
	'maf1=f' => \$maf1,
	'maf2=f' => \$maf2,
	'e=f' => \$e,
	'rp=s' => \$repeat_region,
	'rep=f' => \$rep_rate,
	'ld=s' => \$ld,
	'f=s' => \$ref,
	'l=i' => \$length,
	'-size=i' => \$size,
	'cut' => \$cut,
	'qsub' => \$qsub,
	'mem=f' => \$mem,
	'h|?' => \$help,
);

die `pod2text $0` if(!$outdir || !$frq || !$ld || $help || !$ref);
my ($chr) = ($ref =~ /\/chr(\d+)/);

my $maf = "/ifs1/ST_REHEAL/USER/zengxi/bin/tag-snp/version_final/md3/maf_repeat_masker2.pl";
my $blast = "/ifs1/ST_REHEAL/USER/zengxi/bin/basic_soft/blast/blastall";
my $define_block = "/ifs1/ST_REHEAL/USER/zengxi/bin/tag-snp/version_final/md3/define_block2_mofify5.pl";
my $rm_repeat_forList = "/ifs1/ST_REHEAL/USER/zengxi/bin/tag-snp/version_final/md3/rm_repeatForList_md.pl";
my $rm_repeat = "/ifs1/ST_REHEAL/USER/zengxi/bin/tag-snp/version_final/md3/rm_repeat_md.pl";
my $trans_ID = "/ifs1/ST_REHEAL/USER/zengxi/bin/tag-snp/version_final/md3/trans_ID.pl";
my $ref_index = "/ifs1/ST_REHEAL/USER/PUBLIC_database/database/Homo_sapiens/HG19_noRandom/blast_index/human.fa"; 

make_path("$outdir") if (!-e $outdir);
$outshell = "$outdir/chr$chr.sh";

if(!-d $outdir){
	die "the setting output directory is not a directory! please modify your input\n";
}

chomp (my $lines = `wc -l $frq`);

$lines = (split /\s+/, $lines)[0];

if($cut){
	if($size > $lines){
		open OUT, ">$outshell" or die $!;
		print OUT "date > chr$chr.time\n";
		print OUT "date > maf_chr$chr.time\n";
		print OUT "perl $maf -frq $frq -maf1 $maf1 -maf2 $maf2 -l $length -f $ref -o $outdir/probe$chr.out -re $repeat_region -chr chr$chr\n";
		print OUT "date >> maf_chr$chr.time\n";
		print OUT "date > blast_chr$chr.time\n";
		print OUT "$blast -i $outdir/probe$chr.out -p blastn -m 9 -d $ref_index -o $outdir/blast$chr.out -e $e -b 2\n";
		print OUT "date >> blast_chr$chr.time\n";
		print OUT "perl $rm_repeat -b $outdir/blast$chr.out -o $outdir/tagsnp$chr.out -c chr$chr -m $max_mismatch\n";
		print OUT "perl $trans_ID $frq $outdir/tagsnp$chr.out $outdir/tagsnp$chr.final.out\n";
		print OUT "date > define_block_chr$chr.time\n";
		print OUT "perl $define_block -ld $ld -frq $frq -tagsnp $outdir/tagsnp$chr.final.out -minr $mrsqure -rep $rep_rate -stp $step -o $outdir/final_block$chr.file\n";
		print OUT "date >> define_block_chr$chr.time\n";
		print OUT "date >> chr$chr.time\n";
		my $jid_list = "";
		if($qsub){
			my $id = &autoQsub ($outshell, '2g', $jid_list);
			print STDERR "$id $outshell\n";
		}
	}else{
		chdir "$outdir";
		mkdir "chr$chr\_split" if !-e "chr$chr\_split";
		mkdir "chr$chr\_split/querry" if !-e "chr$chr\_split/querry";
		chdir "chr$chr\_split/querry"; 
		my $list = "$outdir/chr$chr\_split/chr$chr\_blast.list";
		open LIST, ">$list" or die $!;
		`split -$size $frq chr$chr\_frq`;
		chomp (my @split = `ls chr$chr\_frq*`); 
		chdir "$outdir";
		my $count = 0;
		make_path ("$outdir/chr$chr\_split/stp1") if (!-e "$outdir/chr$chr\_split/stp1");
		my $job_id; my @jid;
		for my $i(@split){
			$count++;
			open SH, ">$outdir/chr$chr\_split/stp1/chr$chr\_$count.sh" or die $!;
			my ($tail) = ($i=~/chr$chr\_frq(\w+)/);
#			print SH "date > chr$chr\_$count.time\n";
			print SH "date > maf_chr$chr\_$count.time\n";
			print SH "perl $maf -frq $outdir/chr$chr\_split/querry/$i -maf1 $maf1 -maf2 $maf2 -l $length -f $ref -o $outdir/chr$chr\_split/probe$chr\_$tail.out -re $repeat_region -chr chr$chr\n";
			print SH "date >> maf_chr$chr\_$count.time\n";
			print SH "date > blast_chr$chr\_$count.time\n";
			print SH "$blast -i $outdir/chr$chr\_split/probe$chr\_$tail.out -p blastn -m 9 -d $ref_index -o $outdir/chr$chr\_split/blast$chr\_$tail.out -e $e -b 2\n";		
			print SH "date >> blast_chr$chr\_$count.time\n";

			print LIST "$outdir/chr$chr\_split/blast$chr\_$tail.out\n";
			close SH;
			my $jid_list = "";
	        if($qsub){
    	        $job_id = &autoQsub ("$outdir/chr$chr\_split/stp1/chr$chr\_$count.sh", '2g', $jid_list);
				push @jid, $job_id;
	            print STDERR "$job_id $outdir/chr$chr\_split/stp1/chr$chr\_$count.sh\n";
	        }
		}
		
		make_path("$outdir/chr$chr\_split/stp2") if (!-e "$outdir/chr$chr\_split/stp2");
		open SH2, ">$outdir/chr$chr\_split/stp2/chr$chr.sh" or die $!;
		print SH2 "date > chr$chr.time\n";
		print SH2 "perl $rm_repeat_forList -b $outdir/chr$chr\_split/chr$chr\_blast.list -o $outdir/tagsnp$chr.out -c chr$chr -m $max_mismatch\n";
		print SH2 "perl $trans_ID $frq $outdir/tagsnp$chr.out $outdir/tagsnp$chr.final.out\n";
		print SH2 "date > define_block_chr$chr.time\n";
		print SH2 "perl $define_block -ld $ld -frq $frq -tagsnp $outdir/tagsnp$chr.final.out -minr $mrsqure -rep $rep_rate -stp $step -o $outdir/final_block$chr.file\n";
		print SH2 "date >> define_block_chr$chr.time\n";
		print SH2 "date >> chr$chr.time\n";
		close SH2;
		my $jid_list = join ",", @jid;
		if($qsub){
			my $id = &autoQsub ("$outdir/chr$chr\_split/stp2/chr$chr.sh", "$mem", $jid_list);
			print STDERR "$id $outdir/chr$chr\_split/stp2/chr$chr.sh\n";
		}
	}
}else{
	open OUT, ">$outshell" or die $!;
	print OUT "date > chr$chr.time\n";
    print OUT "date > maf_chr$chr.time\n";
	print OUT "perl $maf -frq $frq -maf1 $maf1 -maf2 $maf2 -l $length -f $ref -o $outdir/probe$chr.out -re $repeat_region -chr chr$chr\n";
	print OUT "date >> maf_chr$chr.time\n";
	print OUT "date > blast_chr$chr.time\n";
	print OUT "$blast -i $outdir/probe$chr.out -p blastn -m 9 -d $ref_index -o $outdir/blast$chr.out -e $e -b 2\n";	
	print OUT "date > blast_chr$chr.time\n";
	print OUT "perl $rm_repeat -b $outdir/blast$chr.out -o $outdir/tagsnp$chr.out -c chr$chr -m $max_mismatch\n";
	print OUT "perl $trans_ID $frq $outdir/tagsnp$chr.out $outdir/tagsnp$chr.final.out\n";
	print OUT "date > define_block_chr$chr.time\n";
	print OUT "perl $define_block -ld $ld -frq $frq -tagsnp $outdir/tagsnp$chr.final.out -minr $mrsqure -rep $rep_rate -stp $step -o $outdir/final_block$chr.file\n";
	print OUT "date >> define_block_chr$chr.time\n";
    print OUT "date >> chr$chr.time\n";
	my $jid_list = "";
	if($qsub){
		my $id = &autoQsub ($outshell, "$mem", $jid_list);
		print STDERR "$id $outshell\n";
	}
}

close OUT;

####################   subroutine   ############################

sub autoQsub {
    my ($script,$ram, $reportDir, $jid_list) = @_;
    my $hold = "";
    $hold = "-hold_jid $jid_list" if $jid_list;
    my $job_id;
    while (1) {
        my $info = `qsub -cwd $hold -l vf=$ram  -o $script`;
        $job_id = ($info =~ /Your job (\d+) /)[0];
        last if defined $job_id;
        sleep 180;
    }
    return $job_id;
}

sub showLog {
    my @t = localtime;
    printf STDERR "[%04d-%02d-%02d %02d:%02d:%02d] \t%s\n", $t[5] + 1900, $t[4] + 1, $t[3], $t[2], $t[1], $t[0], $_[0];
}

