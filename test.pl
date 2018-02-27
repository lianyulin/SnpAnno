=cut
use Bio::SeqIO;
my $in  = Bio::SeqIO->new(-file => "TAIR10_chr_all.fasta" ,-format => 'fasta');  #

while(my $obj=$in->next_seq()){

     my $id=$obj->id;                        #   序列的id

     my $desc=$obj->desc;                    #   序列的描述

     my $seq=$obj->seq;                      #   序列字符串

     my $seqlength=$obj->length;            #   序列的长度

     print  "$id\t$seqlength\n";                     #   以一列id一列序列的方式输出

}

open IN,"Can-o_phaseI.1.intragenic.vcf" or die;

open OUT,">temp1.txt";
while(<IN>){
	chomp;
	next if (/^#/);
	my @a=split/\t/;
	if ($a[-1] =~ /splice5:/){
		next if ($mask{$a[0]."\t".$a[1]});
		print OUT "$_\n";
		$mask{$a[0]."\t".$a[1]}=1;
	}
	
}

my $a="1	7440148	.	A	C,T	734.77	SnpCluster	AC=1,1;AF=0.500,0.500;AN=2;DP=31;Dels=0.00;ExcessHet=3.0103;FS=0.000;HaplotypeScore=0.0000;MLEAC=1,1;MLEAF=0.500,0.500;MQ=36.99;MQ0=0;QD=23.70;SOR=0.756	GT:AD:DP:GQ:PL	1/2:0,21,10:31:99:763,545,518,218,0,158	AT1G21250.1 Intragenic:7439266-7442113|exon:7439266-7440355|cds:7439512-7440355";
my @b=split/\s+/,$a;
(my $info1, my $info)=$a=~/	(\S*[AGD][DTP]:\S*)	(\S+)/;
		my @tag=split/:/,$info1;my @tag_value=split/:/,$info;
		for my $t(0..$#tag){
			if($tag[$t] eq "GT"){$tag_num{'GT'}=$t;}
			if($tag[$t] eq "AD"){$tag_num{'AD'}=$t;}
			if($tag[$t] eq "DP"){$tag_num{'DP'}=$t;}
		}
		@temp=split/,/,$tag_value[$tag_num{'AD'}];
				$temp[1] >= $temp[2] ? ($supported_ratio=$temp[1]/$tag_value[$tag_num{'DP'}]):($supported_ratio=$temp[2]/$tag_value[$tag_num{'DP'}]);
				$temp[1] >= $temp[2] ? $b[4]=~s/,\w*// : $b[4]=~s/\w*,//;
				print "$type1	$type2	$temp[1]	$temp[2]	$b[4]	$supported_ratio\n";



open IN,"Can-o_phaseI.cds.positioning.txt";
while(<IN>){
	if($_!~/^(\S+)	(\S+)(\s*.*)	(\d+)	\((\d+)\)	(\d+)	(.*)([+-])	(\S+)	(\S+)	type=(\S+)/){
		print "$_\n";
	}
}

$a="1	AT1G01070.2	40583	14";
 ($s)=$a=~/(\S+)$/;#$s++;$a=~s/\S+$/$s/;
print "$s\n";

open genome,"Can-o_phaseI.protein.mutation.fasta" or die "Genome file inexistence,please check your Genome file";
while(<genome>){
	chomp;
	if(/>(\S+)/){$id=$1;next;}
	$seq{$id}.=$_;
}
open list,"list.txt";
open OUT,">temp1.txt";
while(<list>){
	chomp;
	print OUT ">$_\n$seq{$_}\n";
}

$a="Can-o.phaseI.intragenic.vcf";
($b)=$a=~/\S+.(\S+).vcf/;
print "$b\n";

$a=5;$b=4;$d=44;$e=55;
$e > $d ? ($c=$b,$c1=$a):($c=$a,$c1=$b);
print "$c	$c1\n";

$a="1	6063	.	T	C	1276.77	PASS	AC=2;AF=1.00;AN=2;DP=42;Dels=0.00;ExcessHet=3.0103;FS=0.000;HaplotypeScore=0.8321;MLEAC=2;MLEAF=1.00;MQ=59.23;MQ0=0;QD=30.40;SOR=0.892	GT:AD:DP:GQ:PL	1/1:0,42:42:99:1305,120,0	AT1G01020.1 Intragenic:5928-8737|exon:5928-6263|UTR3:5928-6263";
@a=split/\t/,$a;
my @tag=split/:/,$a[-3];my @tag_value=split/:/,$a[-2];
			for my $t(0..$#tag){
				if($tag[$t] eq "GT"){$temp_tag.="$tag_value[$t];";}
				if($tag[$t] eq "DP"){$temp_dp=$tag_value[$t];}
				if($tag[$t] eq "AD"){($read_r,$read_v)=$tag_value[$t]=~/(\d*),(\d*)/;}
			}
			if (defined $read_v){$temp_dp ? ($temp_tag.=(sprintf "%.2f",$read_v/$temp_dp * 100 )."%") : ($temp_tag.=(sprintf "%.2f",$read_v/($read_r+$read_v) * 100)."%");}
(my $temp1,my $temp2,my $s,my $e)=$a[-1]=~/^(\S+) .+\|(\S+):(\d*)-(\d*)$/;
print "$temp1	$temp2	$temp_tag\n";

opendir DIR,"$ARGV[0]" or die "There is not the directory\n";
@dir=`find $ARGV[0] -type f`;
print @dir;

for("A".."ZZZZ"){$i++;$hash{$i}=$_;}

for (sort {$a <=> $b} keys %hash){
	print "$_	$hash{$_}\n";
}
=cut
$b="1	CDS	3631	3913	+	AT1G01010.1";
my @a=split/\t/,$b;$region=$a[1];
$$region{$a[-1]}{$a[2]}=$a[2]."-".$a[3];
print "$CDS{'AT1G01010.1'}{'3631'}";