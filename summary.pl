opendir DIR,"$ARGV[0]" or die "Can't find the directory $ARGV[0], check please\n";
for my $sample(readdir DIR){
	next if( "./$ARGV[0]/$sample" eq "./$ARGV[0]/." || "./$ARGV[0]/$sample" eq "./$ARGV[0]/.." ||  -f "./$ARGV[0]/$sample") ;
	opendir DIR1,"./$ARGV[0]/$sample";push (@sample,$sample);
	my @dir=grep/statistical_result/,readdir DIR1;
	for my $f(@dir){
		open result,"$ARGV[0]/$sample/$f" or die;
		while(<result>){
			chomp;
			if(/Total Number of (\S+)\s+: (\d+)/){$num{$1}{$sample}=$2;next;}
			if(/In Total                  : (\d+)/){$num{'Total'}{$sample}=$1;next;}
			if(/Region	/){$tag_region=1;next;}if(/####### Transition/){$tag_region=0;next;}
			if ($tag_region &&  /^\S+/){
				my @a=split/\t/;push (@region_list,$a[0]);
				$num_region{$a[0]}{$sample}{'SNV'}=$a[1];
				$num_region{$a[0]}{$sample}{'Insertion'}=$a[3];
				$num_region{$a[0]}{$sample}{'Deletion'}=$a[5];
			}
			if (/Region: (\S+)/){$temp_region=$1;}if (/Ts\/Tv : (\S+)/){$tstv{$temp_region}{$sample}=$1;}
		}
	}
}
print "Combine statistics result to $ARGV[1]\n";
use Excel::Writer::XLSX;
for("A".."ZZZZ"){$i++;$temp{$i}=$_;}
my $sta = Excel::Writer::XLSX->new("$ARGV[1]");
$format1 = $sta->add_format(
			bold	 => 0,
			align	 => 'centre',
			size	 => '14',
			font	 => 'Times New Roman',
			color	 => 'red'
		     );
$format2 = $sta->add_format(
				bold	 => 1,
				align	 => 'centre',
				size	 => '14',
				font	 => 'Times New Roman',
				 );
$format3 = $sta->add_format(
				align	 => 'centre',
				size	 => '12',
				font	 => 'Times New Roman',
				 );
$xlnum = $sta->add_worksheet("Mutation_Number");
$xlnum->merge_range("A1:C1","#Total Mutation Number",$format1);	

$xlnum->set_column($_,$_,15) for (0);$xlnum->set_column($_,$_,10) for (1..$#sample+2);
$xlnum->write("A2","Type",$format2);my $temp=2;
for my $t(0..$#sample){$t1=$t+2;$xlnum->write("$temp{$t1}2","$sample[$t]",$format2);}
for my $t('SNV','Insertion','Deletion','Total'){
	$temp++;$xlnum->write("A$temp","$t",$format2);
	for my $t1(0..$#sample){
		$num{$t}{$sample[$t1]} ? ($num{$t}{$sample[$t1]}) : ($num{$t}{$sample[$t1]}=0);
		$t2=$t1+2;$xlnum->write("$temp{$t2}$temp","$num{$t}{$sample[$t1]}",$format3);
	}
}
@region_list=grep{++$count{$_} < 2;} @region_list;$temp=9;

for my $r('SNV','Insertion','Deletion','tstv'){
	$temp=8+($n*15);$n++;$temp1=$temp+1;
	if ($r eq "tstv"){$xlnum->merge_range("A$temp:E$temp","#Ts\/Tv in Different Genomic Elements",$format1);}
	else{$xlnum->merge_range("A$temp:E$temp","#$r Number in Different Genomic Elements",$format1);}
	$xlnum->write("A$temp1","Region",$format2);
	for my $t(0..$#sample){$t1=$t+2;$xlnum->write("$temp{$t1}$temp1","$sample[$t]",$format2);}
	for my $t(@region_list){
		$temp1++;$xlnum->write("A$temp1","$t",$format2);
		for my $t1(0..$#sample){
			$num_region{$t}{$sample[$t1]}{$r} ? ($num_region{$t}{$sample[$t1]}{$r}) : ($num_region{$t}{$sample[$t1]}{$r}=0);
			$t2=$t1+2;$tstv{$t}{$sample[$t1]} ? ($tstv{$t}{$sample[$t1]}) : ($tstv{$t}{$sample[$t1]}="NA");
			if ($r eq "tstv"){$xlnum->write("$temp{$t2}$temp1","$tstv{$t}{$sample[$t1]}",$format3);}
			else{$xlnum->write("$temp{$t2}$temp1","$num_region{$t}{$sample[$t1]}{$r}",$format3);}
		}
	}
}