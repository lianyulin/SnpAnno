open IN,"$ARGV[0].summary.INDEL.txt" or die;
while(<IN>){
	chomp;
	next if (/^#Chr	Tran_ID/ || /Non-frameshift/);
	my @a=split/\t/;
	$hash{$a[1]}.="$_ TAG ";
	$count{$a[1]}++;
}
for my $k(sort keys %hash){
	next if ($count{$k} == 1);
	my @b=split/ TAG /,$hash{$k};
	if ($#b == 1){
		my @a1=split/\t/,$b[0];my @a2=split/\t/,$b[1];	
		($r1,$v1)=$a1[3]=~/(\S*)->(\S*)/;($r2,$v2)=$a2[3]=~/(\S*)->(\S*)/;
		$lenr1=length $r1;$lenr2=length $r2;
		$lenv1=length $v1;$lenv2=length $v2;
		if (($lenv1-$lenr1 + $lenv2-$lenr2) % 3 == 0 ){
			$list{"$a1[0]	$a1[1]	$a1[2]/$a2[2]	$a1[3]/$a2[3]	$a1[6]/$a2[6]	$a1[9]/$a2[9]"}=1;
		}
	}
	else{
		for my $item(@b){
			my @a=split/\t/,$item;
			($r,$v)=$a[3]=~/(\S*)->(\S*)/;
			$lenr=length $r;$lenv=length $v;
			$gap=$lenv-$lenr;
			if ($last_gap && ($gap+$last_gap) % 3 == 0){
				$list{"$a[0]	$a[1]	$last_p/$a[2]	$last_indel/$a[3]	$last_p_cds/$a[6]	$last_p_pro/$a[9]"}=1;
			}
			$last_item1=$item;$last_p=$a[2];$last_indel=$a[3];$last_p_cds=$a[6];$last_p_pro=$a[9];$last_gap=$gap;
		}
		undef $last_item1;undef $last_p;undef $last_indel;undef $last_p_cds;undef $last_p_pro;undef $last_gap;
	}
}
open OUT,">$ARGV[0].multiple-indel.txt";
print OUT "Chr	Transcript	Pos1/Pos2	Var1/Var2	Pos1_cds/Pos2_cds	Pos1_pro/Pos2_pro	Sig	Affected_cds	Affected_pro\n";
open co,"$ARGV[1]" or die "Condon file inexistence,please check your condon table";#condon table
while(<co>){if(/(\w+) (\S)/){$codon{$1}=$2;}}
open cds,"$ARGV[0].cds.mutation.fasta" or die;
while(<cds>){
	chomp;
	if (/^>(\S+)/){$id=$1;next;}
	$seq_cds{$id}.=$_;
}
for my $k(sort keys %list){
	my @a=split/\t/,$k;
	(my $p1,my $p2)=$a[4]=~/(\d*)\/(\d*)/;
	$temp_seq=substr($seq_cds{$a[1]},0,$p1);
	while($temp_seq=~/(\S{3,3})/g){if ($codon{$1} eq "*"){$tag_termination=1;}}
	if (!defined $tag_termination){
		if ($p1 % 3 == 1){$s=$p1-1;}
		if ($p1 % 3 == 2){$s=$p1-2;}
		if ($p1 % 3 == 0){$s=$p1-3;}
		if ($p2 % 3 == 1){$e=$p2+2;}
		if ($p2 % 3 == 2){$e=$p2+1;}
		if ($p2 % 3 == 0){$e=$p2;}
		$cds_indel=substr($seq_cds{$a[1]},$s,$e-$s);
		while($cds_indel=~/(\S{3,3})/g){
			$pro_indel.=$codon{$1};
		}
		if ($pro_indel=~/\*/){
			$tag_termination=2;
		}
		else{
			print OUT "$k	Y	$cds_indel	$pro_indel\n";
		}
	}
	if ($tag_termination){
		print OUT "$k	N	N.A.	N.A.\n";
	}
	undef $tag_termination;undef $cds_indel;undef $pro_indel;undef $temp_seq;
}