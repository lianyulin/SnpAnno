use Getopt::Long;
my %opts;
my $usage = <<"Usage";
	Usage:	perl pav.pl annotation_dir <options>
        annotation_dir	Directory of annotated result
	Options:
        -h		Usage
        -G		Report Genotype (GT tag is necessary), default is no
        -R		Report allelic depth and ratio (AD and DP tag are necessary), default is no
        -region	string	Only report PAV results belong to specific regions, split by comma, optional \"intragenic\",\"intergenic\",\"TE\",\"upstream\",\"downstream\",\"cds\",\"all\", default is intragenic
        -sample	string	Only report PAV results belong to specific samples, split by comma, etc sample1,sample2, default is all samples
        -files	string	File contains specific samples, each line contains one sample name
        -o	string	output
Usage
GetOptions (\%opts,"h","G","R","region=s","sample=s","files=s","o=s") or die $usage;
if($opts{h}){die $usage;}
for ("intragenic","intergenic","TE","upstream","downstream","cds"){$tag_region{$_}=1;}
if(!$opts{region}){push (@vcf,"intragenic");}
elsif($opts{region}=~/all/i){push (@vcf,"intragenic","intergenic","TE","upstream","downstream","cds");}
else{
	@vcf=split/,/,$opts{region};
	for my $temp(@vcf){
		if (!$tag_region{$temp}){die "Can't not recognize region: $temp, only \"intragenic\",\"intergenic\",\"TE\",\"upstream\",\"downstream\",\"cds\",\"all\" are optional\n";}
	}
}
if ($opts{sample} && $opts{files}){die "$usage\n-sample and -files can't be selected at the same time\n";}
if ($opts{files}){
	open filesample,"$opts{files}" or die "Can't not file the file $opts{files}, check please\n";
	while(<filesample>){chomp;push (@dir,$_);}
}
elsif ($opts{sample}){@dir=split/,/,$opts{sample};}
else{
	opendir DIR,"$ARGV[0]" or die "Can't not find the directory of $ARGV[0], check please\n";
	@dir =grep/\S+/,readdir DIR;
}
if (!$opts{o}){$opts{o}="$ARGV[0]/PAV.result.txt";}

foreach $sample(@dir){
	next if( "./$ARGV[0]/$sample" eq "./$ARGV[0]/." || "./$ARGV[0]/$sample" eq "./$ARGV[0]/.." ||  -f "./$ARGV[0]/$sample") ;
	if (!(-d "./$ARGV[0]/$sample")){print "Can't not find the directory: $sample, please check it's included in $ARGV[0]\n";next;}
	print "Loading the sample information : $sample\n";push (@spe,$sample);
	for my $f(@vcf){
		if ($f=~/cds/i){$tag_cds=1;next;}
		open vcf,"./$ARGV[0]/$sample/$sample.$f.vcf" or die "Can't find the file ./$ARGV[0]/$sample/$sample.$f.vcf \n";
		while(<vcf>){
			chomp;
			my @a=split/\t/;next if (/^#/);
			next if ($a[-1] =~ /cds/ );
			$len_r=length $a[3];$len_v=length $a[4];
			if($len_r == $len_v && $len_r == 1){$tp='SNV';}
			if($len_r < $len_v ){$tp='Insertion';}
			if($len_r > $len_v ){$tp='Deletion';}
			if($opts{G} || $opts{R}){
				my @tag=split/:/,$a[-3];my @tag_value=split/:/,$a[-2];
				for my $t(0..$#tag){
					if($tag[$t] eq "GT" && $opts{G}){$temp_tag.="\|$tag_value[$t]";}
					if($tag[$t] eq "DP" && $opts{R}){$temp_dp=$tag_value[$t];}
					if($tag[$t] eq "AD" && $opts{R}){($read_r,$read_v)=$tag_value[$t]=~/(\d*),(\d*)/;}
				}
				if (defined $read_v){$temp_dp ? ($temp_tag.="\|$read_v,".(sprintf "%.2f",$read_v/$temp_dp*100 )."%") : ($temp_tag.="\|$read_v,".(sprintf "%.2f",$read_v/($read_r+$read_v)*100)."%");}
			}
			if ($f eq "intergenic"){$temp1="-";$temp2="intergenic";$s="";$e="";}
			elsif($f eq "intragenic"){($temp1,$temp2,$s,$e)=$a[-1]=~/^(\S+) .+\|(\w+):(\d*)-(\d*)$/;}
			else{($temp1,$temp2,$s,$e)=$a[-1]=~/^(\S+) (\S+):(\d*)-(\d*)/;}
			$total{"$a[0]	$a[1]	$temp2	$temp1	$s-$e"}{$sample}="$a[3]->$a[4]\|$tp$temp_tag";
			$sort1_id{"$a[0]	$a[1]	$temp2	$temp1	$s-$e"}=$a[0];$sort2_id{"$a[0]	$a[1]	$temp2	$temp1	$s-$e"}=$a[1];
			undef $temp_tag;undef $temp_dp;undef $read_r;undef $read_v;
		}
	}
	if (!$opts{region} || $opts{region}=~/cds|intragenic/ || $tag_cds ){
		open SNV,"./$ARGV[0]/$sample/$sample.summary.SNV.txt" or die;
		while(<SNV>){
			chomp;next if (/^#/);my @a=split/\t/;
			$cds_value{$a[0]."\t".$a[1]."\t".$a[2]}="$a[3],$a[7]->$a[8],$a[10]->$a[11],$a[12]\|$a[4]";
		}
		open INDEL,"./$ARGV[0]/$sample/$sample.summary.INDEL.txt";
		while(<INDEL>){
			chomp;next if (/^#/);my @a=split/\t/;
			$cds_value{$a[0]."\t".$a[1]."\t".$a[2]}="$a[3],$a[7]->$a[8],$a[10]->$a[11],$a[12]\|$a[4]" if ($a[12] =~ /^Non-frameshift/);
			$cds_value{$a[0]."\t".$a[1]."\t".$a[2]}="$a[3],$a[12]\|$a[4]" if ($a[12] =~/^Frameshift/);
		}
		open cds,"./$ARGV[0]/$sample/$sample.cds.group.vcf" or die;
		while(<cds>){
			chomp;
			my @a=split/\t/;next if (/^#/);
			my @tag=split/:/,$a[-3];my @tag_value=split/:/,$a[-2];
			if($opts{G} || $opts{R}){
				my @tag=split/:/,$a[-3];my @tag_value=split/:/,$a[-2];
				for my $t(0..$#tag){
					if($tag[$t] eq "GT" && $opts{G}){$temp_tag.="\|$tag_value[$t]";}
					if($tag[$t] eq "DP" && $opts{R}){$temp_dp=$tag_value[$t];}
					if($tag[$t] eq "AD" && $opts{R}){($read_r,$read_v)=$tag_value[$t]=~/(\d*),(\d*)/;}
				}
				if (defined $read_v){$temp_dp ? ($temp_tag.="\|$read_v,".(sprintf "%.2f",$read_v/$temp_dp*100 )."%") : ($temp_tag.="\|$read_v,".(sprintf "%.2f",$read_v/($read_r+$read_v)*100)."%");}
			}			
			(my $temp1,my $s,my $e)=$a[-1]=~/^(\S+) .+\|\w+:(\d*)-(\d*)$/;
			if($cds_value{"$a[0]	$temp1	$a[1]"}){
				$total{"$a[0]	$a[1]	CDS	$temp1	$s-$e"}{$sample}=$cds_value{"$a[0]	$temp1	$a[1]"}.$temp_tag;
				$sort1_id{"$a[0]	$a[1]	CDS	$temp1	$s-$e"}=$a[0];$sort2_id{"$a[0]	$a[1]	CDS	$temp1	$s-$e"}=$a[1];
			}
			else{
				$temp_po=$a[1]+1;
				if($cds_value{"$a[0]	$temp1	$temp_po"} && $cds_value{"$a[0]	$temp1	$temp_po"}=~/frameshift/i){
					$total{"$a[0]	$a[1]	CDS	$temp1	$s-$e"}{$sample}=$cds_value{"$a[0]	$temp1	$temp_po"}.$temp_tag;
					$sort1_id{"$a[0]	$a[1]	CDS	$temp1	$s-$e"}=$a[0];$sort2_id{"$a[0]	$a[1]	CDS	$temp1	$s-$e"}=$a[1];
				}
			}
			undef $temp_tag;undef $temp_dp;undef $read_r;undef $read_v;
		}
		undef %cds_value;
	}
}

open OUT,">$opts{o}" or die "Can not creat file $opts{o}\n";
print OUT "Chr	Position	Region	Transcript	Start-End	Ratio";
foreach (sort @spe){print OUT "	$_";}print OUT "\n";
foreach my $key1(sort {$sort1_id{$a} cmp $sort1_id{$b} || $sort2_id{$a} <=> $sort2_id{$b} }keys %sort1_id){
	foreach my $key2(sort @spe){
		$total{$key1}{$key2} ? ($cmd.="	$total{$key1}{$key2}",$exist_sample++) : ($cmd.="	N.A.");
	}
	$percentage_exist=(sprintf "%.2f",$exist_sample/($#spe+1) * 100 )."%";
	print OUT "$key1	$percentage_exist$cmd\n";
	undef $cmd;undef $exist_sample;
}