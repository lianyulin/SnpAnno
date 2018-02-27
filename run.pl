use Getopt::Long;
my %opts;
my $usage = <<"Usage";
  Usage:	perl run.pl -g <genome> -t <annotation> Options
        -g	string	genome file in fasta format, required
        -t	string	genomic annotation in gtf or gff format, required
  Options:
        -h	Usage
        -dir	string	path of variant files, put all files into one fold, default is directory named as "vcf" in the current path
        -o	string	path of output, default is directory named as "Output" in the current path
        -dp	int	filter reads depth equal or less than <int>, default 0
        -script	string	path of script, default is directory named as "script" in the current path
        -flank	int	length of gene flanking region, default is 3000 bp
        -indel		Report Indel annotation, default is no
        -cds		Report post-substituted coding sequence, default is no
        -protein	Report post-substituted protein sequence, default is no
        -mi		Report multiple INDELs prediction, default is no (don't select unless -cds and -indel are applied)
  PAV Options:
        -PAV		Presence/Absence Variation analysis for all variant sites, default is no
        -GT		Report Genotype, default is no (GT tag is required, don't select unless -PAV is applied) 	
        -R		Report allelic depth and ratio, default is no (AD and DP tag are required, don't select unless -PAV is applied)
        -region	string	Only report PAV results that belong to specific regions, splitby comma, optional \"intragenic\",\"intergenic\",\"TE\",\"upstream\",\"downstream\",\"cds\",\"all\", default is intragenic (don't select unless -PAV is applied)	
        -SP	string	Only report PAV results that belong to specific samples, split by comma, etc sample1,sample2, default is all samples (don't select unless -PAV is applied)
        -SPfile	string	File contains specific sample, each line contains one sample name (don't select unless -PAV is applied)
Usage
GetOptions (\%opts,"g=s","t=s","dir=s","o=s","dp=i","script=s","flank=i","h","indel","cds","protein","mi","PAV","GT","R","region=s","SP=s","SPfile=s") or die $usage;
#Check the existence of necessary file
if($opts{h}){die $usage;}
if(!$opts{g} || !(-e $opts{g})){die "$usage\nCan't not find genome file, no such file or directory\n";}
if(!$opts{t} || !(-e $opts{t})){die "$usage\nCan't not find annotation file, no such file or directory\n";}
if(!$opts{dir}){$opts{dir}="\.\/vcf";}
if(!(-e $opts{dir})){die "$usage\nCan't find directory: $opts{dir}, check please!\n";}
if(!$opts{o}){$opts{o}="\.\/Output";}
system "mkdir $opts{o}" if (!(-e $opts{o}));
if (!(-e $opts{o})){die "Can't find directory: $opts{o}, check please!\n";}
if(!$opts{script}){$opts{script}="\.\/script";}
if(!(-e $opts{script})){die "$usage\nCan't find directory: $opts{script}, check please!\n";}
if(!$opts{flank}){$opts{flank}=3000;}
if($opts{mi} && (!$opts{cds} || !$opts{indel})){die "$usage\nOptions '-mi' is executable only if '-cds' and '-indel' are selected\n";}
if(!$opts{PAV}){
	if ($opts{GT}){die "$usage\nOptions 'GT' is executable only if '-PAV' is selected\n";}
	if ($opts{R}){die "$usage\nOptions 'R' is executable only if '-PAV' is selected\n";}
	if ($opts{region}){die "$usage\nOptions 'region' is executable only if '-PAV' is selected\n";}
	if ($opts{SP}){die "$usage\nOptions 'SP' is executable only if '-PAV' is selected\n";}
	if ($opts{SPfile}){die "$usage\nOptions 'SPfile' is executable only if '-PAV' is selected\n";}
}
($spe_name)=$opts{g}=~/(\S+)\.fa\w*$/;
if (!$spe_name){die "Can't determine the name of genome, may be you need to modify the genome file as \"XXX.fa\(fasta\)";}

print "Check the files you input:
     Genome             :	$opts{g}
     Genomic annotation :	$opts{t}
     Variation dir      :	$opts{dir}
     Script dir         :	$opts{script}
     Output dir         :	$opts{o} \n";
#get sequence from genomic annotation
system "perl $opts{script}/get_transcript.pl $opts{g} $opts{t} $opts{o}/$spe_name $opts{script}/Codon.txt" ;
system "perl $opts{script}/get_utr_splice_site.pl $opts{g} $opts{t} $opts{o}/$spe_name";
system "samtools faidx $opts{g}";

#decode variation
opendir vcf,"$opts{dir}" or die "Can't find directory: $opts{dir}, please check!!";
@vcf=grep/.vcf$/,readdir vcf;
if (!@vcf){die "Can't find any vcf file in the directory: $opts{dir}, please check if file suffix is \".vcf\"";}
foreach(@vcf){
	($sample)=$_=~/(\S+)\.vcf$/;
	$cmd="bash $opts{script}/run.sh -i $opts{dir}/$_ -t $opts{t} -g $opts{g} -d $opts{script} -o $opts{o} -f $opts{flank}";
	if ($opts{indel}){$cmd.=" -I";}if ($opts{cds}){$cmd.=" -C";}if ($opts{protein}){$cmd.=" -P";}if ($opts{dp}){$cmd.=" -m $opts{dp}";}
	system $cmd;undef $cmd;
	if ($opts{mi}){system "perl $opts{script}/Multiple_indels.pl $opts{o}/$sample/$sample $opts{script}/Codon.txt";}
}
#summarize statistics
system "perl $opts{script}/summary.pl $opts{o} $opts{o}/statistics_result.xlsx";
#PAV
if ($opts{PAV}){
	$cmd="perl $opts{script}/pav.pl $opts{o} -o $opts{o}/PAV.result.txt";
	if ($opts{GT}){$cmd.=" -G";}if ($opts{R}){$cmd.=" -R";}
	if ($opts{region}){$cmd.=" -region $opts{region}";}if ($opts{SP}){$cmd.=" -sample $opts{SP}";}if ($opts{SPfile}){$cmd.=" -files $opts{SPfile}";}
	system $cmd;undef $cmd;
}