foreach(@ARGV){
	if(/-h/){$help=1;}
	if(/(\S+)=(\S+)/){
		$$1=$2;
	}
}
$genome=$ARGV[0];
if($ARGV[1]=~/gtf|gff/ && !($ARGV[1]=~/=/)){$gtf=$ARGV[1];}
if(!$genome || !($genome=~/\S+.fa/) || $help){
	die "
Usage: perl calling.procedure.reads.pl <genome.fasta> <annotation.gtf> Options
	genome.fasta, string, necessary, genome sequence
	annotation, string, annotation file, in the format of gtf or gff, it's necessary if data type is RNA
	Options:    
	cpu=<INT>, default cpu=1 
	data_type=<string>, <DNA> or <RNA>, default data_type=DNA;If data type is RNA,you must supply gtf or gff file,or there is not need to supply
	calling_type=<string>, <HaplotypeCaller> or <UnifiedGenotyper>, default calling_type=HaplotypeCaller;
	FS=<Float>, default FS=30.0;
	MQ=<Float>, default MQ=40.0;
	QD=<Float>, default QD=2.0;
	cluster=<INT>, default cluster=3;
	windows=<INT>, default windows=35;
	-h, print the Usage\n"
}
if(!$genome){print "Can't not find your genome file $genome or it's not fa/fasta format!Please check!";}
if(!$cpu){$cpu=1;}
if(!$data_type){$data_type="DNA";}
if(!$calling_type){$calling_type="HaplotypeCaller";}
if(!$FS){$FS="30.0";}
if(!$MQ){$MQ="40.0";}
if(!$QD){$QD="2.0";}
if(!$windows){$windows=35;}
if(!$cluster){$cluster=3;}
($index)=$genome=~/(\S+)\.fa/;
`mkdir vcf`;

print " 
	Genome file: $genome 
	Annotation file : $gtf
	Index: $index 
	CPU uasge: $cpu
	Input: $data_type
	Calling  type: $calling_type
	FS_value :$FS  QD_value : $QD  MQ_value : $MQ
	slide windows : $windows   cluster_number : $cluster
";


unless($data_type eq "DNA" or $data_type eq "RNA")
{die "Please check data type, please enter RNA or DNA!!";}
unless($calling_type eq "HaplotypeCaller" or $calling_type eq "UnifiedGenotyper")
{die "Please check calling type, please enter HaplotypeCaller or UnifiedGenotyper!!";}

`samtools faidx $genome`;
`java -jar picard.jar CreateSequenceDictionary R= $genome O= $index.dict`;
`bowtie2-build $genome $index`;

opendir DIR,"./reads/" or die "Can't not find your reads fold, please set up a fold named reads,then put your reads file into this fold!";
@fq=grep/fastq|fq/,readdir DIR;
foreach(@fq)
{
	if(/(.*)_paired.1\.(f.+)/)
	{	
		$isolate=$1;
		$left=$1."_paired.1.".$2;
		$right=$1."_paired.2.".$2;
		die "Be sure that your reads is pair end, we can't find your right reads, check please!!" if (!$right);
		if($data_type eq "DNA"){system("bowtie2 --rg-id $1 --rg \"PL:ILLUMINA\" --rg \"SM:$1\" -x $index -1 ./reads/$left -2 ./reads/$right -p $cpu -S $1.sam");}
		if($data_type eq "RNA"){`mkdir tophat`;system("tophat -p $cpu -g 1 -r 200 -i 20 -I 2200 --no-coverage-search --mate-std-dev 80 --keep-fasta-order -o tophat -G $gtf --rg-sample $1 --rg-id $1 --rg-library $1 --rg-platform illumina $index ./reads/$left ./reads/$right");`mv tophat/accepted_hits.bam ./$isolate.sort.bam`;}
	}
	if(/(.*)_singled\.(f.+)/){
		$isolate=$1;
		if($data_type eq "DNA"){system("bowtie2 --rg-id $1 --rg \"PL:ILLUMINA\" --rg \"SM:$1\" -x $index -U ./reads/$1\_singled.$2 -p $cpu -S $1.sam");}
		if($data_type eq "RNA"){`mkdir tophat`;system("tophat -p $cpu -g 1 -r 200 -i 20 -I 2200 --no-coverage-search --mate-std-dev 80 --keep-fasta-order -o tophat -G $gtf --rg-sample $1 --rg-id $1 --rg-library $1 --rg-platform illumina $index ./reads/$1\_singled.$2");`mv tophat/accepted_hits.bam ./$isolate.sort.bam`;}
	}
	if($data_type eq "DNA"){
		system("samtools view -bS $isolate.sam > $isolate.bam");		
		system("samtools sort $isolate.bam $isolate.sort");		
	}
	system("java -Djava.io.tmpdir=tmp -jar picard.jar MarkDuplicates I=$isolate.sort.bam O=$isolate\_duplication.bam M=$isolate\_metrics VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES=4000\n");		
	system("samtools index $isolate\_duplication.bam\n");		
	if($data_type eq "RNA" || $data_type eq "rna" ){
		system("java -Djava.io.tmpdir=tmp -jar GenomeAnalysisTK.jar -T SplitNCigarReads -R $genome -I $isolate\_duplication.bam -o $isolate\_snc.bam -fixNDN -U ALLOW_N_CIGAR_READS -rf ReassignOneMappingQuality\n");
		$bam="$isolate\_snc.bam";
	}
	if($data_type eq "DNA" || $data_type eq "dna"){
		$bam="$isolate\_duplication.bam";
	}
	system("java -Djava.io.tmpdir=tmp -jar GenomeAnalysisTK.jar -fixNDN -T RealignerTargetCreator -R $genome -I $bam -o $isolate.intervals\n");
	system("java -Djava.io.tmpdir=tmp -jar GenomeAnalysisTK.jar -T IndelRealigner -fixNDN -R $genome -targetIntervals $isolate.intervals -I $bam -o $isolate\_realign.bam\n");
	if($calling_type eq "HaplotypeCaller"){
		system("java -Djava.io.tmpdir=tmp -jar GenomeAnalysisTK.jar -T HaplotypeCaller -fixNDN -filterRNC -filterMBQ -filterNoBases -dontUseSoftClippedBases -rf UnmappedRead -rf BadMate -rf DuplicateRead -rf NotPrimaryAlignment -rf MappingQualityUnavailable -stand_call_conf 20 -stand_emit_conf 20 -ERC GVCF -variant_index_type LINEAR -variant_index_parameter 128000 -o $isolate.raw.vcf -R $genome -I $isolate\_realign.bam\n");
		system("java -jar GenomeAnalysisTK.jar -T GenotypeGVCFs -fixNDN -R $genome -o $isolate.genotype.vcf -V $isolate.raw.vcf\n");
	}
	if($calling_type eq "UnifiedGenotyper"){
		system("java -Djava.io.tmpdir=tmp -jar GenomeAnalysisTK.jar -T UnifiedGenotyper -R $genome -I $isolate\_realign.bam -o $isolate.genotype.vcf -glm BOTH --output_mode EMIT_VARIANTS_ONLY -allowPotentiallyMisencodedQuals");
	}
	system("java -jar GenomeAnalysisTK.jar -T VariantFiltration -R $genome -V $isolate.genotype.vcf -window $windows -cluster $cluster -filterName FS -filter \"FS > $FS\" -filterName QD -filter \"QD < $QD\" -filterName MQ -filter \"MQ < $MQ\" -o $isolate.tag.vcf\n");
	system("perl ./script/filter.tag.pl");
	
}
`rm -rf tmp`;
