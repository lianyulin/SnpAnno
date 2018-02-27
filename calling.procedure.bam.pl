foreach(@ARGV){
	if(/-h/){$help=1;}
	if(/(\S+)=(\S+)/){
		$$1=$2;
	}
	else{
		$genome=$_;
	}
}
if(!$genome || !($genome=~/\S+.fa/) || $help){
	die "
Usage: perl calling.procedure.bam.pl <genome.fasta> Option
	Option:    
	cpu=<INT>, default cpu=1 
	data_type=<string>, <DNA> or <RNA>, default data_type=DNA;
	calling_type=<string>, <HaplotypeCaller> or <UnifiedGenotyper>, default calling_type=HaplotypeCaller;
	FS=<Float>, default FS=30.0;
	MQ=<Float>, default MQ=40.0;
	QD=<Float>, default QD=2.0;
	cluster=<INT>, default cluster=3;
	windows=<INT>, default windows=35;
	-h, print the Usage
	Can't not find your genome file or it's not fa/fasta format!Please check!\n";
}
if(!$genome || !($genome=~/\S+.fa/) || $help){print ""}
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
=cut
`samtools faidx $genome`;
`java -jar picard.jar CreateSequenceDictionary R= $genome O= $index.dict`;
`bowtie2-build $genome $index`;
=cut
opendir DIR,"./bam/" or die "Can't not find your bam fold, please set up a fold named bam,then put your bam file into this fold!";
@BAM=grep/bam$/,readdir DIR;
foreach(@BAM)
{
	if(/(\S+).bam$/){$isolate=$1;}
	system("java -Djava.io.tmpdir=tmp -jar picard.jar MarkDuplicates I=./bam/$isolate.bam O=$isolate\_duplication.bam M=$isolate\_metrics VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES=4000\n");		
	system("samtools index $isolate\_duplication.bam");
	if($data_type eq "RNA" || $data_type eq "rna" ){
		system("java -Djava.io.tmpdir=tmp -jar GenomeAnalysisTK.jar -T SplitNCigarReads -R $genome -I $isolate\_duplication.bam -o $isolate\_snc.bam -fixNDN -U ALLOW_N_CIGAR_READS -rf ReassignOneMappingQuality\n");
		$bam="$isolate\_snc.bam";
		$fixnd=" -fixNDN";
	}
	if($data_type eq "DNA" || $data_type eq "dna"){
		$bam="$isolate\_duplication.bam";$fixnd="";
	}
	
	system("java -Djava.io.tmpdir=tmp -jar GenomeAnalysisTK.jar -nt $cpu -T RealignerTargetCreator -R $genome -I $bam -o $isolate.intervals\n");
	system("java -Djava.io.tmpdir=tmp -jar GenomeAnalysisTK.jar -T IndelRealigner -R $genome -targetIntervals $isolate.intervals -I $bam -o $isolate\_realign.bam\n");
	if($calling_type eq "HaplotypeCaller"){
		system("java -Djava.io.tmpdir=tmp -jar GenomeAnalysisTK.jar -nct $cpu -T HaplotypeCaller -stand_call_conf 20 -stand_emit_conf 20 -ERC GVCF -variant_index_type LINEAR -variant_index_parameter 128000 -o $isolate.raw.vcf -R $genome -I $isolate\_realign.bam\n");
		system("java -jar GenomeAnalysisTK.jar -nt $cpu -T GenotypeGVCFs -R $genome -o $isolate.genotype.vcf -V $isolate.raw.vcf\n");
	}
	if($calling_type eq "UnifiedGenotyper"){
		system("java -Djava.io.tmpdir=tmp -jar -nt $cpu GenomeAnalysisTK.jar -T UnifiedGenotyper -R $genome -I $isolate\_realign.bam -o $isolate.genotype.vcf -glm BOTH --output_mode EMIT_VARIANTS_ONLY -allowPotentiallyMisencodedQuals");
	}
	system("java -jar GenomeAnalysisTK.jar -T VariantFiltration -R $genome -V $isolate.genotype.vcf -window $windows -cluster $cluster -filterName FS -filter \"FS > $FS\" -filterName QD -filter \"QD < $QD\" -filterName MQ -filter \"MQ < $MQ\" -o $isolate.tag.vcf\n");
	system("perl ./script/filter.tag.pl -tag FS:QD");
}
`rm -rf tmp`;
