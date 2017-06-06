##-----------------------------------------------##
##                                               ##
## EVOP 2017 Project :                           ##
## Bad Bromance (SNPs caller testing)            ##
##                                               ##
##-----------------------------------------------##

import glob, os

##-----------------------------------------------##
## Set the number of samples                     ##
##-----------------------------------------------##

#In case of PE reads, you should set only R1 or R2
def _populations():
	Items = list()
	for foo,bar in config["POPULATION"].items():
		if foo in list(config["INDIVIDUALS"]):
			Items.append(str(bar))
	return Items

def _individuals():
	Items = list()
	for foo,bar in config["POPULATION"].items():
		for bas in config["INDIVIDUALS"][foo]:
			Items.append(str(bar)+"/"+str(bas))
	return Items

#def _individuals():
#	Items = list()
#	for foo,bar in config["POPULATION"].items():
#		for bas in config["INDIVIDUALS"][foo]:
#			Items.append(str(bas))
#	return Items

def _mapping_params(individual):
	group=os.path.dirname(individual)
	return config["MAPPER_PARAM"][group]

##-----------------------------------------------##
## Define the Reference location and parameters  ##
##-----------------------------------------------##

REFERENCE=config["REFERENCE"] #REF = "/home/Evop/Z_test_pipe/Data/Ref_ChrEven.fa"
GENOME=config["GENOME"]
NAME=config["NAME"]
WHERE=config["WHERE"]
OTHERSCRIPT=config["OTHERSCRIPT"] # I don't need this 
POPULATION=_populations()
INDIVIDUALS=_individuals() #SAMPLES, = glob_wildcards("/home/Evop/Z_test_pipe/Data/{sample}.R1.fastq.gz") 

##-----------------------------------------------##
## Set software directories                      ##
##-----------------------------------------------##
MAPPER = config["MAPPER"] #"/usr/local/bin/bwa-0.7.12/bwa"
SAMTOOLS = config["SAMTOOLS"] #"/usr/local/bin/samtools-1.3.1/samtools"
JAVA_1_8 = config["JAVA_1_8"] #"/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.121-0.b13.el7_3.x86_64/jre/bin/java"
PICARD = config["PICARD"] #"/usr/local/bin/picard-1.140/dist/picard.jar"
#GATK = config["GATK"] #"/usr/local/bin/GenomeAnalysisTK.jar"
#FREE = config["FREEBAYES"] #"/usr/local/bin/freebayes/bin/freebayes"
#BCFTOOLS = config["BCFTOOLS"] #"/usr/local/bin/bcftools"
#TABIX = config["TABIX"] #"/usr/local/bin/htslib-1.2.1/tabix"
#PLOTS = config["PLOTVCF"] #"/usr/local/bin/plot-vcfstats"
#VCFFILTER = config["VCFFILTER"] #"/usr/local/bin/vcffilter"
#VCFTOOLS = config["VCFTOOLS"] #"/usr/local/bin/vcftools"

# The rg LB PL and other SM... Have to be changed accordingly to your samples.
#RG = "'@RG\\tID:{sample}\\tLB:{sample}\\tPL:ILLUMINA\\tSM:{sample}'" ##There are problems putting this here! (Sincerely I don't know why! I changed to the rule paragraph and it works!)

##------------------------------##
## Define all the ruled         ##
##------------------------------##

rule all:
	input: expand("{where}/Mapped/{individual}.realn.bam", where=WHERE, population=POPULATION, individual=INDIVIDUALS),
		   expand("{where}/Mapped/Statistics/{individual}.idx.txt", where=WHERE, population=POPULATION, individual=INDIVIDUALS),
		   expand("{ref}{gen}/{gen}{name}.bwt",ref=REFERENCE,gen=GENOME, name=NAME),
		   expand("{ref}{gen}/{gen}{name}.fai",ref=REFERENCE,gen=GENOME, name=NAME),
		   expand("{where}/Mapped/{individual}.sort.bam.bai", where=WHERE, population=POPULATION, individual=INDIVIDUALS),
		   expand("{where}/Mapped/Statistics/{individual}.flag.txt", where=WHERE, population=POPULATION, individual=INDIVIDUALS),
		   expand("{where}/Mapped/{individual}.rgmdup.bam.bai",where=WHERE, population=POPULATION, individual=INDIVIDUALS),
		   expand("{where}/Mapped/{individual}.rgmdup.bam.bai",where=WHERE, population=POPULATION, individual=INDIVIDUALS),
		   expand("{where}/Mapped/{individual}.realn.bam.bai",where=WHERE, population=POPULATION, individual=INDIVIDUALS)

rule MapperIndex: 
	input:  expand("{ref}{gen}/{gen}{name}.fa",ref=REFERENCE,gen=GENOME, name=NAME)
	output: expand("{ref}{gen}/{gen}{name}.bwt",ref=REFERENCE,gen=GENOME, name=NAME) 
	message: "Creating the Bwa mem index for {REF}"
	shell: "{MAPPER} index {input}"

rule Alignment:
	input:
	   fwd = "{where}/Raw/{population}/{individual}.R1.fastq",
	   rev = "{where}/Raw/{population}/{individual}.R2.fastq",
	   ref = expand("{ref}{gen}/{gen}{name}.fa",ref=REFERENCE,gen=GENOME, name=NAME)
	output:  "{where}/Mapped/{population}/{individual}.bam"
	message: " Dear user, we are pleased to tell you that you will be spending a sweet time in our company while your aligmnent of {input} is going on... Enjoy and relax... And watch your back!"
	params:  RG = "'@RG\\tID:{individual}\\tLB:{individual}\\tPL:ILLUMINA\\tSM:{individual}'" # I remeber we had some issue here, but can't remember them!!
	log:	"LOGS/mapping_{individual}.log"
	threads: 20
	shell:   "{MAPPER} -M -R {params.RG} -t {threads}  {ref} {input.fwd} {input.rev} | {SAMTOOLS} view -bS - -o {output}"

rule SortingBAM:
	input:
		"{where}/Mapped/{population}/{individual}.bam"
	output:
		"{where}/Mapped/{population}/{individual}.sort.bam"
	threads: 1
	message: "Sorting the {input}"
	shell: "{SAMTOOLS} sort {input} -o {output}"

rule indexing_sort_bam:
	input: 
		"{where}/Mapped/{population}/{individual}.sort.bam"
	output:
		"{where}/Mapped/{population}/{individual}.sort.bam.bai"
	message: " Sorted bam files get indexed for {input}"
	shell:	"{SAMTOOLS} index {input} {output}"

rule idx_statistics_generation:
	input: 
		"{where}/Mapped/{population}/{individual}.sort.bam"
	output:
		"{where}/Mapped/Statistics/{population}/{individual}.idx.txt"
	message: "generating idx statistics for rg bam files"
	shell: "{SAMTOOLS} idxstats {input}"

rule flag_statistics_generation:
	input: 
		"{where}/Mapped/{population}/{individual}.sort.bam"
	output:
		"{where}/Mapped/Statistics/{population}/{individual}.flag.txt"
	message: "generating flag statistics for rg bam files"
	shell: 
		"{SAMTOOLS} flagstat {input}"

rule Mark_duplication:
	input:
		"{where}/Mapped/{population}/{individual}.sort.bam"
	output:
		out = "{where}/Mapped/{population}/{individual}.rgmdup.bam" ,
		metrics = "{where}/Mapped/Statistics/{population}/{individual}.rgmdup.metrics"
	message:"generating metrics statistics"
	shell:
		"{JAVA_1_8} -jar {PICARD} MarkDuplicates I={input} O={output.output} METRICS_FILE={output.metrics} MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000"

rule rmdupBam_index_generation:
	input: 
		"{where}/Mapped/{population}/{individual}.rgmdup.bam"
	output:
		"{where}/Mapped/{population}/{individual}.rgmdup.bam.bai"
	message:"generating bamIndex for .rgmdup.bam files"
	shell:
		"{SAMTOOLS} index {input} {output}"

rule SamtoolsIndex:
	input: 
		expand("{ref}{gen}/{gen}{name}.fa",ref=REFERENCE,gen=GENOME, name=NAME)
	output: 
		expand("{ref}{gen}/{gen}{name}.fai",ref=REFERENCE,gen=GENOME, name=NAME)
	message: 
		"Creating the Bwa mem index for {input}"
	shell: "{SAMTOOLS} faidx {input}"


rule intervals_generation:
	input:
		"{where}/Mapped/{population}/{individual}.rgmdup.bam"
	output:
		"{where}/Mapped/{population}/{individual}.rgmdup.intervals"
	message:
		"generating intervals file for rg bam files FooDuFafa ::Nice Elevator Jazzy song::"
	shell:
		"{JAVA_1_8} -Xmx40g -d64 -jar {GATK} -T RealignerTargetCreator -nt 12 -R {REF} -o {output} -I {input}" 

## No clue why we need it?
rule realigning_rg_bam_files:
	input:
		input = "{where}/Mapped/{population}/{individual}.rgmdup.bam" ,
		intervals = "{where}/Mapped/{population}/{individual}.rgmdup.intervals"
	output:
		"{where}/Mapped/{population}/{individual}.realn.bam"
	message:"Realignment of {input.input}"
	shell:
		"{JAVA_1_8} -Xmx40g -d64 -jar {GATK} -T IndelRealigner -R {REF} -targetIntervals {input.intervals} -o {output} -I {input.input} -compress 0 --maxReadsInMemory 100000"

rule realn_index_generation:
	input: 
		"{where}/Mapped/{population}/{individual}.realn.bam"
	output:
		"{where}/Mapped/{population}/{individual}.realn.bam.bai"
	message: "generating bam index for .realn.bam files"
	shell:
		"{SAMTOOLS} index {input} {output}"

## No clue why we need it?
##rule Fixmate_BAM:
##    input:
##        "0.MapData/{sample}.sort.bam"
##    output:
##        "0.MapData/{sample}.fix.bam"
##    threads: 1
##    message: "And now we keep doing fixmate to the {input}"
##    shell:
##        "{SAMTOOLS} fixmate -O bam {input} -o {output}"
#
##rule indexing_fixed_bam:
##    input: 
##        "0.MapData/{sample}.fix.bam"
##    output:
##        "0.MapData/{sample}.fix.bam.bai"
##    message: " Sorted bam files get indexed for {input}"
##    shell:
##        "{SAMTOOLS} index {input} {output}"
#
## No clue why we need it?
##rule PicardDict:
##	input:
##        "{REF}"
##    output:
##        "{REF}.dict"
##    message: "Creating the Bwa mem index for {REF}"
##    shell:
##        "{JAVA_1_8} -jar {PICARD} CreateSequenceDictionary R={input} O={output} "
#
##rule BaseRecalibrator_generation:
##    input:
##        "0.MapData/{sample}.realn.bam"
##    output:
##        "0.MapData/{sample}.recal.data.table"    
##    message:"Creating a Base Recalibration"
##    shell:
##        "{JAVA_1_8} -Xmx20g -d64 -jar {GATK} -T BaseRecalibrator -nct 12 -I {input} -R {REF} -o {output}"
##
##rule PrintReads_generation:
##    input:
##        input = "0.MapData/{sample}.realn.bam" ,
##        intervals = "0.MapData/{sample}.recal.data.table" 
##    output:
##        "0.MapData/{sample}.recalibrated.bam"    
##    message:"generating printreads file for .realn.bam files"
##    shell:
##        "{JAVA_1_8} -Xmx10g -d64 -jar {GATK} -T PrintReads -nct 12 -I {input.input} -R {REF} -BQSR {input.intervals} -o {output}"
##
##rule indexing_sorted_bam3:
##    input: 
##        "0.MapData/{sample}.recalibrated.bam"
##    output:
##        "0.MapData/{sample}.recalibrated.bam.bai"
##    message: "generating bam index"
##    shell:
##        "{SAMTOOLS} index {input} {output}"
##
##End of common pipeline
#
#
#
###-----------------------------------------------##
###-----------------------------------------------##
#
#
###GATK
#
##rule GATK_STEP_HaplotypeCaller:
##    input:
##        expand("0.MapData/{sample}.recalibrated.bam" ,sample=SAMPLES)
##    output:
##        "GVCF/{sample}.vcf"
##    message:"calling variants for realigned bam files with the GATK pipeline SSSSSsssssSSSSSssSsSSSSsnaake it Baby Bromance"
##    shell:
##        "{JAVA_1_8} -Xmx5g -d64 -jar {GATK} -T HaplotypeCaller -R {REF} -I  --max_alternate_alleles 2 --emitRefConfidence GVCF -variant_index_type LINEAR -variant_index_parameter 128000 -o {output}"
##
### Changes compared with the original haplotype step ## Steph comments :
### stand_emit_conf 20 ## Be careful they mention some stand_emit_conf and stand_call_conf as default for the best practices. 
### stand_call_conf 20 
### mbq 10 
##
##rule joint_Genotyping: ##This command needs to be done to all samples at the same time#
##    input:
##        "GVCF/{sample}.vcf"
##    output:
##        "GVCF/All_GATK.vcf"
##    message:""
##    shell:
##        "ls -1 {input} | {JAVA_1_8} -cp -jar {GATK} -T GenotypeGVCFs -R {REF} --variant  -o {output}"
##
#### PLUS : Here I added the "ls -1 {input} | " step to create a list of samples it should work !
##
##rule GATK_HardFiltering_VCF:
##    input:
##        "GVCF/All_GATK.vcf"
##    output:
##        "GVCF/All_HF_GATK.vcf"
##    message:"Hard filtering the SNPs"
##    shell:
##        "{JAVA_1_8} -jar {GATK} -T VariantFiltration -R {REF} -V {input} --filterExpression QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 --filterName my_snp_filter -o {output}"
##
##rule PASS_SNPs_Filtering_VCF:
##    input:
##        "GVCF/All_HF_GATK.vcf"
##    output:
##        "GVCF/All_HF_SNPsPASS_GATK.recode.vcf"
##    message:"calling variants for realigned bam files with the GATK pipeline SSSSSsssssSSSSSssSsSSSSsnaake it Baby Bromance"
##    shell:
##        "{vcftools} --vcf {input} --remove-indels --remove-filtered-all --recode --recode-INFO-all --out {output}"
##
##rule GATK_recalibration_VCF:
##    input:
##        vcf = "GVCF/All_GATK.vcf" ,
##	SNPdb = "GVCF/All_HF_SNPsPASS_GATK.recode.vcf"	
##    output:
##        recal = "GVCF/recalibrate_SNP.recal" ,
##	tranches = "GVCF/recalibrate_SNP.tranches" ,
##	rscript = "GVCF/recalibrate_SNP_plots.R"
##    message:"calling variants for realigned bam files with the GATK pipeline SSSSSsssssSSSSSssSsSSSSsnaake it Baby Bromance"
##    shell:
##        "{JAVA_1_8} -Xmx5g -d64 -jar {GATK} -T VariantRecalibrator -R {REF} -input {input.vcf} - resource:HARD,known=false,training=true,truth=true,prior=10.0 {input.SNPdb} -an DP -an QD -an FS -an MQRankSum -an ReadPosRankSum -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile {output.recal} -tranchesFile {output.tranches} -rscriptFile {output.rscript}"
##
#### NEED HELP WHY SO MANY BRANCHES ?
##
##rule GATK_apply_recalibration_VCF:
##    input:
##        vcf = "GVCF/All_GATK.vcf" ,
##        recal = "GVCF/recalibrate_SNP.recal" ,
##	tranches = "GVCF/recalibrate_SNP.tranches" ,
##	rscript = "GVCF/recalibrate_SNP_plots.R"
##    output:
##        "GVCF/SNPcall_GATK_recalibrated.vcf"
##    message:"calling variants for realigned bam files with the GATK pipeline SSSSSsssssSSSSSssSsSSSSsnaake it Baby Bromance"
##    shell:
##        "{JAVA_1_8} -jar {GATK} -T ApplyRecalibration -R {REF} -input {input.vcf} -mode SNP --ts_filter_level 90 -recalFile {input.recal} -tranchesFile {input.tranches} -rscriptFile {input.rscript} -o {output}"
##
##
##rule PASS_SNPs_Final_Filtering_VCF:
##    input:
##        "GVCF/SNPcall_GATK_recalibrated.vcf"
##    output:
##        "GVCF/SNPcall_GATK_filtered.recode.vcf"
##    message:"calling variants for realigned bam files with the GATK pipeline SSSSSsssssSSSSSssSsSSSSsnaake it Baby Bromance"
##    shell:
##        "{vcftools} --vcf {input} --remove-indels --remove-filtered-all --recode --recode-INFO-all --out {output}"
##
####Samtools
##
##rule samtools_call:
##    input:
##        bam = expand("0.MapData/{sample}.recalibrated.bam", sample=SAMPLES) ,
##        bai = expand("0.MapData/{sample}.recalibrated.bam.bai", sample=SAMPLES)
##    output:
##        "GVCF/SNPcall_Samtools.vcf.gz"
##    shell:
##        "{SAMTOOLS} mpileup -ugf {REF} {input.bam} | {BCFTOOLS} call -vmO z -o {output}"
##
##rule tabix_prepare:
##    input:
##        "GVCF/SNPcall_Samtools.vcf.gz"
##    output:
##        "GVCF/SNPcall_Samtools.vcf.gz.idx"
##    shell:
##        "{tabix} -p vcf {input}"
##
#### HERE why not just add the pipe to the previous line ?
##
##rule bcftools_stats:
##    input:
##        vcf = "GVCF/SNPcall_Samtools.vcf.gz" ,
##	ref = "{REF}"
##    output:
##        "GVCF/SNPcall_Samtools.vcf.gz.stats"
##    shell:
##        "{BCFTOOLS} stats -F {input.ref} -s - {input.vcf} > {output}"
##
##rule bcftools_plots:
##    input:
##         vcf = "GVCF/SNPcall_Samtools.vcf.gz.stats"
##    output:
##        "plots/plots" #I don't know what will be the output
##    shell:
##        "{plot} -p plots/ {input} > {output}"
##
##rule bcftools_filter_ST:
##    input:
##        "GVCF/SNPcall_Samtools.vcf.gz"
##    output:
##        "GVCF/SNPcall_Samtools_filtered.vcf.gz"
##    shell:
##        "{BCFTOOLS} filter -O z -o {output} -s LOWQUAL -i '%QUAL>10' {input}"
##
##
####FreeBayes
##
##rule FREE_BAYES:
##    input:
##	  expand("0.MapData/{sample}.recalibrated.bam", sample=SAMPLES)
##        #"bam_file_list.txt"   ## Should be a list with the bam localization and names! I don't know which is the best option, could be create a rule to list all the recalibrated files!     
##    output:
##        "GVCF/SNPcall_FreeBayes.vcf"
##    message:"calling variants Free Bayes Bad Bromance "
##    shell:
##        "ls -1 {input} | {FREE} -L  -f {REF} -v {output}"
##
##rule bcftools_filter_FB:
##    input:
##        "GVCF/SNPcall_FreeBayes.vcf"
##    output:
##        "GVCF/SNPcall_FreeBayes_filtered.vcf"
##    shell:
##        "{vcffilter} -f \"QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1\" {input} > {output} "
##
#### Quality Metrics on variants
### Tiv Ratio (2.1 for WGS ~2.8 for exome)
### HapMap concordance
### SNV/Indel Counts
### Rare variant enrichment
### DP
### Q
### GQ 
#
#
#
