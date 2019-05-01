REF=/home/garner1/Work/genomes/Homo_sapiens.GRCh37.dna.primary_assembly.fa/GRCh37.fa
BED=/home/garner1/Work/genomes/Homo_sapiens.GRCh37.dna.primary_assembly.fa/GRCh37.bed
numbproc=24
datadir=/media/garner1/hdd/WES/$(exp)_$(run) 
target_exome=~/Work/dataset/agilent/S07604715_Covered.woChr.bed
target_xgen=~/Work/dataset/GenesSubset/xgen.bed
GTF=~/Work/genomes/Homo_sapiens.GRCh37.87.gtf

usage:
	@echo "Usage: make rawDataDir=fullpath baseCalls=fullpath out=fullpath demultiplex"
	@echo "       make R1=fullpath R2=fullpath R1.fastp.gz=fullpath R2.fastp.gz=fullpath sampleID=string preprocessing"
	@echo "       make exp=string run=string R1=fullpath R2=fullpath align"
	@echo "       make prefix=Prefix2bamfile dedup"
	@echo "       make chrom=chr prefix=fullpathPrefix dedup_parallel"
	@echo "       make bamfilePrefix=fullpathPrefix varcall_mpileup"
	@echo "       make bamfilePrefix=fullpathPrefix varcall_freebayes"
	@echo "       make bamfilePrefix=fullpathPrefix alfred"
	@echo "       make multisample_mpileup"
	@echo "       make inputVcf=fullpath2vcf.norm variant_stats"
	@echo "       make inputVcf=fullpath2vcfPrefix variant_filter"
	@echo "       make inputVcf=fullpath2vcf variant_effect_predictor"
	@echo "       make inputPrefix=fullPath2bam features_count"
	@echo "       make inputbam=fullpath2bam preseq"
	@echo "       make inputbam=fullpath2bam outdir=outputDir qualimap"
	@echo "       make bamfilePrefix=fullpathPrefix varscan"
	@echo "       make samplename=S01 bamfilePrefix=fullpathPrefix vardict"

vardict:
	vardict -G ${REF} -f 0.01 -N $(samplename) -b $(bamfilePrefix).dedup.sorted.bam -z -c 1 -S 2 -E 3 -g 4 ${target_exome} | teststrandbias.R | var2vcf_valid.pl -N $(samplename) -E -f 0.01 > $(bamfilePrefix).vardict.vcf

varscan:
	echo "[`date`]: Start running varscan ... "
	samtools mpileup -f ${REF} $(bamfilePrefix).dedup.sorted.bam > $(bamfilePrefix).mpileup
	varscan mpileup2snp $(bamfilePrefix).mpileup > $(bamfilePrefix).varscan.snp
	varscan mpileup2indel $(bamfilePrefix).mpileup > $(bamfilePrefix).varscan.indel
	varscan filter $(bamfilePrefix).varscan.snp --indel-file $(bamfilePrefix).varscan.indel --output-file $(bamfilePrefix).varscan.snp.filter
	varscan filter $(bamfilePrefix).varscan.indel --output-file $(bamfilePrefix).varscan.indel.filter
	varscan readcounts $(bamfilePrefix).mpileup --output-file $(bamfilePrefix).mpileup.readcounts
	echo "[`date`]: Done with varscan ... "

demultiplex:
	mkdir -p $(out)
	bcl2fastq -i $(baseCalls) -R $(rawDataDir) -o $(out) --no-lane-splitting --tiles s_[3]

preprocessing:
	fastp -i $(R1) -I $(R2) -o $(R1.fastp.gz) -O $(R2.fastp.gz) -h $(sampleID).report.html -j $(sampleID).report.json

align:
	mkdir -p ${datadir}
	bwa mem -v 1 -t ${numbproc} ${REF} $(R1) $(R2) | samtools sort -T $(exp) > ${datadir}/$(exp).bam

dedup:	
	/usr/local/share/anaconda2/bin/umi_tools dedup -I $(prefix).bam  --paired -S $(prefix).dedup.bam --umi-separator=: \
	--edit-distance-threshold 2 -L $(prefix).log --output-stats=$(prefix)
	samtools sort -@ 8 $(prefix).dedup.bam > $(prefix).dedup.sorted.bam

dedup_parallel:				#parallelize wrt chromosomes...much faster
	/usr/local/share/anaconda2/bin/umi_tools dedup --chrom $(chrom) -I $(prefix).bam --paired -S $(prefix).$(chrom).dedup.bam \
	--umi-separator=: --edit-distance-threshold 2 -L $(prefix).$(chrom).log --output-stats=$(prefix).$(chrom)
	# samtools merge -@ 24 $(prefix).dedup.bam $(prefix).$(chrom).dedup.bam
	# samtools sort -m 10G -@ 8 $(prefix).dedup.bam > $(prefix).dedup.sorted.bam

varcall_mpileup:
	samtools mpileup -uvf ${REF} $(bamfilePrefix).dedup.sorted.bam | bcftools call -vm -Ov > $(bamfilePrefix).mpileup.vcf
	bcftools norm -f ${REF} -o $(bamfilePrefix).mpileup.norm.vcf $(bamfilePrefix).mpileup.vcf

varcall_freebayes:
	freebayes -f ${REF} $(bamfilePrefix).dedup.sorted.bam > $(bamfilePrefix).freebayes.vcf
	bcftools norm -f ${REF} -o $(bamfilePrefix).freebayes.norm.vcf $(bamfilePrefix).freebayes.vcf

alfred:
	alfred qc -r ${REF} -b ${target_xgen} -j $(bamfilePrefix).json.gz -o $(bamfilePrefix).tsv.gz $(bamfilePrefix).dedup.sorted.bam

multisample_mpileup:
	cat ${BED} | parallel "bcftools mpileup -Ou -f /home/garner1/Work/genomes/Homo_sapiens.GRCh37.dna.primary_assembly.fa/GRCh37.fa -r {} /media/garner1/hdd/WES/S??_BC145/S??.dedup.sorted.bam | bcftools call -vm -Ov -o /media/garner1/hdd/WES/{}.mpileup.vcf.gz"

	cat ${BED} | parallel "bcftools norm -f /home/garner1/Work/genomes/Homo_sapiens.GRCh37.dna.primary_assembly.fa/GRCh37.fa -o /media/garner1/hdd/WES/{}.mpileup.norm.vcf /media/garner1/hdd/WES/{}.mpileup.vcf.gz"

variant_stats:
	bcftools stats -F ${REF} -s - $(inputVcf) > $(inputVcf).norm.vcf.stats

variant_filter:
	bcftools filter -o $(inputVcf).filtered.vcf -s LOWQUAL -e '%QUAL<20 || DP>250' $(inputVcf)

variant_effect_predictor:
	snpEff -Xmx30G -csvStats $(inputVcf).snpeff.csv GRCh37.75 $(inputVcf) > $(inputVcf).annotated.vcf

features_count:
	featureCounts -p -t exon -g gene_id -a ${GTF} -o $(inputPrefix).counts.txt $(inputPrefix).dedup.sorted.bam

preseq:
	preseq c_curve -v -pe -bam $(inputbam) -o $(inputbam).preseq.txt -l 3000000000

qualimap:
	qualimap bamqc -bam $(inputbam) --java-mem-size=3G -gff ${GTF} -outdir $(outdir)
