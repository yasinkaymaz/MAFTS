
conda activate MedGen
#Fastq QC

fq1='dna1_S1_L001_R1_001.fastq.gz'
fq2='dna1_S1_L001_R2_001.fastq.gz'
sample="dna1"
#Adapter trimming
fastp -i $fq1 -I $fq2 -o "$sample".R1.fq.gz -O "$sample".R2.fq.gz


bowtie2 -x ~/Documents/data/Hg38/GRCh38_idx -1 "$sample".R1.fq.gz -2 "$sample".R2.fq.gz -S "${sample}".sam

samtools view -b "${sample}".sam > "${sample}".bam

samtools sort "${sample}".bam -o "${sample}".sorted.bam

docker run -v ~/Documents:/gatk/my_data -it broadinstitute/gatk:4.1.3.0

ln -s /gatk/my_data/data/Hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.* ./

#samtools faidx hg38_v0_GRCh38.primary_assembly.genome.fa

#gatk --java-options "-Xmx50G" CreateSequenceDictionary -R hg38_v0_GRCh38.primary_assembly.genome.fa

gatk AddOrReplaceReadGroups \
    -I "${sample}".sorted.bam \
    -O "${sample}".sorted_RG.bam \
    -ID 1 \
    -LB "lib1" \
    -PL "ILLUMINA" \
    -PU "unit1" \
    -SM "dna1"

gatk BuildBamIndex \
    -I "${sample}".sorted_RG.bam

gatk Mutect2 \
   -R resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
   -I "${sample}".sorted_RG.bam \
   --f1r2-tar-gz f1r2.tar.gz \
   -O "${sample}".vcf.gz

gatk LearnReadOrientationModel \
    -I f1r2.tar.gz \
    -O read-orientation-model.tar.gz

#Downloaded resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf.gz from GATK bundle
gatk IndexFeatureFile -F resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf.gz

gatk GetPileupSummaries \
    -I "${sample}".sorted_RG.bam \
    -V /gatk/my_data/data/Hg38/hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz \
    -L /gatk/my_data/data/Hg38/hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz \
    -O getpileupsummaries.table

gatk CalculateContamination \
    -I getpileupsummaries.table \
    -tumor-segmentation segments.table \
    -O calculatecontamination.table

gatk FilterMutectCalls \
    -R resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
    -V "${sample}".vcf.gz \
    --min-allele-fraction 0.1 \
    --min-reads-per-strand 5 \
    --unique-alt-read-count 5 \
    --tumor-segmentation segments.table \
    --contamination-table contamination.table \
    --ob-priors read-orientation-model.tar.gz \
    -O "${sample}".filtered.vcf.gz





conda activate DAB

InBam=$1
InBam=DNA001r1.sorted.bam
bamCoverage -b $InBam \
    -o output.sorted.bw \
    --binSize 20 \
    --normalizeUsing CPM \
    --smoothLength 60 \
    --extendReads 150 \
    --centerReads \
    -p 6 2> bamCoverage.log

computeMatrix scale-regions \
    -R target_regions_chr.bed \
    -S output.sorted.bw \
    --skipZeros \
    -p 6 \
    --regionBodyLength 500 \
    -a 500 -b 500 \
    -o "${InBam%.bam}"_sites.gz

plotProfile -m "${InBam%.bam}"_sites.gz \
    -out "${InBam%.bam}"_DeepProfile.png \
    --perGroup  --plotTitle "" \
    --samplesLabel "Targeted Regions" \
    -T "Read alignment sites"  -z "" \
    --startLabel "" \
    --endLabel "" \
    --colors red


snakemake -s ~/Dropbox/codes/MAFTS/Snakefile \
    --configfile ~/Dropbox/codes/MAFTS/config-sub.yml \
    --cores 4 \
    --latency-wait 20 \
    --use-conda \
    --dry-run

snakemake --forceall --rulegraph -s ~/Dropbox/codes/MAFTS/Snakefile --configfile ~/Dropbox/codes/MAFTS/config-sub.yml |dot -Tpdf > dag.pdf


for d in `ls -1|cut -d_ -f1|uniq`;
do
    echo $d;
    ls -al "${d}"_R1.fastq.gz;

    for f in `seq -f "%f" 0.2 0.2 1.0|cut -b 1-3`;
    do
        echo $f;
        seqtk sample "${d}"_R1.fastq.gz $f |gzip > "${d}"_"${f}"_R1.fastq.gz
        seqtk sample "${d}"_R2.fastq.gz $f |gzip > "${d}"_"${f}"_R2.fastq.gz
    done
done


seqtk sample DNA001r2_R1.fastq.gz 0.2 |gzip > DNA001r2_0.2_R1.fastq.gz
seqtk sample DNA001r2_R2.fastq.gz 0.2 |gzip > DNA001r2_0.2_R2.fastq.gz



for i in `ls -1|cut -d_ -f1|uniq|sed 's/dna//g'`;
do
    echo $i;
    mv dna"${i}"_S*_L001_R2_001.fastq.gz DNA0"${i}"r2_R2.fastq.gz;
    mv dna"${i}"_S*_L001_R1_001.fastq.gz DNA0"${i}"r2_R1.fastq.gz;
done

cd ~/Documents/MedGenetics/BioinfoAnalysis/fastq_files
ls -alS|grep DNA|sed 's/M//g'|awk '{if($5 < 7.5) print $9}'|cut -d_ -f1|sort|uniq


for i in `ls -alS BioinfoAnalysis/fastq_files/|grep DNA|sed 's/M//g'|awk '{if($5 < 7.5) print $9}'|cut -d_ -f1|sort|uniq`;do grep -w $i qubit.cons.txt ;done
