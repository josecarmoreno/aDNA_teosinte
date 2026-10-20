##I downloaded the data using the FileZilla software

##Go to the directory where you have your genomic data
##Check if they were downloaded correctly using the checksum files .md5
cd /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/EGP060_26
cat md5/*.md5 | md5sum -c -

#Create directory for the reference genome
mkdir /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/reference_genome
#Download reference genome from MaizeGBD

##Create the other working directories for the rest of the outputs
mkdir /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/2_trimmed_merged
mkdir /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/3_quality_control
mkdir /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/4_mapping
mkdir /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/5_aDNA_characteristics
mkdir /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/6_filtering

##Remove adapters
s_num=1
for i in {001..020}; do
    echo "done with sample ZM${i}"

    ./adapterremoval3 \
        --in-file1 /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/EGP060_26/EGP060_26_${i}_S${s_num}_L003_R1_001.fastq.gz \
        --in-file2 /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/EGP060_26/EGP060_26_${i}_S${s_num}_L003_R2_001.fastq.gz \
        --out-prefix /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/2_trimmed_merged/ZM${i} \
        --merge \
        --out-format fastq.gz \
        --threads 8
        
    s_num=$((s_num + 1))
done

##Performing a quality trim only in the merged reads
##We eliminate bases in the extremes with a Phred score less than 30
##We eliminate reads shorter than 30 bp
for i in {001..020}; do
    echo "Processing sample ZM${i}..."

    ./adapterremoval3 \
        --in-file1 /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/2_trimmed_merged/ZM${i}.merged.fastq.gz \
        --out-prefix /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/2_trimmed_merged/ZM${i}_qfiltered \
        --quality-trimming mott \
        --trim-min-quality 30 \
        --min-length 30 \
        --out-format fastq.gz \
        --threads 8

done



##Check the quality of reads
for i in {001..020}; do
  echo "done with sample ZM${i}..."
  
  ./fastqc -t 8 --outdir /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/3_quality_control/wo_quality_filters \
    /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/2_trimmed_merged/ZM${i}.merged.fastq.gz \
    /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/2_trimmed_merged/ZM${i}.r1.fastq.gz \
    /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/2_trimmed_merged/ZM${i}.r2.fastq.gz
    
done

##Check the quality of filtered reads
for i in {001..020}; do
  echo "done with sample ZM${i}..."
  
  ./fastqc -t 8 --outdir /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/3_quality_control/ \
    /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/2_trimmed_merged/ZM${i}_qfiltered.r1.fastq.gz
    
done

##Indexing the reference genome
bwa index /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/reference_genome/Zx-TIL18-REFERENCE-PanAnd-1.0.fa.gz

samtools faidx /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/reference_genome/Zx-TIL18-REFERENCE-PanAnd-1.0.fa.gz


##Map trimmed and merged reads to the indexed reference genome
for i in {001..020}; do
  echo "done with sample ZM${i}..."
  
  bwa aln -t 8 -l 1024 -f /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/4_mapping/ZM${i}.collapsed.sai \
  /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/reference_genome/Zx-TIL18-REFERENCE-PanAnd-1.0.fa.gz \
  /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/2_trimmed_merged/ZM${i}_qfiltered.r1.fastq.gz
  
done


##Convert mapped reads into a standars aligment format (SAM)
for i in {001..020}; do
  echo "done with sample ZM${i}..."
  
  bwa samse -r "@RG\tID:ZM${i}\tSM:ZM${i}" \
  -f /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/4_mapping/ZM${i}.sam \
  /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/reference_genome/Zx-TIL18-REFERENCE-PanAnd-1.0.fa.gz \
  /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/4_mapping/ZM${i}.collapsed.sai \
  /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/2_trimmed_merged/ZM${i}.merged.fastq.gz
  
done


##Now we use samtools flagstat to calculate the proportion of mapped reads as a proxy for endogenous DNA
for i in {001..020}; do
  echo "flagstats for ZM${i}..."
  samtools flagstat -@ 8 /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/4_mapping/ZM${i}.sam > /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/4_mapping/ZM${i}flagstats.log
done


##Keep the mapped read and create a compressed BAM file
##The -F 4 option is for samtools to leave behind the reads that are not mapped
for i in {001..020}; do
  echo "Transforming ZM${i} to BAM..."
  samtools view -@ 8 -F 4 -Sbh -o /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/4_mapping/ZM${i}.mapped.bam \
  /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/4_mapping/ZM${i}.sam
done


##Sort the BAM file by chromosome and position
for i in {001..020}; do
  echo "done with sample ZM${i}..."
  
  samtools sort -@ 8 -o /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/4_mapping/ZM${i}.mapped.sorted.bam \
  /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/4_mapping/ZM${i}.mapped.bam
  
done


##Remove SAM and unsorted BAM files
##Latorre recommends working with BAM sorted files for optimization of disk storage space
for i in ZM{001..020}; do
  echo "done with files from ZM${i}..."
  
  rm -f /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/4_mapping/${i}.sam
  rm -f /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/4_mapping/${i}.mapped.bam
  
done


##Eliminate pcr duplicates
for i in {001..020}; do
  echo "processing sample ZM${i}..."
  
  dedup \
  -i /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/4_mapping/ZM${i}.mapped.sorted.bam -m \
  -o /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/4_mapping/
  
done



##Analyze ancient DNA associated characteristics
##I had to install these libraries for the bayesina analysis
R install.packages(c("gam", "RcppGSL"))

for i in {001..020}; do
  echo "processing sample ZM${i}..."
  

  mkdir -p /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/5_aDNA_characteristics/ZM${i}
  
 
  mapDamage -i /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/4_mapping/ZM${i}.mapped.sorted_rmdup.bam \
  -r /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/reference_genome/Zx-TIL18-REFERENCE-PanAnd-1.0.fa.gz \
  -d /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/5_aDNA_characteristics/ZM${i}
  
done

##I did this just to rescale the Y axis to be at 0.05, this dont re run the analysis, only does the rescaling
##This doen NOT rerun the bayesian analysis
for i in {001..020}; do
  echo "rescaling ZM${i}..."
  

  mkdir -p /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/5_aDNA_characteristics/ZM${i}
  

  mapDamage -i /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/4_mapping/ZM${i}.mapped.sorted_rmdup.bam \
  -r /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/reference_genome/Zx-TIL18-REFERENCE-PanAnd-1.0.fa \
  -d /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/5_aDNA_characteristics/ZM${i} -y 0.05 --plot-only
  
done


##Filter out the reads of quality mapping less than 30
for i in {001..020}; do
  echo "filtering sample ZM${i}..."
  
  samtools view -@ 8 -b -q 30 \
  -o /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/6_filtering/ZM${i}.filtered_30q.bam \
  /home/genconser/Documentos/JoseCarlos/paleogenomic_Zm/4_mapping/ZM${i}.mapped.sorted_rmdup.bam
  
done
