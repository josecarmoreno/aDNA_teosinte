###This scripts is the one i used to perform my screening of teosinte samples
###It is based on the pipeline supplied by 
#Latorre, S. M., Lang, P. L. M., Burbano, H. A., & Gutaker, R. M.
#(2020). Isolation, library preparation, and bioinformatic analysis of
#historical and ancient plant DNA. Current Protocols in Plant
#Biology, 5, e20121. doi: 10.1002/cppb.20121

#First i had to download the following softwares
#Im using Ubuntu24
#AdapterRemoval v2.3.1 (Schubert, Lindgreen, & Orlando, 2016; https:// github.com/ mikkelschubert/ adapterremoval)
#BWA v. 0.7.17 (Li, 2013; https:// github.com/ lh3/ bwa)
#FastQC v. 0.11.9 (Andrews, n.d.; https:// github.com/ s-andrews/ FastQC)
#samtools v. 1.10 (Li et al., 2009; https:// github.com/ samtools/ samtools)
#DeDup v. 0.12.6 (Peltzer et al., 2016; https:// github.com/ apeltzer/ DeDup)
#mapDamage v. 2.2.1 (Jónsson, Ginolhac, Schubert, Johnson, & Orlando, 2013; https:// github.com/ ginolhac/ mapDamage)

##The reference genome used is the one of Zea mays ssp. mays inbred line B73 Zm-B73-REFERENCE-NAM-5.0 (Assembly accesion GCA_902167145.1). Released in july 2026 by the NAM sequencing consortium
##The reference genome was generated with PacBio long-read sequencing, polished with Illumina reads, and Bionano optical map technology was used for scaffold assembly
#Create directory for the reference genome
mkdir /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/reference_genome
#Download reference genome
wget https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/902/167/145/1/community_cshl/2019_12/genome/unmasked.fa.bgz

##Create the other working directories for the rest of the outputs
mkdir /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/2_trimmed_merged
mkdir /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/3_quality_control
mkdir /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping
mkdir /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/5_aDNA_characteristics


##First I need to remove the adapters and merge pair mates from the raw Illumina reads
##Some of the options were renamed in latter versions in respect with LAtorre pipeline, i just updated the names, but the underlying script is the same
for i in {001..020}; do
  echo "Arrancando con la muestra ZM${i}..."
  
  ./adapterremoval3 \
    --in-file1 /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/data_aDNA_teocintle/ZM${i}_1.fastq.gz \
    --in-file2 /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/data_aDNA_teocintle/ZM${i}_2.fastq.gz \
    --out-prefix /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/2_trimmed_merged/ZM${i} \
    --merge \
    --out-format fastq.gz \
    --threads 8
    
done
  
##EB106 has 12,426 reads
  ./adapterremoval3 \
  --in-file1 /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/data_aDNA_teocintle/EB106_1.fastq.gz \
  --in-file2 /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/data_aDNA_teocintle/EB106_2.fastq.gz \
  --out-prefix /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/2_trimmed_merged/EB106 \
  --merge \
  --out-format fastq.gz \
  --threads 8
  
##EB107 has 34,648 reads  
    ./adapterremoval3 \
  --in-file1 /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/data_aDNA_teocintle/EB107_1.fastq.gz \
  --in-file2 /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/data_aDNA_teocintle/EB107_2.fastq.gz \
  --out-prefix /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/2_trimmed_merged/EB107 \
  --merge \
  --out-format fastq.gz \
  --threads 8

##LB106 has 4,738 reads
  ./adapterremoval3 \
  --in-file1 /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/data_aDNA_teocintle/LB106_1.fastq.gz \
  --in-file2 /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/data_aDNA_teocintle/LB106_2.fastq.gz \
  --out-prefix /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/2_trimmed_merged/LB106 \
  --merge \
  --out-format fastq.gz \
  --threads 8
  
##LB107 has 9,946 reads
  ./adapterremoval3 \
  --in-file1 /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/data_aDNA_teocintle/LB107_1.fastq.gz \
  --in-file2 /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/data_aDNA_teocintle/LB107_2.fastq.gz \
  --out-prefix /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/2_trimmed_merged/LB107 \
  --merge \
  --out-format fastq.gz \
  --threads 8
  

##I need to check the quality of the trimmed and merged and the trimmed and non-merged
reads
for i in {001..020}; do
  echo "Corriendo FastQC para ZM${i}..."
  
  ./fastqc -t 8 -o /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/3_quality_control/ \
    /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/2_trimmed_merged/ZM${i}.merged.fastq.gz \
    /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/2_trimmed_merged/ZM${i}.r1.fastq.gz \
    /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/2_trimmed_merged/ZM${i}.r2.fastq.gz
    
done
  
  
 ./fastqc -t 8 -o /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/3_quality_control/ \
    /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/2_trimmed_merged/EB106.merged.fastq.gz \
    /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/2_trimmed_merged/EB106.r1.fastq.gz \
    /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/2_trimmed_merged/EB106.r2.fastq.gz
    
    
./fastqc -t 8 -o /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/3_quality_control/ \
    /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/2_trimmed_merged/EB107.merged.fastq.gz \
    /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/2_trimmed_merged/EB107.r1.fastq.gz \
    /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/2_trimmed_merged/EB107.r2.fastq.gz
    
./fastqc -t 8 -o /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/3_quality_control/ \
    /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/2_trimmed_merged/LB106.merged.fastq.gz \
    /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/2_trimmed_merged/LB106.r1.fastq.gz \
    /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/2_trimmed_merged/LB106.r2.fastq.gz
    
    ./fastqc -t 8 -o /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/3_quality_control/ \
    /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/2_trimmed_merged/LB107.merged.fastq.gz \
    /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/2_trimmed_merged/LB107.r1.fastq.gz \
    /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/2_trimmed_merged/LB107.r2.fastq.gz
  
  
##Indexing the reference genome
./bwa index /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/reference_genome/unmasked.fa

samtools faidx /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/reference_genome/unmasked.fa

##Map trimmed and merged reads to the indexed reference genome
./bwa aln -t 8 -l 1024 -f /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/ZM001.collapsed.sai \ 
/home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/reference_genome/unmasked.fa \
/home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/2_trimmed_merged/ZM001.merged.fastq.gz

for i in {002..020}; do
  echo "Mapeando muestra ZM${i}..."
  
  ./bwa aln -t 8 -l 1024 -f /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/ZM${i}.collapsed.sai \
  /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/reference_genome/unmasked.fa \
  /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/2_trimmed_merged/ZM${i}.merged.fastq.gz
  
done

./bwa aln -t 8 -l 1024 -f /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/EB106.collapsed.sai \
/home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/reference_genome/unmasked.fa \
/home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/2_trimmed_merged/EB106.merged.fastq.gz

./bwa aln -t 8 -l 1024 -f /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/EB107.collapsed.sai \
/home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/reference_genome/unmasked.fa \
/home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/2_trimmed_merged/EB107.merged.fastq.gz

./bwa aln -t 8 -l 1024 -f /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/LB106.collapsed.sai \
/home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/reference_genome/unmasked.fa \
/home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/2_trimmed_merged/LB106.merged.fastq.gz

./bwa aln -t 8 -l 1024 -f /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/LB107.collapsed.sai \
/home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/reference_genome/unmasked.fa \
/home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/2_trimmed_merged/LB107.merged.fastq.gz

##Convert mapped reads into a standars aligment format (SAM)
./bwa samse -r @RG\\tID:ZM001\\tSM:ZM001 -f /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/ZM001.sam \
/home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/reference_genome/unmasked.fa \
/home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/ZM001.collapsed.sai \
/home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/2_trimmed_merged/ZM001.merged.fastq.gz

for i in {002..020}; do
  echo "Generando archivo SAM para ZM${i}..."
  
  ./bwa samse -r "@RG\tID:ZM${i}\tSM:ZM${i}" \
  -f /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/ZM${i}.sam \
  /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/reference_genome/unmasked.fa \
  /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/ZM${i}.collapsed.sai \
  /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/2_trimmed_merged/ZM${i}.merged.fastq.gz
  
done


for i in EB106 EB107 LB106 LB107; do
  echo "Generando archivo SAM para el blanco ${i}..."
  
  ./bwa samse -r "@RG\tID:${i}\tSM:${i}" \
  -f /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/${i}.sam \
  /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/reference_genome/unmasked.fa \
  /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/${i}.collapsed.sai \
  /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/2_trimmed_merged/${i}.merged.fastq.gz
  
done

##Now we use samtools flagstat to calculate the proportion of mapped reads as a proxy for endogenous DNA
for i in {001..020}; do
  echo "Calculando flagstats para ZM${i}..."
  samtools flagstat -@ 8 /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/ZM${i}.sam > /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/ZM${i}flagstats.log
done

for i in EB106 EB107 LB106 LB107; do
  echo "Calculando flagstats para blanco ${i}..."
  samtools flagstat -@ 8 /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/${i}.sam > /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/${i}flagstats.log
done

##Keep the mapped read and create a compressed BAM file
##The -F 4 option is for samtools to leave behind the reads that are not mapped
samtools view -@ 8 -F 4 -Sbh -o /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/ZM001.mapped.bam \
/home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/ZM001.sam

for i in {002..020}; do
  echo "Filtrando y convirtiendo ZM${i} a BAM..."
  samtools view -@ 8 -F 4 -Sbh -o /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/ZM${i}.mapped.bam \
  /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/ZM${i}.sam
done

for i in EB106 EB107 LB106 LB107; do
  echo "Filtrando y convirtiendo blanco ${i} a BAM..."
  samtools view -@ 8 -F 4 -Sbh -o /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/${i}.mapped.bam \
  /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/${i}.sam
done

##Sort the BAM file by chromosome and position
samtools sort -@ 8 -o /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/ZM001.mapped.sorted.bam \
/home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/ZM001.mapped.bam

for i in {002..020}; do
  echo "Ordenando lecturas para ZM${i}..."
  
  samtools sort -@ 8 -o /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/ZM${i}.mapped.sorted.bam \
  /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/ZM${i}.mapped.bam
  
done

for i in EB106 EB107 LB106 LB107; do
  echo "Ordenando lecturas para el blanco ${i}..."
  
  samtools sort -@ 8 -o /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/${i}.mapped.sorted.bam \
  /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/${i}.mapped.bam
  
done

##Remove SAM and unsorted BAM files
##Latorre recommends working with BAM sorted files for optimization of disk storage space
for i in ZM{001..020} EB106 EB107 LB106 LB107; do
  echo "Liberando espacio de ${i}..."
  
  rm -f /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/${i}.sam
  rm -f /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/${i}.mapped.bam
  
done

##Filter out PCR optical duplicates 

  java -jar /home/jose-carlos-moreno-juarez/Descargas/DeDup-0.12.9.jar \
  -i /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/ZM001.mapped.sorted.bam -m \
  -o /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/
  
 for i in {002..020}; do
  echo "Removiendo duplicados de ZM${i}..."
  
  java -jar /home/jose-carlos-moreno-juarez/Descargas/DeDup-0.12.9.jar \
  -i /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/ZM${i}.mapped.sorted.bam -m \
  -o /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/
  
done

for i in EB106 EB107 LB106 LB107; do
  echo "Removiendo duplicados del blanco ${i}..."
  
  java -jar /home/jose-carlos-moreno-juarez/Descargas/DeDup-0.12.9.jar \
  -i /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/${i}.mapped.sorted.bam -m \
  -o /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/
  
done

##Analyze ancient DNA associated characteristics
##I had to install these libraries for the bayesina analysis
R install.packages(c("gam", "RcppGSL"))

for i in {001..020}; do
  echo "Procesando muestra ZM${i}..."
  
  #Create directories
  mkdir -p /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/5_aDNA_characteristics/ZM${i}
  
  #MapDamage
  mapDamage -i /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/ZM${i}.mapped.sorted_rmdup.bam \
  -r /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/reference_genome/unmasked.fa \
  -d /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/5_aDNA_characteristics/ZM${i}
  
done

for i in EB106 EB107 LB106 LB107; do
  echo "Procesando blanco ${i}..."
  
  # 1. Crear el directorio específico para el blanco
  mkdir -p /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/5_aDNA_characteristics/${i}
  
  # 2. Correr mapDamage y mandar el output (-d) al nuevo directorio
  mapDamage -i /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/${i}.mapped.sorted_rmdup.bam \
  -r /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/reference_genome/unmasked.fa \
  -d /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/5_aDNA_characteristics/${i}
  
done


##I did this just to rescale the Y axis to be at 0.05, this dont re run the analysis, only does the rescaling
##This doen NOT rerun the bayesian analysis

for i in {001..020}; do
  echo "Procesando muestra ZM${i}..."
  
  #Create directories
  mkdir -p /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/5_aDNA_characteristics/ZM${i}
  
  #MapDamage
  mapDamage -i /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/ZM${i}.mapped.sorted_rmdup.bam \
  -r /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/reference_genome/unmasked.fa \
  -d /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/5_aDNA_characteristics/ZM${i} -y 0.05 --plot-only
  
done

for i in EB106 EB107 LB106 LB107; do
  echo "Procesando blanco ${i}..."
  
  # 1. Crear el directorio específico para el blanco
  mkdir -p /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/5_aDNA_characteristics/${i}
  
  # 2. Correr mapDamage y mandar el output (-d) al nuevo directorio
  mapDamage -i /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/${i}.mapped.sorted_rmdup.bam \
  -r /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/reference_genome/unmasked.fa \
  -d /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/5_aDNA_characteristics/${i} -y 0.05 --plot-only
  
done
