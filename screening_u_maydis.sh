mkdir /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping
mkdir /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_aDNA_characteristics

./bwa index /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/reference_genome/U.maydis/GCF_000328475.2_Umaydis521_2.0_genomic.fna

samtools faidx /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/reference_genome/U.maydis/GCF_000328475.2_Umaydis521_2.0_genomic.fna

##Mapping against huitlacoche reference genome
for i in {001..020}; do
  echo "Mapeando muestra ZM${i}..."
  
  ./bwa aln -t 8 -l 1024 -f /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/ZM${i}.collapsed.sai \
  /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/reference_genome/U.maydis/GCF_000328475.2_Umaydis521_2.0_genomic.fna \
  /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/2_trimmed_merged/ZM${i}.merged.fastq.gz
  
done

for i in EB106 EB107 LB106 LB107; do
  echo "Mapeando muestra ${i}..."
  
  ./bwa aln -t 8 -l 1024 -f /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/${i}.collapsed.sai \
  /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/reference_genome/U.maydis/GCF_000328475.2_Umaydis521_2.0_genomic.fna \
  /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/2_trimmed_merged/${i}.merged.fastq.gz
  
done

##Generate sam file 
for i in {001..020}; do
  echo "Generando archivo SAM para ZM${i}..."
  
  ./bwa samse -r "@RG\tID:ZM${i}\tSM:ZM${i}" \
  -f /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/ZM${i}.sam \
  /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/reference_genome/U.maydis/GCF_000328475.2_Umaydis521_2.0_genomic.fna \
  /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/ZM${i}.collapsed.sai \
  /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/2_trimmed_merged/ZM${i}.merged.fastq.gz
  
done

for i in EB106 EB107 LB106 LB107; do
  echo "Generando archivo SAM para ${i}..."
  
  ./bwa samse -r "@RG\tID:ZM${i}\tSM:ZM${i}" \
  -f /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/${i}.sam \
  /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/reference_genome/U.maydis/GCF_000328475.2_Umaydis521_2.0_genomic.fna \
  /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/${i}.collapsed.sai \
  /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/2_trimmed_merged/${i}.merged.fastq.gz
  
done

##Calculate proportion of mapped reads
for i in {001..020}; do
  echo "Calculando flagstats para ZM${i}..."
  samtools flagstat -@ 8 /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/ZM${i}.sam > /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/ZM${i}flagstats.log
done

for i in EB106 EB107 LB106 LB107; do
  echo "Calculando flagstats para ZM${i}..."
  samtools flagstat -@ 8 /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/${i}.sam > /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/${i}flagstats.log
done


## Keep the mapped reads and create a compressed BAM file
## The -F 4 option is for samtools to leave behind the reads that are not mapped

for i in {001..020}; do
  echo "Filtrando y convirtiendo ZM${i} a BAM..."
  samtools view -@ 8 -F 4 -Sbh -o /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/ZM${i}.mapped.bam \
  /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/ZM${i}.sam
done

for i in EB106 EB107 LB106 LB107; do
  echo "Filtrando y convirtiendo ${i} a BAM..."
  samtools view -@ 8 -F 4 -Sbh -o /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/${i}.mapped.bam \
  /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/${i}.sam
done


## Sort the BAM file by chromosome and position

for i in {001..020}; do
  echo "Ordenando lecturas para ZM${i}..."
  samtools sort -@ 8 -o /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/ZM${i}.mapped.sorted.bam \
  /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/ZM${i}.mapped.bam
done

for i in EB106 EB107 LB106 LB107; do
  echo "Ordenando lecturas para ZM${i}..."
  samtools sort -@ 8 -o /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/${i}.mapped.sorted.bam \
  /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/${i}.mapped.bam
done


## Remove SAM and unsorted BAM files
## Latorre recommends working with BAM sorted files for optimization of disk storage space

for i in {001..020}; do
  echo "Liberando espacio de ZM${i}..."
  rm -f /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/ZM${i}.sam
  rm -f /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/ZM${i}.mapped.bam
done

for i in EB106 EB107 LB106 LB107; do
  echo "Liberando espacio de ZM${i}..."
  rm -f /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/${i}.sam
  rm -f /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/${i}.mapped.bam
done


## Filter out PCR optical duplicates 

for i in {001..020}; do
  echo "Removiendo duplicados de ZM${i}..."
  java -jar /home/jose-carlos-moreno-juarez/Descargas/DeDup-0.12.9.jar \
  -i /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/ZM${i}.mapped.sorted.bam -m \
  -o /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/
done

for i in EB106 EB107 LB106 LB107; do
  echo "Removiendo duplicados de ZM${i}..."
  java -jar /home/jose-carlos-moreno-juarez/Descargas/DeDup-0.12.9.jar \
  -i /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/${i}.mapped.sorted.bam -m \
  -o /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/
done


## Analyze ancient DNA associated characteristics
## I had to install these libraries for the bayesian analysis

for i in {001..020}; do
  echo "Procesando muestra ZM${i} en mapDamage..."
  
  # Genera las subcarpetas de cada muestra dentro de tu directorio ya creado
  mkdir -p /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_aDNA_characteristics/ZM${i}
  
  # mapDamage apuntando al genoma del huitlacoche
  mapDamage -i /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/ZM${i}.mapped.sorted_rmdup.bam \
  -r /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/reference_genome/U.maydis/GCF_000328475.2_Umaydis521_2.0_genomic.fna \
  -d /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_aDNA_characteristics/ZM${i}
done

for i in EB106 EB107 LB106 LB107; do
  echo "Procesando muestra ZM${i} en mapDamage..."
  
  # Genera las subcarpetas de cada muestra dentro de tu directorio ya creado
  mkdir -p /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_aDNA_characteristics/${i}
  
  # mapDamage apuntando al genoma del huitlacoche
  mapDamage -i /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/${i}.mapped.sorted_rmdup.bam \
  -r /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/reference_genome/U.maydis/GCF_000328475.2_Umaydis521_2.0_genomic.fna \
  -d /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_aDNA_characteristics/${i}
done


## Rerun for rescaling of the plots

for i in {001..020}; do
  echo "Procesando muestra ZM${i} en mapDamage..."
  
  # Genera las subcarpetas de cada muestra dentro de tu directorio ya creado
  mkdir -p /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_aDNA_characteristics/ZM${i}
  
  # mapDamage apuntando al genoma del huitlacoche
  mapDamage -i /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/ZM${i}.mapped.sorted_rmdup.bam \
  -r /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/reference_genome/U.maydis/GCF_000328475.2_Umaydis521_2.0_genomic.fna \
  -d /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_aDNA_characteristics/ZM${i} -y 0.05 --plot-only
done

for i in EB106 EB107 LB106 LB107; do
  echo "Procesando muestra ZM${i} en mapDamage..."
  
  # Genera las subcarpetas de cada muestra dentro de tu directorio ya creado
  mkdir -p /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_aDNA_characteristics/${i}
  
  # mapDamage apuntando al genoma del huitlacoche
  mapDamage -i /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/${i}.mapped.sorted_rmdup.bam \
  -r /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/reference_genome/U.maydis/GCF_000328475.2_Umaydis521_2.0_genomic.fna \
  -d /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_aDNA_characteristics/${i} -y 0.05 --plot-only
done



##Filter out the reads of quality mapping less than 1
for i in {001..020}; do
  echo "Filtrando muestra ZM${i}..."
  
  samtools view -@ 8 -b -q 1 \
  -o /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/ZM${i}.filtered.bam \
  /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/ZM${i}.mapped.sorted_rmdup.bam
  
done


for i in EB106 EB107 LB106 LB107; do
  echo "Filtrando blanco ${i}..."
  
  samtools view -@ 8 -b -q 1 \
  -o /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/${i}.filtered.bam \
  /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/${i}.mapped.sorted_rmdup.bam
  
done


##Create directory with flagstats of filtered by mapping quality genome
mkdir /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/Filtered_mapping_quality




##Calculate the flagstats of the data sets with the filter by mapping quality
for i in {001..020}; do
  echo "Corriendo flagstat en ZM${i}..."
  
  samtools flagstat /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/ZM${i}.filtered.bam > /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/Filtered_mapping_quality/ZM${i}.filtered.flagstat.txt
  
done


for i in EB106 EB107 LB106 LB107; do
  echo "Corriendo flagstat en el blanco ${i}..."
  
  samtools flagstat /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/${i}.filtered.bam > /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/U.maydis_mapping/Filtered_mapping_quality/${i}.filtered.flagstat.txt
  
done

























































##Keep the mapped read and create a compressed BAM file
##The -F 4 option is for samtools to leave behind the reads that are not mapped

for i in {001..020}; do
  echo "Filtrando y convirtiendo ZM${i} a BAM..."
  samtools view -@ 8 -F 4 -Sbh -o /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/ZM${i}.mapped.bam \
  /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/ZM${i}.sam
done

##Sort the BAM file by chromosome and position

for i in {001..020}; do
  echo "Ordenando lecturas para ZM${i}..."
  
  samtools sort -@ 8 -o /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/ZM${i}.mapped.sorted.bam \
  /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/ZM${i}.mapped.bam
  
done


##Remove SAM and unsorted BAM files
##Latorre recommends working with BAM sorted files for optimization of disk storage space
for i in ZM{001..020} EB106 EB107 LB106 LB107; do
  echo "Liberando espacio de ${i}..."
  
  rm -f /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/${i}.sam
  rm -f /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/${i}.mapped.bam
  
done


##Filter out PCR optical duplicates 

 for i in {001..020}; do
  echo "Removiendo duplicados de ZM${i}..."
  
  java -jar /home/jose-carlos-moreno-juarez/Descargas/DeDup-0.12.9.jar \
  -i /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/ZM${i}.mapped.sorted.bam -m \
  -o /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/
  
done

##Analyze ancient DNA associated characteristics
##I had to install these libraries for the bayesina analysis

for i in {001..020}; do
  echo "Procesando muestra ZM${i}..."
  
  #Create directories
  mkdir -p /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/5_aDNA_characteristics/ZM${i}
  
  #MapDamage
  mapDamage -i /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/4_mapping/ZM${i}.mapped.sorted_rmdup.bam \
  -r /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/reference_genome/unmasked.fa \
  -d /home/jose-carlos-moreno-juarez/Documentos/Maestria/screening_teocintle/5_aDNA_characteristics/ZM${i}
  
done
