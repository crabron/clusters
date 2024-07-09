checkm2 predict --input ../../rnf_c/h_mags \
	--allmodels \
	-t 100 \
	--database_path \
	~/somebases/checkm2/CheckM2_database/uniref100.KO.1.dmnd \
	-o checkm2 \
	-x fa


medaka_consensus -d 1/assembly.fasta -i 1/all.fastq.gz  -t 20 -m r941_min_hac_g507 -b 20 -o  1/medaka_consensus &
medaka_consensus -d 2/assembly.fasta -i 2/all.fastq.gz  -t 20 -m r941_min_hac_g507  -b 20-o  2/medaka_consensus &
medaka_consensus -d 3/assembly.fasta -i 3/all.fastq.gz  -t 20 -m r941_min_hac_g507  -b 20 -o  3/medaka_consensus &
medaka_consensus -d 4/assembly.fasta -i 4/all.fastq.gz  -t 20 -m r941_min_hac_g507  -b 20 -o  4/medaka_consensus &
medaka_consensus -d 5/assembly.fasta -i 5/all.fastq.gz  -t 20 -m r941_min_hac_g507  -b 20 -o  5/medaka_consensus 



for fname in ; do
    medaka_consensus \
        -d ${fname}/assembly.fasta \
        -i ../called_sup/${fname}.*/all.fastq.gz \
        -t 100 \
        -m r941_min_sup_g507 \
        -o  ${fname}/medaka_consensus 
done;


for fname in gb_3_path.txt; do
    medaka_consensus \
        -d ${fname}/assembly.fasta \
        -i ../called_sup/${fname}.*/all.fastq.gz \
        -t 100 \
        -m r941_min_sup_g507 \
        -o  ${fname}/medaka_consensus 
done;

for fname in read gb_3_path.txt; do
	echo $fname
done;

cat gb_3_path.txt | while read line; do
        du -h $line
        reads=$(basename "$line" .fa) | grep -Po '.*(?=_)'
        echo $reads
done;

cat gb_3_path.txt | while read line; do
        du -h $line
        fname=$(basename "$line" .fa | grep -Po '.*(?=_)')
        id=$(basename "$line" .fa)
        medaka_consensus \
		-d $line \
		-i in/reads/${fname}* \
		-t 100 \
		-b 30 \
		-m r941_min_hac_g507 \
		-o  medaka_consensus_3/${id}
done;
cat gb_4_path.txt | while read line; do
        du -h $line
        fname=$(basename "$line" .fa | grep -Po '.*(?=_)')
        id=$(basename "$line" .fa)
        medaka_consensus \
		-d $line \
		-i in/reads/${fname}* \
		-t 100 \
		-b 30 \
		-m r941_min_hac_g507 \
		-o  medaka_consensus_4/${id}
done;

for dir in semibin_results_4/bins/*; do
        fname=$(basename "$dir")
      	cp $dir/consensus.fasta h_bins_3/${fname}.fasta
done;

for dir in semibin_results_4/bins/c1_bin.*.fa; do
	echo $dir
done;

for file in semibin_results_4/bins/*.fa; do 
	fname=$(basename "$dir" .fa)
	medaka_consensus \
		-d $file \
		-i in/reads/${fname}.fq.gz \
		-t 15 \
		-b 8 \
		-m r941_min_hac_g507 \
		-o  medaka_consensus_4_all/$fname &
	wait -n 10
done; 

for file in bins_medium/*.fasta; do 
	fname=$(basename "$file" .fasta)
	barrnap \
		--threads 30 \
		--outseq barrnap/${fname}.fasta \
		$file
done; 


		
#!/bin/bash		

for file in semibin_results_4/bins/*.fa; do
  if [ $(jobs -r | wc -l) -ge 10 ]; then
    wait $(jobs -r -p | head -1)
  fi

	echo Begin processing $file         
	fname=$(basename "$dir" .fa)
	medaka_consensus \
		-d $file \
		-i in/reads/${fname}.fq.gz \
		-t 15 \
		-b 8 \
		-m r941_min_hac_g507 \
		-o  medaka_consensus_4_all/$fname &
	echo End processing $file
done
wait


for file in ../in/sorted_bam_2/*.bam; do
	fname=$(basename $file .bam)
	echo $fname
	jgi_summarize_bam_contig_depths --outputDepth  ${fname}_depth.txt $file &
done;

minimap2 -x map-ont -t 50 -a assembly.fasta ../assembled_by_repeates/1/all.fastq.gz | samtools sort -o al.bam --write-index -

for file in ass/*.fasta; do
	fname=$(basename $file .fasta)
	echo $fname
	minimap2 -x map-ont -t 30 -a $file reads/${fname}.fq.gz | samtools sort -o sorted_bam_2/${fname}.bam --write-index - &
done;


for file in semibin_results_4/bins/*.fa; do
  if [ $(jobs -r | wc -l) -ge 2 ]; then
    wait $(jobs -r -p | head -1)
  fi

	echo Begin processing $file         
	fname=$(basename "$dir" .fa)
	medaka_consensus \
		-d $file \
		-i in/reads/${fname}.fq.gz \
		-t 15 \
		-b 8 \
		-m r941_min_hac_g507 \
		-o  medaka_consensus_4_all/$fname &
		
	gtdbtk classify_wf \
		--genome_dir bins_medium \
		--out_dir gtdb_medium \
		-x fasta \
		--mash_db ~/somebases/gtdb/mash \
		--cpus 100 \
		--full_tree \
		--pplacer_cpus 50 \
		--write_single_copy_genes 
	
	echo End processing $file
done
wait


while read -r line; do
	cp bins_4_all/${line}.fasta bins_medium
done < bins_medium.txt

grep ">" l3_bin.16.fa | sed 's/>/l3:/' | sed 's/$/,/' | grep -f  -  ../samples/l3.bam_5_data_cov.csv | sed 's/^.*,//' | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }'

