#!/usr/bin/env bash

# the raw reads should be in directory called "1-rawreads"
# reference databases should be in the path:
# ITS2_DBs/Sickel_et_al_2015/viridiplantae_all_2014.sintax.fa
# ITS2_DBs/PLANiTS_29-03-2020/ITS2.SINTAX_format.fas

mkdir 2-merged/
mkdir 3-cleaned/
mkdir 4-filtered/
mkdir 5-dereplicated/
mkdir reads-classified_sintax
THREADS=4 # change according to the number of threads available on your machine

#merge reads
for fwd in 1-rawreads/*_R1_001.fastq.gz
do
rev=${fwd/R1/R2}
out=$(basename $fwd | sed 's/_R1_001.fastq.gz//g')
pearRM -j $THREADS -f $fwd -r $rev -o 2-merged/$out  2>&1 | tee 2-merged/$out.log
done

#remove primer sequences
for f in 2-merged/*.assembled.fastq
do
out=$(basename $f .assembled.fastq)
cutadapt -g ATGCGATACTTGGTGTGAAT...GCATATCAATAAGCGGAGGA \
         --discard-untrimmed \
         --revcomp \
         --overlap 42 \
         --report minimal \
         --cores $THREADS \
         -o 3-cleaned/$out.trimmed.fastq $f 2>&1 | tee 3-cleaned/$out.trimmed.log
done

#filter reads
for f in 3-cleaned/*.trimmed.fastq
do
out=$(basename $f| sed 's/.trimmed.fastq//g')
vsearch --threads $THREADS \
        --fastq_filter  $f \
        --fastq_maxee 1 \
        --fastq_minlen 250 \
        --fastq_maxns 0 \
        --fastaout 4-filtered/$out.filtered.fasta \
        --fasta_width 0 2>&1 | tee 4-filtered/$out.filtered.log
done

#dereplicate, discard reads which appear once, rename reads with sample id
for f in 4-filtered/*.filtered.fasta
do
out=$(basename $f| sed 's/.filtered.fasta//g')
label=${out/-/}
vsearch --threads $THREADS \
        --derep_fulllength $f \
        --minuniquesize 2 \
        --strand plus \
        --output 5-dereplicated/$out.derep.fasta \
        --sizeout \
        --uc 5-dereplicated/$out.derep.uc \
        --relabel $label. \
        --fasta_width 0 2>&1 | tee 5-dereplicated/$out.derep.log
done

#merge and dereplicate again
mkdir processing
cat 5-dereplicated/*.derep.fasta > processing/all.fasta

vsearch --derep_fulllength processing/all.fasta \
        --threads $THREADS \
        --sizein \
        --sizeout \
        --fasta_width 0 \
        --uc processing/all.derep.uc \
        --output processing/derep.fasta



mkdir output_pool_otus0.98/
mkdir output_pool_otus0.99/
mkdir output_pool_zotus/

#cluster at 98%
vsearch --cluster_size processing/derep.fasta \
        --threads $THREADS \
        --id 0.98 \
        --strand plus \
        --sizein \
        --sizeout \
        --fasta_width 0 \
        --centroids output_pool_otus0.98/otus.centroids.fasta

#cluster at 99%
vsearch --cluster_size processing/derep.fasta \
        --threads $THREADS \
        --id 0.99 \
        --strand plus \
        --sizein \
        --sizeout \
        --fasta_width 0 \
        --centroids output_pool_otus0.99/otus.centroids.fasta

#denoise
vsearch --cluster_unoise processing/derep.fasta \
        --threads $THREADS \
        --minsize 4 \
        --strand plus \
        --sizein \
        --sizeout \
        --fasta_width 0 \
        --centroids output_pool_zotus/zotus.centroids.fasta

#sort
vsearch --sortbysize output_pool_otus0.98/otus.centroids.fasta \
        --threads $THREADS \
        --sizein \
        --sizeout \
        --fasta_width 0 \
        --minsize 2 \
        --output output_pool_otus0.98/otus.sorted.fasta

vsearch --sortbysize output_pool_otus0.99/otus.centroids.fasta \
        --threads $THREADS \
        --sizein \
        --sizeout \
        --fasta_width 0 \
        --minsize 2 \
        --output output_pool_otus0.99/otus.sorted.fasta

vsearch --sortbysize output_pool_zotus/zotus.centroids.fasta \
        --threads $THREADS \
        --sizein \
        --sizeout \
        --fasta_width 0 \
        --minsize 2 \
        --output output_pool_zotus/zotus.sorted.fasta

# de novo chimera detection
vsearch --uchime_denovo output_pool_otus0.98/otus.sorted.fasta \
        --sizein \
        --sizeout \
        --fasta_width 0 \
        --qmask none \
        --nonchimeras output_pool_otus0.98/otus.nonchimeras.fasta

vsearch --uchime_denovo output_pool_otus0.99/otus.sorted.fasta \
        --sizein \
        --sizeout \
        --fasta_width 0 \
        --qmask none \
        --nonchimeras output_pool_otus0.99/otus.nonchimeras.fasta

vsearch --uchime3_denovo output_pool_zotus/zotus.sorted.fasta \
        --sizein \
        --sizeout \
        --fasta_width 0 \
        --qmask none \
        --nonchimeras output_pool_zotus/zotus.nonchimeras.fasta


#relabel OTUs
vsearch --fastx_filter output_pool_otus0.98/otus.nonchimeras.fasta \
        --threads $THREADS \
        --sizein \
        --sizeout \
        --fasta_width 0 \
        --relabel OTU_ \
        --fastaout output_pool_otus0.98/otus.fasta

vsearch --fastx_filter output_pool_otus0.99/otus.nonchimeras.fasta \
        --threads $THREADS \
        --sizein \
        --sizeout \
        --fasta_width 0 \
        --relabel OTU_ \
        --fastaout output_pool_otus0.99/otus.fasta

vsearch --fastx_filter output_pool_zotus/zotus.nonchimeras.fasta \
        --threads $THREADS \
        --sizein \
        --sizeout \
        --fasta_width 0 \
        --relabel OTU_ \
        --fastaout output_pool_zotus/zotus.fasta

#map sequences to OTUs
vsearch --usearch_global processing/all.fasta \
        --threads $THREADS \
        --db output_pool_otus0.98/otus.fasta \
        --id 0.98 \
        --strand plus \
        --sizein \
        --sizeout \
        --fasta_width 0 \
        --qmask none \
        --dbmask none \
        --otutabout output_pool_otus0.98/otutab.txt 2>&1 | tee output_pool_otus0.98/otutab.log

vsearch --usearch_global processing/all.fasta \
        --threads $THREADS \
        --db output_pool_otus0.99/otus.fasta \
        --id 0.99 \
        --strand plus \
        --sizein \
        --sizeout \
        --fasta_width 0 \
        --qmask none \
        --dbmask none \
        --otutabout output_pool_otus0.99/otutab.txt 2>&1 | tee output_pool_otus0.99/otutab.log

vsearch --usearch_global processing/all.fasta \
        --threads $THREADS \
        --db output_pool_zotus/zotus.fasta \
        --id 0.98 \
        --strand plus \
        --sizein \
        --sizeout \
        --fasta_width 0 \
        --qmask none \
        --dbmask none \
        --otutabout output_pool_zotus/zotutab.txt 2>&1 | tee output_pool_zotus/zotutab.log

# classify OTUs
mkdir classification_planits/

vsearch --sintax output_pool_otus0.98/otus.fasta \
        --db ITS2_DBs/PLANiTS_29-03-2020/ITS2.SINTAX_format.fas \
        --tabbedout classification_planits/otus0.98.sintax.txt \
        --notrunclabels \
        --strand both

vsearch --sintax output_pool_otus0.99/otus.fasta \
        --db ITS2_DBs/PLANiTS_29-03-2020/ITS2.SINTAX_format.fas \
        --tabbedout classification_planits/otus0.99.sintax.txt \
        --notrunclabels \
        --strand both

vsearch --sintax output_pool_zotus/zotus.fasta \
        --db ITS2_DBs/PLANiTS_29-03-2020/ITS2.SINTAX_format.fas \
        --tabbedout classification_planits/zotus.sintax.txt \
        --notrunclabels \
        --strand both


mkdir classification_sickel/

vsearch --sintax output_pool_otus0.98/otus.fasta \
        --db ITS2_DBs/Sickel_et_al_2015/viridiplantae_all_2014.sintax.fa \
        --tabbedout classification_sickel/otus0.98.sintax.txt \
        --notrunclabels \
        --strand both

vsearch --sintax output_pool_otus0.99/otus.fasta \
        --db ITS2_DBs/Sickel_et_al_2015/viridiplantae_all_2014.sintax.fa \
        --tabbedout classification_sickel/otus0.99.sintax.txt \
        --notrunclabels \
        --strand both

vsearch --sintax output_pool_zotus/zotus.fasta \
        --db ITS2_DBs/Sickel_et_al_2015/viridiplantae_all_2014.sintax.fa \
        --tabbedout classification_sickel/zotus.sintax.txt \
        --notrunclabels \
        --strand both
