#!/usr/bin/env bash

mkdir PCA/

# filter the VCF file:
vcftools --gzvcf vcardui.vcf.gz \
         --max-missing 0.9 \
         --recode \
         --recode-INFO-all \
         --out PCA/vcardui_maxmiss0.9

cd PCA/

# run PLINK:
plink --vcf vcardui_maxmiss0.9.recode.vcf \
      --double-id \
      --allow-extra-chr \
      --set-missing-var-ids @:# \
      --indep-pairwise 50 10 0.1 \
      --out vcardui_maxmiss0.9

plink --vcf vcardui_maxmiss0.9.recode.vcf \
      --double-id \
      --allow-extra-chr \
      --set-missing-var-ids @:# \
      --extract vcardui_maxmiss0.9.prune.in \
      --make-bed \
      --pca \
      --out vcardui_maxmiss0.9
