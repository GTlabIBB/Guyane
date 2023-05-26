#!/usr/bin/env bash

# vcardui.alleles file exported with finerad_input.py from
# https://github.com/edgardomortiz/fineRADstructure-tools:

python finerad_input.py -i vcardui.alleles \
                 -t ipyrad \
                 -n 2 \
                 -o fineRADstr/vcardui.renamed.alleles.min2.finerad

# sample names are then renamed according to fineRADstructure
# requirements as in samplenames-fineRADstructure.txt


fineRADstructure/RADpainter paint fineRADstr/vcardui.renamed.alleles.min2.finerad


fineRADstructure/finestructure -x 100000 \
              -y 100000 \
              -z 1000 \
              fineRADstr/vcardui.renamed.alleles.min2_chunks.out \
              fineRADstr/vcardui.renamed.alleles.min2_chunks.mcmc.xml


fineRADstructure/finestructure -m T \
              -x 10000 \
              fineRADstr/vcardui.renamed.alleles.min2_chunks.out \
              fineRADstr/vcardui.renamed.alleles.min2_chunks.out_chunks.mcmc.xml \
              fineRADstr/vcardui.renamed.alleles.min2_chunks.out_chunks.mcmcTree.xml
