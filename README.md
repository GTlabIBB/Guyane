# Guyane

## Metabarcoding analyses

The raw metabarcoding sequences are stored in the European Nucleotide Archive project PRJEB57731 (https://www.ebi.ac.uk/ena/data/view/PRJEB57731).

The script to process and classify the raw reads is in `Metabarcoding/metabarcoding_script.sh` and takes the data stored in the directory called `1-rawreads`. Reference databases should be in the path: `ITS2_DBs/Sickel_et_al_2015/viridiplantae_all_2014.sintax.fa` and `ITS2_DBs/PLANiTS_29-03-2020/ITS2.SINTAX_format.fas`.

The script generates SINTAX classifications using two different reference databases:
1. Sickel et al. (2015) - output in `classification_sickel` directory,
2. the PLANiTS database - output in `classification_planits` directory,

of the reads clustered at:
1. 98% sequence similarity (`otus0.98.sintax.txt` in the abovementioned directories),
2. 99% sequence similarity (`otus0.99.sintax.txt` in the abovementioned directories),
3. denoised using UNOISE algorithm implemented in VSEARCH (`zotus.sintax.txt` in the abovementioned directories).


For the details please refer to the accompanying publication.


The R scripts used to summarize the classification of the processed reads are in `Metabarcoding/metabarcoding_summarize-classification-sickel.R` and `Metabarcoding/metabarcoding_summarize-classification-planits.R` and use the data generated by `Metabarcoding/metabarcoding_script.sh` in `classification_sickel/zotus.sintax.txt` and  `classification_planits/zotus.sintax.txt`, respectively. The script can be changed to process `otus0.98.sintax.txt` and `otus0.99.sintax.txt` files from the above directories as well.


## Population genetics analyses

Demultiplexed and adapter-trimmed RAD-sequencing reads are available in the European Nucleotide Archive project PRJEB57763 (https://www.ebi.ac.uk/ena/data/view/PRJEB57763).

Allele output of the ipyrad analysis (`vcardui.alleles`) is processed in fineRADstructure using commands provided in `Population_genetics/fineRADstructure_script.sh`. Please note that for the R script provided by the fineRADstructure to plot the results the samples must be renamed to follow "Country+number" format (eg. Spain1, Spain2, etc.).

The VCF file form ipyrad (`vcardui.vcf.gz`) is filtered using vcftools and used to run PCA in PLINK with the `Population_genetics/PCA_script.sh`.

## Wind trajectory analyses

The scripts for the wind trajectory calculations can be found in the `Hysplit_wind_analysis` folder. The script `guyane_wind_trajectories-one-day_blocks.sh` was used to generate the daily trajectory plots, whereas the script `guyane_wind_trajectories_two-day_blocks.sh` can be used to generate the two-day plots. To run these scripts, you will need to install [this GitHub repository](https://github.com/etd530/Hysplit_R_interface) as well as [GNU Parallel](https://www.gnu.org/software/parallel/) (the latter is optional, but the scripts need to be slightly modified to be ran without it).

After obtaining the trajectory calculations (`.RData` files), the scripts `wind_speed_plots.R` and `mean_speed_per_traj.R` can be used to obtain the figures and table on wind speed and mean trajectory speeds, respectively. Lastly, `Rdata_to_shp.R` and `get_intersection_hours.R` can be used to obtain the table containing times to intersect the African coastline and the proportion of trajectories that intersect it. This step requires a shapefile available in the Source Data of the article.

## Isotope-based geographic assignment

The R markdown file, `isotope_script_guiana_FINAL`, that contains the script for the isotope-based geographic assignment can be found in the `Isotope_Rscript` folder. 
