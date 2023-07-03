# gutprotist-search
Snakemake workflow to search for gut protists without genomes in short read metagenomics data. This is currently limited to Iodamoeba, Entamoeba coli, Dientamoeba fragilis, Isospora belli, Retortamonas, Pentatrichomonas hominis, Cystoisospora belli, Cystoisospora belli, Endolimax nana, Enteromonas hominis, Entamoeba hartmanni, and Chilomastix mesnili.

This is an offshoot of EukDetect but is substantially less tested. A fairly substantial number of gut eukaryotes do not have sequenced genomes and are not detectable by EukDetect. To get around this, this tool uses a database of sequences deposited in NCBI for a set of given taxids. There is no way to use this quantitatively and this can be used only for presence/absence determination only. The database has been somewhat cleaned of bacterial sequence but it's not a guarantee that it's free of any contamination, and it's also not a guarantee that all of these sequences were correctly assigned by those who deposited them in NCBI.

## Usage

**Edit the config file**

Copy the `default_configfile.yml` to `your_configfile.yml`. Change all parameters in the config file as described.

If you don't know what the length of your reads is, this is a handy one-liner to estimate it: `gzip -dc {file.fastq.gz} | head -n 10000 | awk '{ if (NR%4==2){count++; bases += length}} END{printf "%3.0f\n", bases/count}'`

**Activate the conda environment**

This pipeline was devleoped alongside EukDetect. You can either install eukdetect completely as described in the [github](https://github.com/allind/EukDetect) or just the conda environment. Regardless, you need to activate the environment.

`conda activate eukdetect`

**Run snakemake directly**

This is a snakemake worflow. Options for running are `runall`, `aln`, or `filter` as the target rule.

Examples:
```
snakemake --snakefile nogenome_protist.rules --configfile [config file] --cores [cores] runall
```

## Database info

This database contains a variable number of genes corresponding to <i> Chilomastix mesnili, Endolimax nana, Entamoeba coli, Entamoeba hartmanni, Enteromonas hominis, Iodamoeba sp., Cystoisospora belli, Pentatrichomonas hominis, Retortamonas sp., and Dientamoeba fragilis</i>, which as of September 2022 do not have sequenced genomes deposited in NCBI. All of these are either gut commensals or of uncertain pathogenicity. Making the database involved downloading nucleotide sequences, de-deduplicating, and then using a number of methods to try to remove bacterial contamination (aligning simulated bacterial reads, aligning real gut metagenomic data and filtering out obviously spurious reads).

Caveats:
- I cannot guarantee there is no possibility bacterial-derived sequences will align to this database
- I cannot guarantee that the sequences deposited in NCBI definitively came from the organism - errors made by the depositor are possible
- There is variable coverage across all of these species, and it's not possible to use a lack of hits to this database as proof a sample does not contain the protist
- This information can't be use to quantitate the amount of protist present in a sample
