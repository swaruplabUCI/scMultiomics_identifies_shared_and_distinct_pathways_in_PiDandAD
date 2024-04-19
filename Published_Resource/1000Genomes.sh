#=====================
## 1000 Genomes
#=====================

#=====================
#==== GRCh37 assembly
# website::
https://zenodo.org/record/3359882#.Y5D57OzMITt

# NOTEs from: The 1000 Genomes Project Consortium
# The 1000 Genomes Project Consortium. (2019). 1000 Genomes Project (1.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.3359882
# The final phase of the 1000 Genomes Project was phase 3 and represents 2504 samples on GRCh37.
# The data from phase three of the 1000 Genomes Project was subsequently reanalysed on GRCh38.
# Following this work, the samples have been resequenced to high-coverage, with additional related samples being sequenced, bringing the total number of samples up to 3,202. This data was analysed on GRCh38.

#=====================
#==== GRCh38 assembly
#########
# Paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7059836/
# Variant calling on the GRCh38 assembly with the data from phase three of the 1000 Genomes Project

# http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/
#	 	- 20181203_biallelic_SNV
cd 20181203_biallelic_SNV

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/20181203_biallelic_SNV_README.txt
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/20181203_biallelic_SNV_manifest.txt
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr1.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr1.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr2.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr2.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr3.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr3.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr4.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr4.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr5.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr5.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr6.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr6.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr7.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr7.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr8.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr8.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr9.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr9.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr10.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr10.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr11.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr11.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr12.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr12.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr13.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr13.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr14.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr14.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr15.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr15.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr16.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr16.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr17.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr17.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr18.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr18.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr19.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr19.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr20.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr20.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr21.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr21.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chrX.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chrX.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.vcf.gz.tbi
