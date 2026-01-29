#!/bin/bash
#SBATCH --job-name=krakenDB
#SBATCH --output=krakenDB.o%j
#SBATCH --error=krakenDB.e%j
#SBATCH --partition=longq7
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16
#SBATCH --time=8:00:00

# Load modules if needed
module load mamba/latest
source activate get_data
module load kraken2-2.0.9-beta-gcc-12.1.0

# Database path
DB=/scratch/kcarls36/projects/data/krakendb
SYM=/scratch/kcarls36/projects/data/symGenomes

# Download taxonomy and standard libraries
kraken2-build --download-taxonomy --db $DB
kraken2-build --download-library archaea --threads 16 --db $DB
kraken2-build --download-library bacteria --threads 16 --db $DB
kraken2-build --download-library viral --threads 16 --db $DB
kraken2-build --download-library human --threads 16 --db $DB
kraken2-build --download-library fungi --threads 16 --db $DB
kraken2-build --download-library protozoa --threads 16 --db $DB
kraken2-build --download-library UniVec_Core --threads 16 --db $DB

# Add Symbiodiniaceae genomes with taxids
sed '/>/ s/$/|kraken:taxid|2951/' $SYM/Symbiodinium_microadriacticum_genome.scaffold.fasta > $SYM/S_microadriacticum.fa
sed '/>/ s/$/|kraken:taxid|2499525/' $SYM/Breviolum_minutum.v1.0.genome.fa > $SYM/B_minutum.fa
sed '/>/ s/$/|kraken:taxid|2562237/' $SYM/Cladocopium_goreaui_Genome.Scaffolds.fasta > $SYM/C_goreaui.fa
sed '/>/ s/$/|kraken:taxid|1381693/' $SYM/102_symbd_genome_scaffold.fa > $SYM/D_trenchii.fa

# Add to library
kraken2-build --add-to-library $SYM/S_microadriacticum.fa --db $DB
kraken2-build --add-to-library $SYM/B_minutum.fa --db $DB
kraken2-build --add-to-library $SYM/C_goreaui.fa --db $DB
kraken2-build --add-to-library $SYM/D_trenchii.fa --db $DB

# Build database
kraken2-build --build --threads 16 --db $DB

