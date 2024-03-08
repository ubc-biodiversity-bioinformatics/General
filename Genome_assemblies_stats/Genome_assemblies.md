
<!-- Genome Assemblies -->
## Genome assemblies

Here are some general notes on how to make a reference genome assembly from PacBio or Nanopore sequence data.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- HELPFUL RESOURCES -->
## Helpful resources

Assembly quality
* Assembly stats (N50): https://github.com/sanger-pathogens/assembly-stats
* BUSCO: 
* kmer estimates: https://kat.readthedocs.io/en/latest/

Base calling from ONT Nanopore fast5 files
* https://github.com/nanoporetech/bonito 

Contig assembly 
* PacBio: Hifiasm https://github.com/chhylp123/hifiasm
* Nanopore: Shasta https://github.com/chanzuckerberg/shasta 

Scaffold assembly and polish
* Visualizing de novo assembly graph: https://rrwick.github.io/Bandage/
* Medaka https://github.com/nanoporetech/medaka 
* Longstitch https://github.com/bcgsc/LongStitch 

HiC Scaffolding
* nf-core-hic https://svn.bcgsc.ca/bitbucket/projects/CAN/repos/hic_scaffolding/browse 
* Salsa https://github.com/marbl/SALSA 

HiC Map editing
* 3D DNA 
* Juicer 
* Juicebox - online version cannot edit https://aidenlab.org/juicebox/ 

Methylation calling 
* Nanopore Remora

<!-- Assembly Quality  -->
## Evaluating Genome Assembly Quality

### General stats
There are a number of statistics that can be used to evaluate a genome assembly. This github package provides a number of them for any input fasta file.

To install:
   ```
cd ~
git clone https://github.com/sanger-pathogens/assembly-stats.git
cd assembly-stats
mkdir build
cd build
cmake -DINSTALL_DIR:PATH=/home/msandler/assembly-stats/ ..
make
make test
make install
   ```
Running once installed
   ```
module load StdEnv/2020 
module load gcc/9.3.0
module load r-bundle-bioconductor/3.14
cd ~/assembly-stats/

./assembly-stats /home/celphin/projects/def-henryg/celphin/Oxyria/Oxyria1_Sept6_haplotigs.p_ctg.fa

stats for /home/celphin/projects/def-henryg/celphin/Oxyria/Oxyria1_Sept6_haplotigs.p_ctg.fa
sum = 590 373 804, n = 974, ave = 606 133.27, largest = 100 326 276
N50 = 56 351 117, n = 4
N60 = 39 295 833, n = 6
N70 = 34 734 330, n = 7
N80 = 28 174 138, n = 9
N90 = 10 901 797, n = 13
N100 = 14349, n = 974
N_count = 0
Gaps = 0

   ```

### BUSCO

BUSCO is a software that finds conserved genes in your assembly and determines based on the number of genes expected how complete the assembly is.
BUSCO is also used to evaluate genome annotations to determine how many genes are missing from an annotation compared to the assembly.

Installation:
   ```
cd ~

# Setup
#module purge
module load StdEnv/2020
module load gcc/9.3.0 python augustus hmmer blast+ metaeuk prodigal r bbmap

virtualenv ~/busco_env

source ~/busco_env/bin/activate

pip install biopython pandas busco==5.1.2 --no-index

# to close
deactivate

# to activate the environment
module load StdEnv/2020
module load gcc/9.3.0 python augustus hmmer blast+ metaeuk prodigal r bbmap
source ~/busco_env/bin/activate

# check installed bbmap
bbduk.sh --version

java -ea -Xmx3235m -Xms3235m -cp /cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/bbmap/38.86/current/ jgi.BBDuk --version
Picked up JAVA_TOOL_OPTIONS: -Xmx2g
BBMap version 38.86
For help, please run the shellscript with no parameters, or look in /docs/.

cd ~
   ```

There are a number of BUSCO libraries. Here we use the eudicots_odb10 but any library can be downloaded from here: https://busco-data.ezlab.org/v5/data/

The data for each library must follow the folder structure of the website above.

Library downloads:
   ```
mkdir BUSCO_downloads; cd BUSCO_downloads

mkdir lineages; cd lineages
wget https://busco-data.ezlab.org/v5/data/lineages/viridiplantae_odb10.2020-09-10.tar.gz
wget https://busco-data.ezlab.org/v5/data/lineages/eudicots_odb10.2020-09-10.tar.gz
wget https://busco-data.ezlab.org/v5/data/lineages/embryophyta_odb10.2020-09-10.tar.gz
cd ..

mkdir information; cd information
wget https://busco-data.ezlab.org/v5/data/information/lineages_list.2021-12-14.txt.tar.gz
cd ..
mkdir placement_files; cd placement_files
wget https://busco-data.ezlab.org/v5/data/placement_files/list_of_reference_markers.eukaryota_odb10.2019-12-16.txt.tar.gz
wget https://busco-data.ezlab.org/v5/data/placement_files/mapping_taxid-lineage.eukaryota_odb10.2019-12-16.txt.tar.gz
wget https://busco-data.ezlab.org/v5/data/placement_files/mapping_taxids-busco_dataset_name.eukaryota_odb10.2019-12-16.txt.tar.gz
wget https://busco-data.ezlab.org/v5/data/placement_files/supermatrix.aln.eukaryota_odb10.2019-12-16.faa.tar.gz
wget https://busco-data.ezlab.org/v5/data/placement_files/tree.eukaryota_odb10.2019-12-16.nwk.tar.gz
wget https://busco-data.ezlab.org/v5/data/placement_files/tree_metadata.eukaryota_odb10.2019-12-16.txt.tar.gz
cd..

cd lineages
tar -xvf ./viridiplantae_odb10.2020-09-10.tar.gz
tar -xvf ./eudicots_odb10.2020-09-10.tar.gz
tar -xvf ./embryophyta_odb10.2020-09-10.tar.gz
cd ..
cd information
tar -xvf ./lineages_list.2021-12-14.txt.tar.gz
cd ..
cd placement_files
tar -xvf ./list_of_reference_markers.eukaryota_odb10.2019-12-16.txt.tar.gz
tar -xvf ./mapping_taxid-lineage.eukaryota_odb10.2019-12-16.txt.tar.gz
tar -xvf ./mapping_taxids-busco_dataset_name.eukaryota_odb10.2019-12-16.txt.tar.gz
tar -xvf ./supermatrix.aln.eukaryota_odb10.2019-12-16.faa.tar.gz
tar -xvf ./tree.eukaryota_odb10.2019-12-16.nwk.tar.gz
tar -xvf ./tree_metadata.eukaryota_odb10.2019-12-16.txt.tar.gz

cd ..
   ```

Running BUSCO
   ```
salloc -c10 --time 2:55:00 --mem 120000m --account def-rieseber
module load StdEnv/2020
module load gcc/9.3.0 python augustus hmmer blast+ metaeuk prodigal r bbmap
source ~/busco_env/bin/activate
#conda deactivate

#Oxyria digyna (our version):
busco --offline --in  Oxyria_digyna.fasta  \
--out  BUSCO_Oxyria_Assembly_eudicots  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/

   ```

 --------------------------------------------------
        |Results from dataset eudicots_odb10              |
        --------------------------------------------------
        |C:93.4%[S:85.2%,D:8.2%],F:1.8%,M:4.8%,n:2326     |
        |2171   Complete BUSCOs (C)                       |
        |1981   Complete and single-copy BUSCOs (S)       |
        |190    Complete and duplicated BUSCOs (D)        |
        |43     Fragmented BUSCOs (F)                     |
        |112    Missing BUSCOs (M)                        |
        |2326   Total BUSCO groups searched               |
        --------------------------------------------------
		
		
<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- Contig assembly-->
## Contig assembly

https://github.com/nanoporetech/bonito 

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- Contig assembly-->
## Contig assembly

Hifiasm
   ```
# need ccs reads
tmux new-session -s bam
tmux attach-session -t bam

cd /home/celphin/projects/def-henryg/celphin/Oxyria/PacBio/PacBio_Aug2022

module load StdEnv/2020
module load samtools/1.15.1

samtools bam2fq Sequel.RunS218_S2.003.Oxyria_C.ccs.bam > Oxyria_0003_ccs.fastq 
[M::bam2fq_mainloop] discarded 0 singletons
[M::bam2fq_mainloop] processed 1138973 reads

samtools bam2fq Sequel.RunS218_S2.004.Oxyria_C.ccs.bam > Oxyria_0004_ccs.fastq 
[M::bam2fq_mainloop] discarded 0 singletons
[M::bam2fq_mainloop] processed 1237599 reads

# zip up data
salloc -c1 --time 2:50:00 --mem 1510G --account def-henryg
tar -cvzf HiFi-ccs_reads.tar.gz Oxyria_0003_ccs.fastq Oxyria_0004_ccs.fastq 

# try again with ccs reads
salloc -c32 --time 23:50:00 --mem 1510G --account def-henryg

/home/celphin/projects/def-henryg/celphin/Oxyria/hifiasm/hifiasm
 -o /home/celphin/projects/def-henryg/celphin/Oxyria/Oxyria1_Sept6.asm
 -t32 --h1 /home/celphin/projects/def-henryg/celphin/Oxyria/Hi-C/NS.1871.003.IDT_i7_217---IDT_i5_217.Oxyria_Hi-C_R1.fastq.gz
 --h2 /home/celphin/projects/def-henryg/celphin/Oxyria/Hi-C/NS.1871.003.IDT_i7_217---IDT_i5_217.Oxyria_Hi-C_R2.fastq.gz
 /home/celphin/projects/def-henryg/celphin/Oxyria/PacBio/PacBio_Aug2022/HiFi-ccs_reads.tar.gz
   ```

Shasta
   ```

   ```

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- Scaffold assembly -->
## Polishing your assembly

Medaka for polishing your assembly
   ```
# August 2022
# Medaka - Contig polish
# https://github.com/nanoporetech/medaka

tmux new-session -s Plect
tmux attach-session -t Plect

cd /home/celphin/projects/def-cronk/celphin/Plectritis_analysis/Assembly/Polish/Medaka

virtualenv medaka --python=python3 --prompt "(medaka) "
source /home/celphin/projects/def-cronk/celphin/Plectritis_analysis/Assembly/Polish/Medaka/medaka/bin/activate
pip install medaka

#-------------------------------
source /home/celphin/projects/def-cronk/celphin/Plectritis_analysis/Assembly/Polish/Medaka/medaka/bin/activate
module load StdEnv/2020
module load python/3.10.2
module load scipy-stack/2022a
module load gcc/9.3.0
module load cuda/11.7
module load hdf5/1.12.1
module load parasail/2.5
module load arrow/8.0.0
module load samtools/1.15.1
module load minimap2/2.17
module load tabix/0.2.6
module load htslib/1.15.1

medaka -h
medaka tools list\_models
-m r941_min_high_g303


Available: r103_fast_g507, r103_fast_snp_g507, r103_fast_variant_g507, r103_hac_g507, r103_hac_snp_g507, r103_hac_variant_g507, r103_min_high_g345, r103_min_high_g360, 
r103_prom_high_g360, r103_prom_snp_g3210, r103_prom_variant_g3210, r103_sup_g507, r103_sup_snp_g507, r103_sup_variant_g507, r104_e81_fast_g5015, r104_e81_hac_g5015,
 r104_e81_sup_g5015, r10_min_high_g303, r10_min_high_g340, r941_min_fast_g303, r941_min_fast_g507, r941_min_fast_snp_g507, r941_min_fast_variant_g507, r941_min_hac_g507,
 r941_min_hac_snp_g507, r941_min_hac_variant_g507, r941_min_high_g303, r941_min_high_g330, r941_min_high_g340_rle, r941_min_high_g344, r941_min_high_g351, 
r941_min_high_g360, r941_min_sup_g507, r941_min_sup_snp_g507, r941_min_sup_variant_g507, r941_prom_fast_g303, r941_prom_fast_g507, r941_prom_fast_snp_g507, 
r941_prom_fast_variant_g507, r941_prom_hac_g507, r941_prom_hac_snp_g507, r941_prom_hac_variant_g507, r941_prom_high_g303, r941_prom_high_g330, r941_prom_high_g344, 
r941_prom_high_g360, r941_prom_high_g4011, r941_prom_snp_g303, r941_prom_snp_g322, r941_prom_snp_g360, r941_prom_sup_g507, r941_prom_sup_snp_g507,
 r941_prom_sup_variant_g507, r941_prom_variant_g303, r941_prom_variant_g322, r941_prom_variant_g360
Default consensus:  r941_min_hac_g507
Default variant:  r941_min_hac_variant_g507

medaka_consensus 
-i /home/celphin/projects/def-cronk/celphin/Plectritis_analysis/Assembly/fastq/basecalls_Apr*.fastq  
-d /home/celphin/projects/def-cronk/celphin/Plectritis_analysis/Assembly/fasta/Assembly_May14.fasta 
-o /home/celphin/projects/def-cronk/celphin/Plectritis_analysis/Assembly/Polish/Medaka/medaka_consensus_aug14
-t 30 
-m r941_min_sup_g507


#-------------------------------
salloc -c30 --time 23:50:00 --mem 120000m --account def-cronk

source /home/celphin/projects/def-cronk/celphin/Plectritis_analysis/Assembly/Polish/Medaka/medaka/bin/activate
module load StdEnv/2020
module load python/3.10.2
module load scipy-stack/2022a
module load gcc/9.3.0
module load cuda/11.7
module load hdf5/1.12.1
module load parasail/2.5
module load arrow/8.0.0
module load samtools/1.15.1
module load minimap2/2.17
module load tabix/0.2.6
module load htslib/1.15.1
module load bcftools/1.11


medaka_consensus -i /home/celphin/projects/def-cronk/celphin/Plectritis_analysis/Assembly/fastq/basecalls_Apr10.fastq -d /home/celphin/projects/def-cronk/celphin/Plectritis_analysis/Assembly/fasta/Assembly_May14.fasta -o /home/celphin/projects/def-cronk/celphin/Plectritis_analysis/Assembly/Polish/Medaka/medaka_consensus_aug14 -t 30 -m r941_min_sup_g507

Checking program versions
This is medaka 1.6.0
...
Polished assembly written to /home/celphin/projects/def-cronk/celphin/Plectritis_analysis/Assembly/Polish/Medaka/medaka_consensus_aug14/consensus.fasta, have a nice day.

   ```
   
Longstitch for polishing long read data

   ```
# Using singularity
salloc --ntasks-per-node=2 --mem=8000M --time 00:50:00 --account def-cronk

module load singularity
singularity build longstich-1.0.1.sif docker://chrishah/longstitch:v1.0.1

singularity shell longstich-1.0.1.sif

longstitch --help

#-------------------
# Running Longstitch

salloc --ntasks-per-node=30 --mem=12000M --time 23:00:00 --account def-cronk

module load singularity

singularity exec -B $PWD:/home/celphin \
longstich-1.0.1.sif \
longstitch run draft=consensus reads=basecalls_Apr10 G=14e8 k_ntLink=24

<!-- tigmint-make tigmint-long draft=consensus reads=basecalls_Apr10 cut=250 t=8 G=14e8 span=auto dist=auto
make[1]: Entering directory '/home/celphin'
sh -c 'gunzip -c basecalls_Apr10.fq.gz | \
/opt/conda/bin/share/tigmint-1.2.5-1/bin/tigmint_estimate_dist.py - -n 1000000 -o basecalls_Apr10.tigmint-long.params.tsv'
/opt/conda/bin/share/tigmint-1.2.5-1/bin/../src/long-to-linked-pe -l 250 -m2000 -g14e8 -s -b basecalls_Apr10.barcode-multiplicity.tsv --bx -t8 --fasta -f basecalls_Apr10.tigmint-long.params.tsv basecalls_Apr10.fq.gz | \
minimap2 -y -t8 -x map-ont --secondary=no consensus.fa - | \
/opt/conda/bin/share/tigmint-1.2.5-1/bin/tigmint_molecule_paf.py -q0 -s2000 -p basecalls_Apr10.tigmint-long.params.tsv - | sort -k1,1 -k2,2n -k3,3n > consensus.basecalls_Apr10.cut250.molecule.size2000.distauto.bed
long-to-linked-pe v1.2.5: Using more than 6 threads does not scale, reverting to 6.
... -->

   ```
   
   
<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- HiC Scaffolding -->
## HiC Scaffolding with salsa

Pipeline here: https://svn.bcgsc.ca/bitbucket/scm/can/hic_scaffolding

First need to map HiC read to the genome (HiC Mapping)

Second use these HiC reads for scaffolding (Salsa)

   ```
   #Notes still to come
   ```

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- HiC Map editing -->
## HiC Map editing

https://github.com/aidenlab/Juicebox/wiki
https://github.com/aidenlab/juicer/wiki
https://aidenlab.org/juicebox/

   ```
# Juicer
# https://github.com/aidenlab/juicer/wiki
#-----------------------------
# Questions about code on CC please contact: cassandra.elphinstone@shaw.ca
#-----------------------
# Fasta file: /home/celphin/projects/def-henryg/celphin/Oxyria/Oxy_dig_1.scaffolds_FINAL.fasta
# HiC files: /home/celphin/projects/def-henryg/celphin/Oxyria/Hi-C_raw_data/NS.1871.003.IDT_i7_217---IDT_i5_217.Oxyria_Hi-C_R*.fastq.gz

#-------------------
# Dependancies
    # For alignment and creation of the Hi-C pairs file merged_nodups.txt:
        # GNU CoreUtils
        # Burrows-Wheeler Aligner (BWA)
    # For .hic file creation and Juicer tools analysis:
        # Java 1.7 or 1.8 JDK. (Alternative link for Ubuntu/LinuxMint). Minimum system requirements for running Java can be found at http://java.com/en/download/help/sysreq.xml
        # Latest Juicer Tools jar
    # For peak calling:
        # CUDA and an NVIDIA GPU
        # The native libraries included with Juicer are compiled for CUDA 7. Other versions of CUDA can be used, but you will need to download the respective native libraries from JCuda.
        # For best performance, use a dedicated GPU. You may also be able to obtain access to GPU clusters through Amazon Web Services or a local research institution.

#You should have a Juicer directory containing scripts/, references/, and optionally restriction_sites/, and a different working directory containing fastq/. 
#You should download the Juicer tools jar and install it in your scripts/ directory.

#---------------------------
# Dependancies
module load StdEnv/2020
module load nixpkgs/16.09 
module load bwa/0.7.15
module load java/1.8.0_192 
module load gcc/4.8.5
module load cuda/7.5.18

#------------------------
# Setup
# https://github.com/theaidenlab/juicer/wiki/Running-Juicer-on-a-cluster
cd /home/celphin/projects/def-henryg/celphin/Oxyria/
mkdir Juicer
cd Juicer
mkdir references/

# copy over reference fasta
cp -v /home/celphin/projects/def-henryg/celphin/Oxyria/Oxy_dig_1.scaffolds_FINAL.fasta ./references/

# Get Juicer
git clone https://github.com/theaidenlab/juicer.git

# symbolic link to scripts
ln -s juicer/SLURM/scripts/ scripts
cd scripts
# Download Juicer tools jar
wget https://github.com/aidenlab/Juicebox/releases/download/v2.18.00/juicer_tools.2.18.00.jar
ln -s juicer_tools.2.18.00.jar juicer_tools.jar

cd ..
mkdir fastq; cd fastq
# copy over the HiC files into fastq
cp -v /home/celphin/projects/def-henryg/celphin/Oxyria/Hi-C_raw_data/NS.1871.003.IDT_i7_217---IDT_i5_217.Oxyria_Hi-C_R*.fastq.gz .
cd ../..

#-----------------------------
# copy original scripts and replace with scripts for Cedar SLURM specifically

mv /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/scripts/jucier.sh /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/scripts/jucier-copy.sh
mv /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/scripts/generate_site_positions.py /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/scripts/generate_site_positions-copy.py
mv /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/scripts/juicer_arrowhead.sh /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/scripts/juicer_arrowhead-copy.sh
mv /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/scripts/juicer_postprocessing.sh /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/scripts/juicer_postprocessing-copy.sh
mv /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/scripts/mega.sh /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/scripts/mega-copy.sh

# replace jucier scripts with new CC specific juicer scripts 
git clone https://github.com/rieseberglab/Juicer_Compute_Canada_scripts.git

#----------------------------
# set default SLURM submission account 
# https://docs.alliancecan.ca/wiki/Running_jobs#Accounts_and_projects
export SLURM_ACCOUNT=def-rieseber
export SBATCH_ACCOUNT=$SLURM_ACCOUNT
export SALLOC_ACCOUNT=$SLURM_ACCOUNT

#################################
# Add genome to Juicer
# https://github.com/aidenlab/juicer/wiki/Usage

tmux new-session -s Juicer # Ctrl+B D to dettach
tmux attach-session -t Juicer

# get allocation to start
salloc -c1 --time 2:50:00 --mem 120000m --account def-rieseber

cd /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer

module load StdEnv/2020
module load nixpkgs/16.09 
module load bwa/0.7.15
module load java/1.8.0_192 
module load gcc/4.8.5
module load cuda/7.5.18

# index fasta 
bwa index ./references/Oxy_dig_1.scaffolds_FINAL.fasta
# [main] Real time: 725.084 sec; CPU: 711.967 sec

#------------------------
# add restriction sites to reference

# copy over python script
cp -v /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/juicer/miscgenerate_site_positions.py /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/scripts

module load StdEnv/2020
module load python
python ./scripts/generate_site_positions.py DpnII Oxyria '/project/6064374/celphin/Oxyria/Juicer/references/Oxy_dig_1.scaffolds_FINAL.fasta'
mkdir restriction_sites
mv Oxyria_DpnII.txt restriction_sites

#---------------------
# get chromosome sizes list
awk 'BEGIN{OFS="\t"}{print $1, $NF}' ./restriction_sites/Oxyria_DpnII.txt > ./restriction_sites/Oxyria.chrom.sizes

##################################
# try running

module load StdEnv/2020
module load nixpkgs/16.09 
module load bwa/0.7.15
module load java/1.8.0_192 
module load gcc/4.8.5
module load cuda/7.5.18

cd /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer

./scripts/juicer.sh \
-D /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer \
-z ./references/Oxy_dig_1.scaffolds_FINAL.fasta \
-p ./restriction_sites/Oxyria.chrom.sizes \
-y ./restriction_sites/Oxyria_DpnII.txt


(-: Looking for fastq files...fastq files exist
(-: Aligning files matching /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/fastq/*_R*.fastq*
 in queue default to genome ./references/Oxy_dig_1.scaffolds_FINAL.fasta with no fragment delimited maps.
(-: Created /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/splits and /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/aligned.
(-: Starting job to launch other jobs once splitting is complete
(-: Finished adding all jobs... Now is a good time to get that cup of coffee... Last job id 60897757

####################
# make merged_nodups.txt from bam 

tmux attach-session -t Juicer

salloc -c10 --time 2:50:00 --mem 120000m --account def-rieseber

module load StdEnv/2020 samtools/1.16.1

samtools view -@ 8 -O SAM -F 1024 /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/aligned/merged_dedup.*am | \
awk -v mnd=1 -f /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/scripts/sam_to_pre.awk > \
/home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/aligned/merged_nodups.txt

###################################
# need to run through 3D DNA to get movable scaffolds .assembly file
# https://aidenlab.org/assembly/manual_180322.pdf 

# downloading 3D DNA
# https://github.com/aidenlab/3d-dna

cd /home/celphin/projects/def-henryg/celphin/Oxyria/3D_DNA
git clone https://github.com/aidenlab/3d-dna.git

chmod -R 777 /home/celphin/projects/def-henryg/celphin/Oxyria/3D_DNA/3d-dna/*

#--------------------------
# Prerequisites

    # LastZ (version 1.03.73 released 20150708) – for diploid mode only
    # Java version >=1.7
    # Bash >=4
    # GNU Awk >=4.0.2
    # GNU coreutils sort >=8.11
    # Python >=2.7 - for chromosome number-aware splitter module only
    # scipy numpy matplotlib - for chromosome number-aware splitter module only

module load StdEnv/2020 lastz/1.04.03 java/1.8.0_192 python/3.10.2 scipy-stack/2022a

#-----------------------
# get .assembly file from  merged_nodups.txt and fasta
tmux attach-session -t Juicer

cd /home/celphin/projects/def-henryg/celphin/Oxyria/3D_DNA
awk -f ./3d-dna/utils/generate-assembly-file-from-fasta.awk /home/celphin/projects/def-henryg/celphin/Oxyria/Oxy_dig_1.scaffolds_FINAL.fasta > Oxyria_draft.assembly

#-------------
# run 3D DNA

salloc -c10 --time 2:50:00 --mem 120000m --account def-rieseber
module load StdEnv/2020 lastz/1.04.03 java/1.8.0_192 python/3.10.2 scipy-stack/2022a
cd /home/celphin/projects/def-henryg/celphin/Oxyria/3D_DNA
/home/celphin/projects/def-henryg/celphin/Oxyria/3D_DNA/3d-dna/run-asm-pipeline-post-review.sh \
-r /home/celphin/projects/def-henryg/celphin/Oxyria/3D_DNA/Oxyria_draft.assembly \
/home/celphin/projects/def-henryg/celphin/Oxyria/Oxy_dig_1.scaffolds_FINAL.fasta \
/home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/aligned/merged_nodups.txt

#########################
# open in Juicebox 
# https://aidenlab.org/juicebox/

# Oxy_dig_1.scaffolds_FINAL.final.assembly
# Oxy_dig_1.scaffolds_FINAL.final.hic

# use desktop version to edit manually
# https://www.youtube.com/watch?v=Nj7RhQZHM18
# need .assembly and .hic files

#################################
# after editing 
# convert back to fasta: juicebox_assembly_converter.py
# copy review.assembly into Cedar
# https://github.com/phasegenomics/juicebox_scripts

git clone https://github.com/phasegenomics/juicebox_scripts.git
chmod -R 770 /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/juicebox_scripts/*

tmux new-session -s Juicer 
tmux attach-session -t Juicer

salloc -c1 --time 2:50:00 --mem 120000m --account def-rieseber
module load StdEnv/2020 lastz/1.04.03 java/1.8.0_192 python/3.10.2 scipy-stack/2022a

cd /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/juicebox_scripts/

python /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/juicebox_scripts/juicebox_scripts/juicebox_assembly_converter.py \
-a /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/juicebox_scripts/data/Oxy_dig_1.scaffolds_FINAL.final.review.assembly  \
-f /home/celphin/projects/def-henryg/celphin/Oxyria/Oxy_dig_1.scaffolds_FINAL.fasta 


########################

   ```

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- RagTag genome assembly -->
## RagTag genome assembly

If you want to map your not yet chromosome scale scaffolds to a chromosome scale assembly you can do this using RagTag.
This program will maintain the synteny within your scaffolds but join them to match the species you are mapping to. 

   ```
# RagTag - with order and join scaffolds based on reference sequence - will not change the scaffolds though
# https://github.com/malonge/RagTag
# https://github.com/malonge/RagTag/wiki/scaffold 

# Needs
    # Minimap2
    # Python 3 (with the following auto-installed packages)
        # numpy
        # intervaltree
        # pysam
        # networkx

tmux new-session -s Tox
tmux attach-session -t Tox

cd /home/celphin/scratch/Annotation/Poison_Oak/

python3 -m venv ~/scratch/Annotation/Poison_Oak/RagTag/RagTag_vir_env
source ~/scratch/Annotation/Poison_Oak/RagTag/RagTag_vir_env/bin/activate
# deactivate # to close

module load StdEnv/2020
module load python/3.10.2
module load minimap2/2.24
module load scipy-stack/2022a

pip install RagTag

pip list --local

# test
ragtag.py --version
#v2.1.0

#--------------------------------
cd ~/scratch/Annotation/Poison_Oak/RagTag/
wget https://github.com/malonge/RagTag/archive/refs/heads/master.zip
unzip master.zip
chmod -R 755 .

#--------------------------------
salloc -c32 --time 23:55:00 --mem 120000m --account def-rieseber

cd /home/celphin/scratch/Annotation/Poison_Oak/
source ~/scratch/Annotation/Poison_Oak/RagTag/RagTag_vir_env/bin/activate

module load StdEnv/2020
module load python/3.10.2
module load minimap2/2.24
module load scipy-stack/2022a

pip install RagTag

/home/celphin/scratch/Annotation/Poison_Oak/RagTag/RagTag-master/ragtag.py \
scaffold Mango.fasta T-diversilobum.fasta \
-t 28 -o ./T-diversilobum_ragtag_output/

/home/celphin/scratch/Annotation/Poison_Oak/RagTag/RagTag-master/ragtag.py \
scaffold Mango.fasta T-radicans.fasta \
-t 28 -o ./T-radicans_ragtag_output/

cd 
more ./T-radicans_ragtag_output/ragtag.scaffold.stats
placed_sequences        placed_bp       unplaced_sequences      unplaced_bp       gap_bp     gap_sequences
212                     371 934 398        844                    23 180 968        15 100   151


more ./T-diversilobum_ragtag_output/ragtag.scaffold.stats
placed_sequences        placed_bp       unplaced_sequences         unplaced_bp       gap_bp    gap_sequences
634                     361 529 572       3065                     90 977 449        51700      517


#-------------------------------
# check stats - same as above again
# https://github.com/sanger-pathogens/assembly-stats
#-----------------------
module load StdEnv/2020 
module load gcc/9.3.0
module load r-bundle-bioconductor/3.14
cd ~/assembly-stats/

./assembly-stats /home/celphin/scratch/Annotation/Poison_Oak/T-radicans_ragtag_output/ragtag.scaffold.fasta

stats for /home/celphin/scratch/Annotation/Poison_Oak/T-radicans_ragtag_output/ragtag.scaffold.fasta
sum = 395130466, n = 905, ave = 436608.25, largest = 46349665
N50 = 26424856, n = 6
N60 = 25406786, n = 7
N70 = 20827032, n = 9
N80 = 13093294, n = 11
N90 = 2357007, n = 17
N100 = 23, n = 905
N_count = 15100
Gaps = 151

   ```
   
   
<!-- DNA methylation   -->
## DNA methylation calling 

Nanopore fast5 files contain DNA methylation information. You can use a pre-trained machine learning model to call the methylated cytosines from your raw fast5 reads.

### Remora
   ```


   ```


<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- CONTACT -->
## Contact

If you have questions please contact:

Cassandra Elphinstone - cassandra.elphinstone@gmail.com


<p align="right">(<a href="#readme-top">back to top</a>)</p>

