################################################
#File Name: run.sh
#Author: Up Lee 
#Mail: uplee@pku.edu.cn
#Created Time: Tue 29 Nov 2022 10:10:39 PM CST
################################################

#!/bin/sh

#### 2022-11-29 ####

###########################
#### Script preparation
###########################

cp /home/user/data3/uplee/projects/spatialTransWeb/bin/temp/merged/023-nicheAnchor-cellInteraction/build_contact_map.py ./
cp /home/user/data3/uplee/projects/spatialTransWeb/bin/temp/merged/023-nicheAnchor-cellInteraction/cell_interaction.sh ./
cp /home/user/data3/uplee/projects/spatialTransWeb/bin/temp/merged/023-nicheAnchor-cellInteraction/cellphoneDB_utils.py ./
cp /home/user/data3/uplee/projects/spatialTransWeb/bin/temp/merged/023-nicheAnchor-cellInteraction/celltrek_utils.R ./
cp /home/user/data3/uplee/projects/spatialTransWeb/bin/temp/merged/023-nicheAnchor-cellInteraction/MIA_align.py ./
cp /home/user/data3/uplee/projects/spatialTransWeb/bin/temp/merged/023-nicheAnchor-cellInteraction/MIA_utils.py ./
cp /home/user/data3/uplee/projects/spatialTransWeb/bin/temp/merged/023-nicheAnchor-cellInteraction/nicheAnchor-cellInteraction.sh ./
cp /home/user/data3/uplee/projects/spatialTransWeb/bin/temp/merged/023-nicheAnchor-cellInteraction/nicheAnchor.sh ./
cp /home/user/data3/uplee/projects/spatialTransWeb/bin/temp/merged/023-nicheAnchor-cellInteraction/orthologs_mouse2human.json ./
cp /home/user/data3/uplee/projects/spatialTransWeb/bin/temp/merged/023-nicheAnchor-cellInteraction/plot_interaction.R ./
cp /home/user/data3/uplee/projects/spatialTransWeb/bin/temp/merged/023-nicheAnchor-cellInteraction/prepare_data.py ./
cp /home/user/data3/uplee/projects/spatialTransWeb/bin/temp/merged/023-nicheAnchor-cellInteraction/prepare_data.sc_registered.py ./
cp /home/user/data3/uplee/projects/spatialTransWeb/bin/temp/merged/023-nicheAnchor-cellInteraction/prepare_utils.py ./
cp /home/user/data3/uplee/projects/spatialTransWeb/bin/temp/merged/023-nicheAnchor-cellInteraction/process_sc.py ./
cp /home/user/data3/uplee/projects/spatialTransWeb/bin/temp/merged/023-nicheAnchor-cellInteraction/__pycache__ ./
cp /home/user/data3/uplee/projects/spatialTransWeb/bin/temp/merged/023-nicheAnchor-cellInteraction/run_celltrek.R ./
cp /home/user/data3/uplee/projects/spatialTransWeb/bin/temp/merged/023-nicheAnchor-cellInteraction/ST_screening.sh ./
cp ../023-nicheAnchor-cellInteraction/compress_results.sh  ./

###########################
#### Revise run_celltrek.R & celltrek_utils.R (revision1)
###########################

mkdir -p revision/revision1

#### EO 2022-11-29 ####


#### 2022-11-30 ####

###########################
#### Example (STW-M-Brain-Stereo-seq-1_Mouse-corticogenesis)
###########################

mkdir -p STW-M-Brain-Stereo-seq-1_Mouse-corticogenesis

## ST_screening.sh

mkdir -p STW-M-Brain-Stereo-seq-1_Mouse-corticogenesis/rawdata/
mkdir -p STW-M-Brain-Stereo-seq-1_Mouse-corticogenesis/log

ln -s /home/user/data3/uplee/projects/spatialTransWeb/scRNA-seq/preparation/Mouse-corticogenesis/counts.csv.gz STW-M-Brain-Stereo-seq-1_Mouse-corticogenesis/rawdata/
ln -s /home/user/data3/uplee/projects/spatialTransWeb/scRNA-seq/preparation/Mouse-corticogenesis/labels.csv.gz STW-M-Brain-Stereo-seq-1_Mouse-corticogenesis/rawdata/

count=STW-M-Brain-Stereo-seq-1_Mouse-corticogenesis/rawdata/counts.csv.gz
label=STW-M-Brain-Stereo-seq-1_Mouse-corticogenesis/rawdata/labels.csv.gz
key_celltype=cell_type
datasets="STW-M-Brain-Slide-seq2-1,STW-M-Brain-ST-1,STW-M-Brain-ST-3,STW-M-Brain-Stereo-seq-1,STW-M-Brain-Visium-1,STW-M-Brain-Visium-3,STW-M-Brain-Visium-6,STW-M-Embryo-DBiT-seq-1,STW-M-Embryo-DBiT-seq-2,STW-M-Kidney-Visium-4,STW-M-Liver-ST-1"
sections="GSM5173925_OB1_01,Rep10_MOB,Hip_adpolb_rep1,coronal_1,GSM5621972,ST8059048,WT,GSM4096261_10t,FFPE-1,f12hr_140_processed,CN65-D2"
n_threads=100
outDir=STW-M-Brain-Stereo-seq-1_Mouse-corticogenesis

time (bash ST_screening.sh \
 --count $count  \
 --label $label \
 --key_celltype $key_celltype \
 --datasets $datasets \
 --sections $sections \
 --n_threads $n_threads  \
 --outDir $outDir) &>$outDir/log/ST_screening.log

## nicheAnchor.sh + cell_interaction.sh

sc_h5ad=STW-M-Brain-Stereo-seq-1_Mouse-corticogenesis/sc.h5ad
dataset=STW-M-Brain-Stereo-seq-1
section=coronal_1
divergence_cutoff=0.5
band_width=20
n_bootstrap=20
species='Mus musculus'

time (bash nicheAnchor-cellInteraction.sh \
 --sc_h5ad $sc_h5ad \
 --key_celltype $key_celltype \
 --dataset $dataset \
 --section $section \
 --divergence_cutoff $divergence_cutoff \
 --band_width $band_width \
 --n_bootstrap $n_bootstrap \
 --species "$species" \
 --n_threads $n_threads \
 --outDir $outDir ) &>$outDir/log/nicheAnchor-cellInteraction.log


###########################
#### Add support of tsv/txt (revision2)
###########################

mkdir -p revision/revision2

# See revision/revision2/run.sh

###########################
#### Debug for duplicated cell name
###########################

mkdir -p revision/revision3

# See revision/revision3/run.sh

# Paired-Tag_seq_RNA_filtered_matrix_H3K27ac_Mouse_res_h5ad.zip: from qijt

###########################
#### Debug for cellphonedb error
###########################

mkdir -p revision/revision4

# See revision/revision4/run.sh

#### EO 2022-11-30 ####

#### 2022-12-02 ####
###########################
#### Debug for cellphonedb error
###########################

mkdir -p revision/revision5

# See revision/revision5/run.sh

#### EO 2022-12-02 ####

#### 2022-12-07 ####

###########################
#### Revise log
###########################

mkdir -p revision/revision6

# See revision/revision6/run.sh

###########################
#### Update scrips
###########################

# NOTE: Update scripts to revise log
# (1) Update ST_screening.sh
# (2) Update nicheAnchor-cellInteraction.sh, run_celltrek.R, celltrek_utils.R
# (3) Update cellphoneDB_utils.py ( node_color = np.array(plt.cm.YlOrRd(0.6)).reshape(1,-1) ). No need to update!

#### Backup scripts

mkdir -p backup

cp ST_screening.sh backup/ST_screening.bak20221207.sh
cp nicheAnchor-cellInteraction.sh  backup/nicheAnchor-cellInteraction.bak20221207.sh
cp run_celltrek.R backup/run_celltrek.bak20221207.R
cp celltrek_utils.R backup/celltrek_utils.bak20221207.R

cp revision/revision6/ST_screening.sh  ./ # Change scriptDir to this dir
cp revision/revision6/nicheAnchor-cellInteraction.sh  ./ # Change scriptDir to this dir
cp revision/revision6/run_celltrek.R ./ # Change source path of celltrek_utils.R
cp revision/revision6/celltrek_utils.R ./ # No need to do additional change

##############################
#### Example data (Pancreas_cancer)
##############################

mkdir -p Pancreas_cancer

## ST_screening.sh

mkdir -p Pancreas_cancer/rawdata/
mkdir -p Pancreas_cancer/log

ln -s /home/user/data3/uplee/projects/spatialTransWeb/scRNA-seq/preparation/pancreas_cancer/counts.csv.gz Pancreas_cancer/rawdata/
ln -s /home/user/data3/uplee/projects/spatialTransWeb/scRNA-seq/preparation/pancreas_cancer/labels.csv.gz Pancreas_cancer/rawdata

count=Pancreas_cancer/rawdata/counts.csv.gz
label=Pancreas_cancer/rawdata/labels.csv.gz
key_celltype=cell_type
datasets="STW-H-Pancreas-ST-1,STW-H-Pancreas-ST-1"
sections="GSM3036911_PDAC-A-ST1,GSM3405534_PDAC-B-ST1"
n_threads=100
outDir=Pancreas_cancer

bash ST_screening.sh \
 --count $count  \
 --label $label \
 --key_celltype $key_celltype \
 --datasets $datasets \
 --sections $sections \
 --n_threads $n_threads  \
 --outDir $outDir 1>$outDir/log/ST_screening.log 2>$outDir/log/ST_screening.err

## nicheAnchor.sh + cell_interaction.sh

sc_h5ad=Pancreas_cancer/sc.h5ad
dataset=STW-H-Pancreas-ST-1
section=GSM3036911_PDAC-A-ST1
divergence_cutoff=0.2
band_width=20
n_bootstrap=20
species='Homo sapiens'

bash nicheAnchor-cellInteraction.sh \
 --sc_h5ad $sc_h5ad \
 --key_celltype $key_celltype \
 --dataset $dataset \
 --section $section \
 --divergence_cutoff $divergence_cutoff \
 --band_width $band_width \
 --n_bootstrap $n_bootstrap \
 --species "$species" \
 --n_threads $n_threads \
 --outDir $outDir 1>$outDir/log/nicheAnchor-cellInteraction.log 2>$outDir/log/nicheAnchor-cellInteraction.err

## compress_results.sh

bash compress_results.sh --outDir $outDir 1>$outDir/log/compress_results.log 2>$outDir/log/compress_results.err

#### EO 2022-12-07 ####
