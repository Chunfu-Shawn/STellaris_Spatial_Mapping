#!/bin/sh 

usage(){
  echo "Usage: bash $(basename $0) --sc_h5ad sc.h5ad --key_celltype cell_type --dataset ST_dataset_id --section ST_section_id --outDir output_directory [-h]"
  echo "Author: Kevin Lee"
  echo "Description: This script anchor single cells to spatial niche based on ST data, and then portrays cell-cell interaction in spatial niche for single cells."
  echo "Date: 2022-11-28"
  echo "------------------------------------------------"
  echo "OPTIONS"
  echo -e "     --sc_h5ad \t\tPath to h5ad of scRNA-seq (sc.h5ad) [Required]"
  echo -e "     --key_celltype \t\tColumn name of cell type [default: cell_type]"
  echo -e "     --dataset \t\tST dataset ID chosen by user for spatial mapping [Required]"
  echo -e "     --section \t\tST section ID chosen by user for spatial mapping [Required]"
  echo -e "     --divergence_cutoff \t\tFilter out divergence lower than cutoff percentile, ranging from 0 to 1. [default: 0.5]"
  echo -e "     --band_width \t\tBandwidths for x and y directions, more details in KernelDensity function in sklearn.neighbors package. [default: 20]"
  echo -e "     --n_bootstrap \t\tNumber of bootstrapping iterations. [default: 20]"
  echo -e "     --species \t\tQuery species, could be 'Mus musculus' or 'Homo sapiens'"  
  echo -e "     --n_threads \t\tNumber of threads available to perform random forest prediction [default: 30]"
  echo -e "     --outDir \t\tPath to output directory [Required]"
  echo -e "     -h|--help \t\tprint this help page"
  exit 1
}

## Default argument
key_celltype=cell_type
divergence_cutoff=0.5
band_width=20
n_bootstrap=20
n_threads=30

while [[ $# -gt 0 ]]; do
    case $1 in
        --sc_h5ad)            sc_h5ad=$2;shift;;
        --key_celltype)       key_celltype=$2;shift;;
        --dataset)            dataset=$2;shift;;
        --section)            section=$2;shift;;
        --divergence_cutoff)  divergence_cutoff=$2;shift;;
        --band_width)         band_width=$2;shift;;
        --n_bootstrap)        n_bootstrap=$2;shift;;
        --species)            species=$2;shift;;        
        --n_threads)          n_threads=$2;shift;;
        --outDir)             outDir=$2;shift;;
        -h)                   usage;exit 1;;
        --)                   shift; break;;
        *)                    usage; echo -e "\n[ERR] $(date -u) Unkonwn option: $1"; exit 1;;
   esac
    shift
done

## Check necessary arguments
if [ -z $sc_h5ad ];then
   echo "The sc.h5ad must be provided!" && usage
fi

if [ -z $dataset ];then
   echo "No ST datasets were provided!" && usage
fi

if [ -z $section ];then
   echo "No sections of ST datasets were retrieved!" && usage
fi

if [ -z "$species" ];then
   echo "Query species was not provided!" && usage
fi

if [ -z $outDir ];then
   echo "Please specify output directory attributed to the current job!" && usage
fi

## Print argument
#echo -e "*** Arguments"
#echo -e "sc_h5ad\t$sc_h5ad"
#echo -e "key_celltype\t$key_celltype"
#echo -e "dataset\t$dataset"
#echo -e "section\t$section"
#echo -e "divergence_cutoff\t$divergence_cutoff"
#echo -e "band_width\t$band_width"
#echo -e "n_bootstrap\t$n_bootstrap"
#echo -e "species\t$species"
#echo -e "n_threads\t$n_threads"
#echo -e "outDir\t$outDir"

## Configuration
scriptDir=/home/user/data3/uplee/projects/spatialTransWeb/bin/temp/merged/0123-overall
source /home/user/BGM/uplee/anaconda3/bin/activate spatialWeb
mkdir -p $outDir/log $outDir/out/table $outDir/out/json $outDir/out/pdf 

echo -e "*** Execution"

## Run celltrek
echo -e "`date -u`\tRun spatial mapping..."

Rscript $scriptDir/run_celltrek.R \
 -sc_h5ad=$sc_h5ad \
 -key_celltype=$key_celltype \
 -dataset=$dataset \
 -section=$section \
 -n_threads=$n_threads \
 -outDir=$outDir

if [ ! $? -eq 0  ];then
    echo -e "`date -u`\tSpatial mapping failed for some reason, please check log files!" >&2
    exit 1
else
    echo -e "`date -u`\tSuccessfully performed spatial mapping!"
fi

## Convert celltrek results to jsonl for visualization
echo -e "`date -u`\tPrepare visualization..."
 
# Original scRNA-seq
echo -e "`date -u +'%H:%M:%S'` Original cells..."

python $scriptDir/prepare_data.py \
 --dataset $outDir/sc_reduction.h5ad \
 --name sc_reduction \
 --group $key_celltype \
 --outDir $outDir 1>$outDir/log/prepare_data.sc_reduction.log

if [ ! $? -eq 0  ];then
    echo -e "`date -u`\tVisualization preparation of original scRNA-seq failed for some reason, please check log files!" >&2
    exit 1
else
    echo -e "`date -u`\tSuccessfully performed visualization preparation of original scRNA-seq!"
fi

# Registered cells
echo -e "`date -u +'%H:%M:%S'` Registered cells..."

python $scriptDir/prepare_data.sc_registered.py \
 --dataset $dataset \
 --section $section \
 --sc_registered $outDir/sc_registered.h5ad \
 --name sc_registered \
 --outDir $outDir 1>$outDir/log/prepare_data.sc_registered.log

if [ ! $? -eq 0  ];then
    echo -e "`date -u`\tVisualization preparation of registered cells failed for some reason, please check log files!" >&2
    exit 1
else
    echo -e "`date -u`\tSuccessfully performed visualization preparation of registered cells!"
fi

## Build contact map
echo -e "`date -u`\tBuild cell-cell contact map..."

sc_coordinate=$outDir/out/table/sc_coordinate.csv

python $scriptDir/build_contact_map.py \
 --sc_coordinate $sc_coordinate  \
 --divergence_cutoff $divergence_cutoff  \
 --band_width $band_width  \
 --n_bootstrap $n_bootstrap  \
 --n_threads $n_threads \
 --outDir $outDir 1>$outDir/log/build_contact_map.log

if [ ! $? -eq 0  ];then
    echo -e "`date -u`\tBuilding cell-cell contact map failed for some reason, please check log files!" >&2
    exit 1
else
    echo -e "`date -u`\tSuccessfully performed building cell-cell contact map!"
fi

## Run cellphonedb
echo -e "`date -u`\tDetect ligand-receptor interactions in microenvironment..."

source /home/user/BGM/uplee/anaconda3/bin/activate cellphonedb

meta=$outDir/meta.txt
counts=$outDir/sc_registered.h5ad
micro=$outDir/out/table/microenvironment.csv
output_path_cpdb=$outDir/out/table

cellphonedb method statistical_analysis $meta $counts \
        --subsampling \
        --subsampling-log True \
        --counts-data hgnc_symbol \
        --microenvs $micro \
        --species "$species" \
        --threads $n_threads \
        --output-path $output_path_cpdb 1>$outDir/log/cellphonedb.log

if [ ! $? -eq 0  ];then
    echo -e "`date -u`\tDetecting cell-cell ligand-receptor interactions failed for some reason, please check log files!" >&2
    exit 1
else
    echo -e "`date -u`\tSuccessfully performed detecting cell-cell ligand-receptor interactions!"
fi

## Plot cell-cell interactions
echo -e "`date -u`\tPlot cell-cell interactions..."

Rscript $scriptDir/plot_interaction.R -w $outDir 1>$outDir/log/plot_interaction.log

if [ ! $? -eq 0  ];then
    echo -e "`date -u`\tPlotting cell-cell interactions failed for some reason, please check log files!" >&2
    exit 1
else
    echo -e "`date -u`\tSuccessfully performed plotting cell-cell interactions!"
fi

echo -e "`date -u`\tAll finished!"
