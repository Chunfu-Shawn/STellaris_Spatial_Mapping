#!/bin/sh 

usage(){
  echo "Usage: bash $(basename $0) --sc_h5ad sc.h5ad --key_celltype cell_type --dataset ST_dataset_id --section ST_section_id --outDir output_directory [-h]"
  echo "Author: Kevin Lee"
  echo "Description: This script anchor single cells to spatial niche based on ST data."
  echo "Date: 2022-11-23"
  echo "------------------------------------------------"
  echo "OPTIONS"
  echo -e "     --sc_h5ad \t\tPath to h5ad of scRNA-seq (sc.h5ad) [Required]"
  echo -e "     --key_celltype \t\tColumn name of cell type [default: cell_type]"
  echo -e "     --dataset \t\tST dataset ID chosen by user for spatial mapping [Required]"
  echo -e "     --section \t\tST section ID chosen by user for spatial mapping [Required]"
  echo -e "     --n_threads \t\tNumber of threads available to perform random forest prediction [default: 30]"
  echo -e "     --outDir \t\tPath to output directory [Required]"
  echo -e "     -h|--help \t\tprint this help page"
  exit 1
}

## Default argument
key_celltype=cell_type
n_threads=30

while [[ $# -gt 0 ]]; do
    case $1 in
        --sc_h5ad)            sc_h5ad=$2;shift;;
        --key_celltype)       key_celltype=$2;shift;;
        --dataset)            dataset=$2;shift;;
        --section)            section=$2;shift;;
        --n_threads)          n_threads=$2;shift;;
        --outDir)             outDir=$2;shift;;
        -h)                   usage;exit 1;;
        --)                   shift; break;;
        *)                    usage; echo -e "\n[ERR] $(date) Unkonwn option: $1"; exit 1;;
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

if [ -z $outDir ];then
   echo "Please specify output directory attributed to the current job!" && usage
fi

## Print argument
echo -e "*** Arguments"
echo -e "sc_h5ad\t$sc_h5ad"
echo -e "key_celltype\t$key_celltype"
echo -e "dataset\t$dataset"
echo -e "section\t$section"
echo -e "n_threads\t$n_threads"
echo -e "outDir\t$outDir"

## Configuration
scriptDir=/home/user/data3/uplee/projects/spatialTransWeb/bin/temp/merged/0123-overall
source /home/user/BGM/uplee/anaconda3/bin/activate spatialWeb
mkdir -p $outDir/log $outDir/out/table $outDir/out/json $outDir/out/pdf 

echo -e "*** Execution"

## Run celltrek
echo -e "`date`\tRun spatial mapping..."
time (Rscript $scriptDir/run_celltrek.R \
 -sc_h5ad=$sc_h5ad \
 -key_celltype=$key_celltype \
 -dataset=$dataset \
 -section=$section \
 -n_threads=$n_threads \
 -outDir=$outDir ) &>$outDir/log/run_celltrek.log

if [ ! $? -eq 0  ];then
    echo -e "`date`\tSpatial mapping failed for some reason, please check log files!" >&2
    exit 1
else
    echo -e "`date`\tSuccessfully performed spatial mapping!"
fi


## Convert celltrek results to jsonl for visualization
echo -e "`date`\tPrepare visualization..."
 
# Original scRNA-seq
echo -e "Original cells..."
time (python $scriptDir/prepare_data.py \
 --dataset $outDir/sc_reduction.h5ad \
 --name sc_reduction \
 --group $key_celltype \
 --outDir $outDir) &>$outDir/log/prepare_data.sc_reduction.log

if [ ! $? -eq 0  ];then
    echo -e "`date`\tVisualization preparation of original scRNA-seq failed for some reason, please check log files!" >&2
    exit 1
else
    echo -e "`date`\tSuccessfully performed visualization preparation of original scRNA-seq!"
fi


# Registered cells
echo -e "Registered cells..."
time (python $scriptDir/prepare_data.sc_registered.py \
 --dataset $dataset \
 --section $section \
 --sc_registered $outDir/sc_registered.h5ad \
 --name sc_registered \
 --outDir $outDir) &>$outDir/log/prepare_data.sc_registered.log

if [ ! $? -eq 0  ];then
    echo -e "`date`\tVisualization preparation of registered cells failed for some reason, please check log files!" >&2
    exit 1
else
    echo -e "`date`\tSuccessfully performed visualization preparation of registered cells!"
fi

echo -e "`date`\tAll finished!"
