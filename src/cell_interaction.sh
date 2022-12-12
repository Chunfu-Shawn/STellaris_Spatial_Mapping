#!/bin/sh 

usage(){
  echo "Usage: bash $(basename $0) --sc_coordinate sc_coordinate.csv --species species --outDir output_directory [-h]"
  echo "Author: Kevin Lee"
  echo "Description: This script portrays cell-cell interaction in spatial niche for single cells."
  echo "Date: 2022-11-28"
  echo "------------------------------------------------"
  echo "OPTIONS"
  echo -e "     --sc_coordinate \t\tPath to spatial coordinate of registered single cell (sc_coordinate.csv) [Required]"
  echo -e "     --divergence_cutoff \t\tFilter out divergence lower than cutoff percentile, ranging from 0 to 1. [default: 0.5]"
  echo -e "     --band_width \t\tBandwidths for x and y directions, more details in KernelDensity function in sklearn.neighbors package. [default: 20]"
  echo -e "     --n_bootstrap \t\tNumber of bootstrapping iterations. [default: 20]"
  echo -e "     --species \t\tQuery species, could be 'Mus musculus' or 'Homo sapiens'"
  echo -e "     --n_threads \t\tNumber of threads available to perform MIA [default: 30]"
  echo -e "     --outDir \t\tPath to output directory [Required]"
  echo -e "     -h|--help \t\tprint this help page"
  exit 1
}

## Default argument
divergence_cutoff=0.5
band_width=20
n_bootstrap=20
n_threads=30

while [[ $# -gt 0 ]]; do
    case $1 in
        --sc_coordinate)           sc_coordinate=$2;shift;;
        --divergence_cutoff)       divergence_cutoff=$2;shift;;
        --band_width)              band_width=$2;shift;;
        --n_bootstrap)             n_bootstrap=$2;shift;;
        --species)                 species=$2;shift;;
        --n_threads)               n_threads=$2;shift;;
        --outDir)                  outDir=$2;shift;;
        -h)                        usage;exit 1;;
        --)                        shift; break;;
        *)                         usage; echo -e "\n[ERR] $(date) Unkonwn option: $1"; exit 1;;
   esac
    shift
done

## Check necessary arguments
if [ -z $sc_coordinate ];then
   echo "The coordinate of registered single cells must be provided!" && usage
fi

if [ -z "$species" ];then
   echo "Query species was not provided!" && usage
fi

if [ -z $outDir ];then
   echo "Please specify output directory attributed to the current job!" && usage
fi

## Print argument
echo -e "*** Arguments"
echo -e "sc_coordinate\t$sc_coordinate"
echo -e "divergence_cutoff\t$divergence_cutoff"
echo -e "band_width\t$band_width"
echo -e "n_bootstrap\t$n_bootstrap"
echo -e "species\t$species"
echo -e "n_threads\t$n_threads"
echo -e "outDir\t$outDir"

## Configuration
scriptDir=/home/user/data3/uplee/projects/spatialTransWeb/bin/temp/merged/0123-overall
source /home/user/BGM/uplee/anaconda3/bin/activate spatialWeb
mkdir -p $outDir/log $outDir/out/table $outDir/out/json $outDir/out/pdf 

echo -e "*** Execution"

## Build contact map
echo -e "`date`\tBuild cell-cell contact map..."
time (python $scriptDir/build_contact_map.py \
 --sc_coordinate $sc_coordinate  \
 --divergence_cutoff $divergence_cutoff  \
 --band_width $band_width  \
 --n_bootstrap $n_bootstrap  \
 --outDir $outDir \
 --n_threads $n_threads ) &>$outDir/log/build_contact_map.log

if [ ! $? -eq 0  ];then
    echo -e "`date`\tBuilding cell-cell contact map failed for some reason, please check log files!" >&2
    exit 1
else
    echo -e "`date`\tSuccessfully performed building cell-cell contact map!"
fi

## Run cellphonedb
echo -e "`date`\tDetect ligand-receptor interactions in microenvironment..."

source /home/user/BGM/uplee/anaconda3/bin/activate cellphonedb

meta=$outDir/meta.txt
counts=$outDir/sc_registered.h5ad
micro=$outDir/out/table/microenvironment.csv
output_path_cpdb=$outDir/out/table

time (cellphonedb method statistical_analysis $meta $counts \
        --subsampling \
        --subsampling-log True \
        --counts-data hgnc_symbol \
        --microenvs $micro \
        --output-path $output_path_cpdb \
        --species "$species" \
        --threads $n_threads ) &>$outDir/log/cellphonedb.log

if [ ! $? -eq 0  ];then
    echo -e "`date`\tDetecting cell-cell ligand-receptor interactions failed for some reason, please check log files!" >&2
    exit 1
else
    echo -e "`date`\tSuccessfully performed detecting cell-cell ligand-receptor interactions!"
fi

## Plot cell-cell interactions

time (Rscript $scriptDir/plot_interaction.R -w $outDir) &>$outDir/log/plot_interaction.log

if [ ! $? -eq 0  ];then
    echo -e "`date`\tPlotting cell-cell interactions failed for some reason, please check log files!" >&2
    exit 1
else
    echo -e "`date`\tSuccessfully performed plotting cell-cell interactions!"
fi

echo -e "`date`\tAll finished!"
