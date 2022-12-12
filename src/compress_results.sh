#!/bin/sh 

usage(){
  echo "Usage: bash $(basename $0) --outDir output_directory [-h]"
  echo "Author: Kevin Lee"
  echo "Description: This script compresses spatial mapping and cell-cell interaction results."
  echo "Date: 2022-11-29"
  echo "------------------------------------------------"
  echo "OPTIONS"
  echo -e "     --outDir \t\tPath to output directory [Required]"
  echo -e "     -h|--help \t\tprint this help page"
  exit 1
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --outDir)             outDir=$2;shift;;
        -h)                   usage;exit 1;;
        --)                   shift; break;;
        *)                    usage; echo -e "\n[ERR] $(date -u) Unkonwn option: $1"; exit 1;;
   esac
    shift
done

## Check necessary arguments

if [ -z $outDir ];then
   echo "Please specify output directory attributed to the current job!" && usage
fi

## Print argument
#echo -e "*** Arguments"
#echo -e "outDir\t$outDir"

echo -e "*** Execution"

## Compressing
echo -e "`date -u`\tCompressing results..."

# pdf
tar -czf $outDir/figures.tar.gz  $outDir/out/pdf/*.pdf & 

# table
tar -czf $outDir/results.tar.gz  $outDir/out/table/* & 

# pdf + table + sc_registered.h5ad
tar -czf $outDir/all_results.tar.gz $outDir/out/pdf/*pdf $outDir/out/table/* $outDir/sc_registered.h5ad &

wait;wait;wait

echo -e "`date -u`\tAll finished!"
