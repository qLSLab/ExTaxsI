#!/bin/bash
# change the package name to the existing PyPi package you would like to build
pkg='extaxsi'
# adjust the Python versions you would like to build
array=( 3.5 3.6 3.7 3.8 3.9 )
echo "Building conda package ..."
#conda skeleton pypi $pkg
cd $pkg
conda update conda
conda update conda-build
cd /home/abrusati/Documents/SCRIPTS/TMP
# building conda packages
for i in "${array[@]}"
do
	conda-build --python $i $pkg -c conda-forge -c etetoolkit -c plotly
done
# convert package to other platforms
cd ~
platforms=( osx-64 linux-32 linux-64 win-32 win-64 )
find $HOME/conda-bld/linux-64/ -name *.tar.bz2 | while read file
do
    echo $file
    #conda convert --platform all $file  -o $HOME/conda-bld/
    for platform in "${platforms[@]}"
    do
       conda convert --platform $platform $file  -o $HOME/conda-bld/
    done    
done
# upload packages to conda
find $HOME/conda-bld/ -name *.tar.bz2 | while read file
do
    echo $file
    anaconda upload $file
done
echo "Building conda package done!"
