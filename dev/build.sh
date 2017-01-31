#!/bin/bash

Rscript -e "remove.packages('poisbinom')"

echo "--------------------------"
echo "Script to compile 'poisbinom' and run test code"
echo

if [ -e "poisbinom" ]&& [ -d "poisbinom" ]; then
    rm -rf poisbinom
    echo
    echo "Old 'poisbinom' folder is removed."
    echo "--------------------------"
    echo
fi


Rscript -e "library(Rcpp); Rcpp.package.skeleton('poisbinom')"

echo
echo
echo "Skeleton package is created."
echo "--------------------------"
echo
echo

cp -r src/* poisbinom/src/
cp -r man/* poisbinom/man/
rm poisbinom/src/rcpp_hello_world.cpp

echo
echo
echo "C++ code and manual are copied."
echo "--------------------------"
echo
echo

cd poisbinom
Rscript -e "library(Rcpp); compileAttributes(verbose=TRUE)"
rm man/poisbinom-package.Rd
rm man/rcpp_hello_world.Rd

echo
echo
echo "'poisbinom' is compiled."
echo "--------------------------"
echo
echo

cd ..
R CMD INSTALL poisbinom

echo
echo
echo "'poisbinom' is installed."
echo "--------------------------"
echo
echo
echo "Run test code"
echo "--------------------------"
echo
echo

Rscript test/test_poisbinom.R --verbose


