#!/bin/bash

Rscript -e "remove.packages('pbquick')"

echo "--------------------------"
echo "Script to compile 'pbquick' and run test code"
echo

if [ -e "pbquick" ]&& [ -d "pbquick" ]; then
    rm -rf pbquick
    echo
    echo "Old 'pbquick' folder is removed."
    echo "--------------------------"
    echo
fi


Rscript -e "library(Rcpp); Rcpp.package.skeleton('pbquick')"

echo
echo
echo "Skeleton package is created."
echo "--------------------------"
echo
echo

cp src/dpbquick.cpp pbquick/src/
rm pbquick/src/rcpp_hello_world.cpp

echo
echo
echo "C++ codes are copied."
echo "--------------------------"
echo
echo

cd pbquick
Rscript -e "library(Rcpp); compileAttributes(verbose=TRUE)"
rm man/pbquick-package.Rd

echo
echo
echo "'pbquick' is compiled."
echo "--------------------------"
echo
echo

cd ..
R CMD INSTALL pbquick

echo
echo
echo "'pbquick' is installed."
echo "--------------------------"
echo
echo
echo "Run test code"
echo "--------------------------"
echo
echo

Rscript test/test_pbquick.R --verbose


