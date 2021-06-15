#!/bin/bash

libDir=lib
binDir=bin
srcDir=src
common=../library

totVaccine=$1
strategy=$2

name=V${totVaccine}S${strategy}

function debugBuild {
	g++ -std=c++17 -Wall -g \
        -I ${common} -I ${libDir}\
	    ${srcDir}/main-meanField.cpp\
        -o ${binDir}/${name}
}

function build {
	g++ -std=c++17 -O3 -flto -march=native\
        -I ${common} -I ${libDir}\
		${srcDir}/main-meanField.cpp\
        -o ${binDir}/${name}
}

#* Compile the source files
# debugBuild
build

#* Run
./${binDir}/${name} ${totVaccine} ${strategy}
rm ./${binDir}/${name}




