#!/bin/bash

START=`date +%s.%N`
./lib/julia-*/bin/julia *.jl
END=`date +%s.%N`

echo $END - $START | bc
