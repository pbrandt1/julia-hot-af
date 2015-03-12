#!/bin/bash

START=`date +%s.%N`
./lib/julia-*/bin/julia julia_gravity_test.jl
END=`date +%s.%N`

echo $END - $START | bc
