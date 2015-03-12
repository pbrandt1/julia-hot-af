#!/bin/bash

START=`date +%s.%N`
./lib/julia-*/bin/julia parallel_test_serial.jl
END=`date +%s.%N`

echo "serial test took"
echo $END - $START | bc


START=`date +%s.%N`
./lib/julia-*/bin/julia parallel_test.jl
END=`date +%s.%N`

echo "parallel test took"
echo $END - $START | bc
