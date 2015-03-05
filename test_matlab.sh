#!/bin/bash

START=`date +%s.%N`
./lib/octave*/run-octave *.m
END=`date +%s.%N`

echo $END - $START | bc
