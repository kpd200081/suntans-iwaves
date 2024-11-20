#!/bin/bash

cd rundata
octave ../../suntans/mfiles/onedgrid.m
cd ../
make clobber
