#!/bin/sh
g++ -O2 gauss.cpp -o gauss.exe -I ../../
./gauss.exe gauss.data
python3 read_headmap.py 4 5 gauss.data

