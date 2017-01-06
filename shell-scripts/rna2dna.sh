#!/bin/bash

for file in $(ls *.fa); do
    filename=$(basename $file .fa).dna.fa
    grep '^>' $file > $filename
    grep '^[^>]' $file | tr Uu Tt >> $filename
done
