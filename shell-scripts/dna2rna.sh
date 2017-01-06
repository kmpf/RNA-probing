#!/bin/bash

for file in $(ls *.fa)
do
    grep '^>' $file > $(basename $file .fa).rna.fa
    grep '^[^>]' $file | tr Tt Uu >> $(basename $file .fa).rna.fa
done
