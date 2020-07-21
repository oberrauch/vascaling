#! /bin/bash

for file in $(find . -name "*.tex" ! -iname "thesis.tex"); do
	echo "Processing $file"
	sed '1d' $file > $file.mod
	mv $file.mod $file
	cat insert.txt $file > $file.mod
	mv $file.mod $file
done
