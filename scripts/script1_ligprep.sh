#!/bin/bash

# Copyright (C) 2023, Artificial Intelligence for the Sciences, Freie Universitaet Berlin, Germany
#
# This script is free software; you can redistribute it and/or modify
# it under the terms of the MIT license. See LICENSE for details.


SCHRODINGER=/import/ag_cmb/software/schrodinger2020-2
LIG_PREP=$SCHRODINGER/ligprep
MOL2_CONV=$SCHRODINGER/utilities/mol2convert

OUT_DIR=$1
SMI_FILE=$2

while read -r name smiles
do
	echo $smiles > $OUT_DIR/$name.smi
	$LIG_PREP -WAIT -ismi $OUT_DIR/$name.smi -omae $OUT_DIR/$name.mae -epik -We,-ph,7.05,-ms,1 -s 1 -nt -g
	rm $OUT_DIR/$name.smi
	$MOL2_CONV -imae $OUT_DIR/$name.mae -omol2 $OUT_DIR/$name.mol2
done < $SMI_FILE
