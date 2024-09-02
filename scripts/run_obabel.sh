#!/bin/bash

# Copyright (C) 2023, Artificial Intelligence for the Sciences, Freie Universitaet Berlin, Germany
#
# This script is free software; you can redistribute it and/or modify
# it under the terms of the MIT license. See LICENSE for details.


DOCK_DIR=$1
ligand=$2

cd $DOCK_DIR/$ligand
for pdbqtfile in *.pdbqt
do
	if [ ! -f $pdbqtfile.pdb ]; then
		echo $pdbqtfile >> obabelout.log
		obabel -i pdbqt $pdbqtfile -o pdb -O tmp.pdb >> obabelout.log 2>&1
		grep -v "MODEL   " tmp.pdb > $pdbqtfile.pdb
		rm tmp.pdb
	fi
done
cd -
