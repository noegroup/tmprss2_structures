#!/bin/bash

# Copyright (C) 2023, Artificial Intelligence for the Sciences, Freie Universitaet Berlin, Germany
#
# This script is free software; you can redistribute it and/or modify
# it under the terms of the MIT license. See LICENSE for details.


LIG_DIR=$1
REC_DIR=$2
DOCK_DIR=$3

ligand=$4
receptor=$5

if [ ! -f $DOCK_DIR/$ligand/$receptor.pdbqt ]; then

mkdir -p $DOCK_DIR/$ligand
chmod -R g+rw $DOCK_DIR/$ligand

x=$(grep 'CA  SER A 186' $REC_DIR/$receptor.pdbqt | awk '{print $7}')
y=$(grep 'CA  SER A 186' $REC_DIR/$receptor.pdbqt | awk '{print $8}')
z=$(grep 'CA  SER A 186' $REC_DIR/$receptor.pdbqt | awk '{print $9}')
smina \
	-r $REC_DIR/$receptor.pdbqt \
	-l $LIG_DIR/$ligand.pdbqt \
	--center_x $x \
	--center_y $y \
	--center_z $z \
	--size_x 30 \
	--size_y 30 \
	--size_z 30 \
	--scoring 'vinardo' \
	-o $DOCK_DIR/$ligand/$receptor.pdbqt \
	--log $DOCK_DIR/$ligand/$receptor.log \
	--seed 42 \
	--exhaustiveness 10 \
	--cpu 1 \
	--flexres A:44,A:45,A:180,A:183,A:206 \
	--out_flex $DOCK_DIR/$ligand/$receptor\_flex.pdbqt
fi
