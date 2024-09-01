#!/bin/bash

# Copyright (C) 2023, Artificial Intelligence for the Sciences, Freie Universitaet Berlin, Germany
#
# This script is free software; you can redistribute it and/or modify
# it under the terms of the MIT license. See LICENSE for details.

LIG_DIR=$1
LIG_LIST=$2

while read line
do
        echo Preparing ligand $line
        prepare_ligand4.py \
                -l $LIG_DIR/$line.mol2 \
                -o $LIG_DIR/$line.pdbqt
        echo
done < $LIG_LIST
