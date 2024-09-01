#!/bin/bash

# Copyright (C) 2023, Artificial Intelligence for the Sciences, Freie Universitaet Berlin, Germany
#
# This script is free software; you can redistribute it and/or modify
# it under the terms of the MIT license. See LICENSE for details.

REC_DIR=$1
REC_LIST=$2

while read line
do
        echo Preparing receptor $line
        prepare_receptor4.py \
                -r $REC_DIR/$line.pdb \
                -o $REC_DIR/$line.pdbqt
        echo
done < $REC_LIST
