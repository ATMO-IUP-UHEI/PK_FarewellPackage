#!/bin/bash

if [ -z $(squeue -u $USER -o %.22j | grep $1) ]; then echo 'Finished'; else echo 'running'; fi

