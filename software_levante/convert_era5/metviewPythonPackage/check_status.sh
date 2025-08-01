#!/bin/bash

if [ -z $(squeue -u b382762 -o %.22j | grep $1) ]; then echo 'Finished'; else echo 'running'; fi

