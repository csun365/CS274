#!/bin/bash

for iteration_number in 100 500 1000; do
    output_file="${iteration_number}.log"
    for seed in $(seq 1 100); do
        python pvalue.py -r ${seed} -n ${iteration_number} data/drugs.csv data/targets.csv P54577 Q7RTX0 >> ${output_file}
    done
done