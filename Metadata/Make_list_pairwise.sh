#!/bin/bash

pop_list="pop15inds.list"

# This script produces a list of all pairwise comparisons from a list
# The following command does so for population names (listed in the first column):
set -- $(cat ${pop_list} | cut -d$'\t' -f1)
for a; do
    shift
    for b; do
        printf "%s\t%s\n" "$a" "$b"
#        printf "%s\t%s\n" "$a","$b"
#               printf "$a - $b\n"
    done
done > pop_name_pairs
