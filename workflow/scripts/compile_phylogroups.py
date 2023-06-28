#!/usr/bin/env python3

from logging.config import fileConfig
import pandas as pd

# define wildcard from input

# read file
input_list = list(set(snakemake.input))

for file in input_list:
    try:
        df = pd.read_table(file, sep='\t', header=None) 
        df.to_csv(snakemake.output[0], mode='a', header=False, index=False) # append to output csv
    except pd.errors.EmptyDataError:
        print("No data")