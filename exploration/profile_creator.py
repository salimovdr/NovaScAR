import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import poisson
from statsmodels.stats.multitest import multipletests
import numpy as np
import json
from tqdm import tqdm
import time
from datetime import datetime
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

path = '~/MyMACS2/data/ATAC/RG_aligns/data_GSE' 
name = 'GSM7790861'

db_top100_cells = pd.read_csv(f'{path}/{name}_top_100_rg_coordinates_of_reads.csv')
db_top100_cells.columns = ['RN', 'chr', 'x', 'y', 'RG', 'strand1', 'strand2']
db_top100_cells = db_top100_cells[['RG', 'chr', 'x', 'y']]

max_x = max(db_top100_cells['y'])
chrs = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
db_top100_cells = db_top100_cells[db_top100_cells['chr'].isin(chrs)]

db_peaks = pd.read_csv(f'{path}/{name}_peaks.xls', sep='\t', comment='#')
db_peaks = db_peaks.rename(columns={'start': 'x', 'end': 'y', 'abs_summit': 'x_summit', 'name': 'id'})
db_peaks = db_peaks[db_peaks['chr'].isin(chrs)]

max_x = max(max_x, max(db_peaks['y']))

RGs = list(db_top100_cells['RG'].unique())

dict_peaks = {}
for chr_value, db_chr in db_peaks.groupby('chr'):
    dict_peaks[chr_value] = {RG: db_chr.sort_values(by='y') for RG in RGs}

dict_top100_cells = {}
for chr_value, db_chr in db_top100_cells.groupby('chr'):
    dict_top100_cells[chr_value] = {RG: db_rg.sort_values(by='y') for RG, db_rg in db_chr.groupby('RG')}

def count_intensity(RG, chr):
    start_end_count_in_chr = [[0, 0] for _ in range(max(max(dict_peaks[chr][RG]['y']), max(dict_top100_cells[chr][RG]['y'])) + 1)]
    
    dict_top100_cells[chr][RG]['length'] = dict_top100_cells[chr][RG]['y'] - dict_top100_cells[chr][RG]['x']
    dict_top100_cells[chr][RG] = dict_top100_cells[chr][RG][dict_top100_cells[chr][RG]['length'] < 5000]
    
    for index, row in dict_top100_cells[chr][RG][['x', 'y']].iterrows():
        x = row['x']
        y = row['y']
        start_end_count_in_chr[x][0] += 1
        start_end_count_in_chr[y][1] += 1

    intensity_nucleotides = [0] * len(start_end_count_in_chr)

    current_intensity = 0
    for i in range(len(start_end_count_in_chr)):
        if i != 0:
            current_intensity += (start_end_count_in_chr[i][0] - start_end_count_in_chr[i - 1][1])
        else:
            current_intensity += start_end_count_in_chr[i][0]
        intensity_nucleotides[i] = current_intensity

    intensity_nucleotides.append(0)
    return (RG, chr, intensity_nucleotides)

def process_RG(RG):
    result = {}
    for chr in chrs:
        _, _, intensity = count_intensity(RG, chr)
        result[chr] = intensity
    return RG, result

intensity_nucleotides = {}
with ProcessPoolExecutor() as executor:
    futures = {executor.submit(process_RG, RG): RG for RG in RGs}
    for future in tqdm(futures):
        RG, result = future.result()
        intensity_nucleotides[RG] = result

with open(f'profile/{name}_intensity_nucleotides.json', 'w') as file:
    json.dump(intensity_nucleotides, file)