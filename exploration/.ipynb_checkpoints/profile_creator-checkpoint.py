import os
import json
import pandas as pd
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse

path = '~/MyMACS2/data/ATAC/RG_aligns/data_GSE' 

os.makedirs('profile', exist_ok=True)
if not os.path.exists(f'logs/'):
    os.makedirs(f'logs/')

parser = argparse.ArgumentParser()
parser.add_argument("--name", type=str, required=True)
parser.add_argument("--start", type=int, required=True)
parser.add_argument("--end", type=int, required=True)
parser.add_argument("--numcore", type=int, required=True)

args = parser.parse_args()
name = args.name
start = args.start
end = args.end
numcore = args.numcore


db_top100_cells = pd.read_csv(f'{path}/{name}_top_100_rg_coordinates_of_reads.csv')
db_top100_cells.columns = ['RN', 'chr', 'x', 'y', 'RG', 'strand1', 'strand2']
db_top100_cells = db_top100_cells[['RG', 'chr', 'x', 'y']]

#chrs = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
chrs = ['chr1']
db_top100_cells = db_top100_cells[db_top100_cells['chr'].isin(chrs)]

RGs = list(db_top100_cells['RG'].unique())

dict_top100_cells = {}
for chr_value, db_chr in db_top100_cells.groupby('chr'):
    dict_top100_cells[chr_value] = {RG: db_rg.sort_values(by='y') for RG, db_rg in db_chr.groupby('RG')}

def write_log(log):
    print(log)
    with open(f'logs/log_profile_{name}.txt', "a") as file:
        file.write(log + '\n')

def count_intensity(RG, chr):
    log = f'Начал считать интенсивность: RG={RG}, chr={chr}.'
    write_log(log)
    start_end_count_in_chr = [[0, 0] for _ in range(max(dict_top100_cells[chr][RG]['y']) + 1)]
    
    dict_top100_cells[chr][RG]['length'] = dict_top100_cells[chr][RG]['y'] - dict_top100_cells[chr][RG]['x']
    dict_top100_cells[chr][RG] = dict_top100_cells[chr][RG][dict_top100_cells[chr][RG]['length'] < 2000]
    
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
    log = f'Закончил считать интенсивность: RG={RG}, chr={chr}.'
    write_log(log)
    return intensity_nucleotides

def save_to_json(RG, chr, intensities):
    file_path = os.path.join('profile', f"{name}_{RG}_{chr}.json")
    data = {RG: {chr: intensities}}
    with open(file_path, 'w') as f:
        json.dump(data, f)

def merge_json_files():
    all_data = {}

    for file_name in os.listdir('profile'):
        if file_name.startswith(name):
            file_path = os.path.join('profile', file_name)
            
            with open(file_path, 'r') as f:
                data = json.load(f)
                all_data.update(data)

    with open(f'{name}_profile.json', 'w') as f:
        json.dump(all_data, f, indent=4)

def process_rg_chr(RG, chr):
    intensity_nucleotides = count_intensity(RG, chr)
    save_to_json(RG, chr, intensity_nucleotides)
    

def parallel_processing(chrs, RGs):
    tasks = []
    with ProcessPoolExecutor(max_workers=numcore) as executor:
        for RG in RGs[start:end + 1]:
            for chr in chrs:
                tasks.append(executor.submit(process_rg_chr, RG, chr))

        for _ in tqdm(as_completed(tasks), total=len(tasks), desc="Processing RGs and chromosomes"):
            pass

parallel_processing(chrs, RGs)

merge_json_files()
