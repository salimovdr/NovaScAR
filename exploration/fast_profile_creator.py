import os
import pandas as pd
import numpy as np
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse
import json

# Путь к данным
path = os.path.expanduser('~/MyMACS2/data/ATAC/RG_aligns/data_GSE')

# Создание директорий
os.makedirs('profile', exist_ok=True)
os.makedirs('logs', exist_ok=True)

# Параметры командной строки
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

# Чтение данных
db_top100_cells = pd.read_csv(f'{path}/{name}_top_100_rg_coordinates_of_reads.csv')
db_top100_cells.columns = ['RN', 'chr', 'x', 'y', 'RG', 'strand1', 'strand2']
db_top100_cells = db_top100_cells[['RG', 'chr', 'x', 'y']]

# Ограничение по хромосомам
chrs = ['chr1', 'chr2']  # Например, только chr1
db_top100_cells = db_top100_cells[db_top100_cells['chr'].isin(chrs)]

# Группировка данных
RGs = list(db_top100_cells['RG'].unique())

print(f'Всего {len(RGs)} групп')

dict_top100_cells = {
    chr_value: {RG: db_rg.sort_values(by='y') for RG, db_rg in db_chr.groupby('RG')}
    for chr_value, db_chr in db_top100_cells.groupby('chr')
}

# Логирование
def write_log(log):
    print(log)
    with open(f'logs/log_profile_{name}.txt', "a") as file:
        file.write(log + '\n')

# Подсчёт интенсивности
def count_intensity(RG, chr):
    log = f'Начал считать интенсивность: RG={RG}, chr={chr}.'
    write_log(log)
    
    # Инициализация
    max_y = max(dict_top100_cells[chr][RG]['y']) + 1
    start_end_count_in_chr = np.zeros((max_y, 2), dtype=int)

    # Фильтрация и расчёт длины
    dict_top100_cells[chr][RG]['length'] = dict_top100_cells[chr][RG]['y'] - dict_top100_cells[chr][RG]['x']
    dict_top100_cells[chr][RG] = dict_top100_cells[chr][RG][dict_top100_cells[chr][RG]['length'] < 2000]

    # Заполнение интенсивности
    for _, row in dict_top100_cells[chr][RG][['x', 'y']].iterrows():
        x = row['x']
        y = row['y']
        start_end_count_in_chr[x, 0] += 1
        start_end_count_in_chr[y, 1] += 1

    # Накопительная сумма интенсивности
    intensity_nucleotides = np.cumsum(start_end_count_in_chr[:, 0] - np.pad(start_end_count_in_chr[:, 1][1:], (1, 0), 'constant'))

    log = f'Закончил считать интенсивность: RG={RG}, chr={chr}.'
    write_log(log)
    return intensity_nucleotides.tolist()  # Преобразуем в список для JSON

# Сохранение в JSON
def save_to_json(RG, chr, intensities):
    file_path = os.path.join('profile', f"{name}_{RG}_{chr}.json")
    data = {RG: {chr: {}}}
    data[RG][chr] = intensities
    with open(file_path, 'w') as f:
        json.dump(data, f, indent=4)
    write_log(f"Данные сохранены в {file_path} для RG={RG}, chr={chr}")

# Обработка RG и chr
def process_rg_chr(RG, chr):
    intensity_nucleotides = count_intensity(RG, chr)
    save_to_json(RG, chr, intensity_nucleotides)

# Параллельная обработка
def parallel_processing(chrs, RGs):
    tasks = []
    with ProcessPoolExecutor(max_workers=numcore) as executor:
        for RG in RGs[start:end + 1]:
            for chr in chrs:
                tasks.append(executor.submit(process_rg_chr, RG, chr))

        for _ in tqdm(as_completed(tasks), total=len(tasks), desc="Processing RGs and chromosomes"):
            pass

# Выполнение кода
parallel_processing(chrs, RGs)
write_log("Обработка завершена.")
