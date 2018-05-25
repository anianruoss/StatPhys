import configparser
from os import listdir, path

import matplotlib.pyplot as plt
import numpy as np
from natsort import natsorted

data_files = listdir('data')

temperature_files = [
    file for file in data_files if file.startswith('coords') and '_t' in file
]

c = configparser.ConfigParser()
c.read('mdpar.ini')
num_atom = int(c['md']['NumberAtoms'])

temperature_data = {'Base': {}, 'Harmonic': {}, 'Shake': {}}

for temperature_file in natsorted(temperature_files):
    simtype = temperature_file.split('coords_')[1].split('_t')[0]

    with open(path.join('data', temperature_file)) as file:
        temp = file.name.split('_t')[1].split('.')[0]
        temperature_data[simtype][temp] = []

        for line in file:
            if not line.startswith('MODEL'):
                break

            start_atom = file.readline()
            for i in range(num_atom - 2):
                file.readline()
            end_atom = file.readline()
            file.readline()

            start_atom = start_atom.split('UNK     1     ')[1]
            start_atom = start_atom.split('  1.00  0.00')[0]
            end_atom = end_atom.split('UNK     1     ')[1]
            end_atom = end_atom.split('  1.00  0.00')[0]

            start_coords = np.asarray(start_atom.split(' '), dtype=float)
            end_coords = np.asarray(end_atom.split(' '), dtype=float)

            temperature_data[simtype][temp].append(
                np.linalg.norm(start_coords - end_coords))

for simtype, data in temperature_data.items():
    fig, ax = plt.subplots()

    ax.set_title(simtype)
    ax.set_xlabel('NumMDSteps')
    ax.set_xlabel('NumMDSteps / TrajectoryOutputInterval')
    ax.set_ylabel('End-to-End Distance [$\AA$]')

    for temp, lengths in data.items():
        ax.plot(lengths, label=f'{temp}K')

    fig.savefig(f'{simtype}_temperature.png', bbox_inches='tight')
    plt.close(fig)
