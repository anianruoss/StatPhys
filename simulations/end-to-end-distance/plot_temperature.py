from os import listdir, path

import matplotlib.pyplot as plt
import numpy as np

data_files = listdir('data')

temperature_files = [
    file for file in data_files if file.startswith('coords') and '_t' in file
]

temperature_data = {'Base': {}, 'Harmonic': {}, 'Shake': {}}

for temperature_file in sorted(temperature_files):
    simtype = temperature_file.split('coords_')[1].split('_t')[0]

    with open(path.join('data', temperature_file)) as file:
        temp = file.name.split('_t')[1].split('.')[0]
        temperature_data[simtype][temp] = []

        for line in file:
            if not line.startswith('MODEL'):
                break

            start_atom = file.readline()
            for i in range(6):
                file.readline()
            end_atom = file.readline()
            file.readline()

            start_atom = start_atom.split('1      ')[1].split('  1.00')[0]
            end_atom = end_atom.split('1      ')[1].split('  1.00')[0]

            start_atom = start_atom.split(' ')
            end_atom = end_atom.split(' ')

            start_coords = np.asarray(
                [start_atom[0], start_atom[2], start_atom[4]],
                dtype=float
            )
            end_coords = np.asarray(
                [end_atom[0], end_atom[2], end_atom[4]],
                dtype=float
            )

            temperature_data[simtype][temp].append(
                np.linalg.norm(start_coords - end_coords))

for simtype, data in temperature_data.items():
    fig, ax = plt.subplots()

    ax.set_title(simtype)
    ax.set_xlabel('NumMDSteps')
    ax.set_ylabel('End-to-End Distance')

    for temp, lengths in data.items():
        ax.plot(lengths, label=f'{temp}K')

    ax.legend()

    fig.show()
    fig.savefig(f'{simtype}_temperature.pdf', bbox_inches='tight')
