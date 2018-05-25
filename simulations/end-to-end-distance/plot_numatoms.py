from os import listdir, path

import matplotlib.pyplot as plt
import numpy as np
from natsort import natsorted

data_files = listdir('data')

numatom_files = [
    file for file in data_files if file.startswith('coords') and '_n' in file
]

numatom_data = {'Base': {}, 'Harmonic': {}, 'Shake': {}}

for numatom_file in natsorted(numatom_files):
    simtype = numatom_file.split('coords_')[1].split('_n')[0]

    with open(path.join('data', numatom_file)) as file:
        num_atoms = int(file.name.split('_n')[1].split('.')[0])
        numatom_data[simtype][num_atoms] = []

        for line in file:
            if not line.startswith('MODEL'):
                break

            start_atom = file.readline()
            for i in range(num_atoms - 2):
                file.readline()
            end_atom = file.readline()
            file.readline()

            start_atom = start_atom.split('UNK     1     ')[1]
            start_atom = start_atom.split('  1.00  0.00')[0]
            end_atom = end_atom.split('UNK     1     ')[1]
            end_atom = end_atom.split('  1.00  0.00')[0]

            start_coords = np.asarray(start_atom.split(' '), dtype=float)
            end_coords = np.asarray(end_atom.split(' '), dtype=float)

            numatom_data[simtype][num_atoms].append(
                np.linalg.norm(start_coords - end_coords))

for simtype, data in numatom_data.items():
    fig, ax = plt.subplots()

    ax.set_title(simtype)
    ax.set_xlabel('NumMDSteps / TrajectoryOutputInterval')
    ax.set_ylabel('End-to-End Distance [$\AA$]')

    for temp, lengths in data.items():
        ax.plot(lengths, label=f'{temp} Atoms')

    fig.savefig(f'{simtype}_numatom.png', bbox_inches='tight')
    plt.close(fig)
