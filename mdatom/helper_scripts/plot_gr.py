import sys

import matplotlib.pyplot as plt

with open(sys.argv[1], "r") as gr_file:
    plt.plotfile(gr_file, delimiter=' ', cols=(0, 1), skiprows=1, checkrows=0)

    plt.title('Radial Distribution Function')
    plt.xlabel('$r/\sigma$')
    plt.ylabel('$g(r)$', rotation=0, labelpad=20)

    plt.savefig(sys.argv[2])
