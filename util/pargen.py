#!/usr/bin/env python3
# generate hideous .inp file from
# readable configuration file

import configparser
from sys import argv

if len(argv) < 2:
    print('Usage: python pargen.py config_file {arg=val}')
    exit(-1)

c = configparser.ConfigParser()
c.read(argv[1])

for i in range(2, len(argv)):
    arg, val = argv[i].split('=')
    c['md'][arg.lower()] = val


config = [
    ('Title',),
    ('NumberAtoms', 'AtomicMass', 'MDType'),
    ('BoxSize(x)', 'BoxSize(y)', 'BoxSize(z)'),
    ('NumberMDSteps', 'InitialTime', 'TimeStep'),
    ('InitialTemperature', 'TargetTemperature', 'TemperatureCouplingTime'),
    ('RandomSeed',),
    ('XVInitialization', 'CoordInitializationSpread', 'FinalXVOutput'),
    ('NAtomsOnBoxEdge(x)', 'NAtomsOnBoxEdge(y)', 'NAtomsOnBoxEdge(z)'),
    ('EpsilonLJ', 'SigmaLJ', 'InteractionCutoffRadius'),
    ('PropertyPrintingInterval', 'NumberRadialDistrPoints',
        'RadialDistrCutoffRadius'),
    ('TrajectoryOutput', 'TrajectoryOutputFormat', 'TrajectoryOutputInterval')
]


def pad(string):
    return (27-len(string)) * " " + string


for line in config:
    text = '# '
    for option in line:
        text += pad(option)
    print(text)
    values = ' ' * 2
    for option in line:
        values += pad(c['md'].get(option))
    print(values)
