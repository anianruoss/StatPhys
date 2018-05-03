#!/usr/bin/python
# convert coords.final file to xyz file


from sys import argv
import re

# import numpy as np
atom = 'Ar'

coords = re.compile(
        '^ +[0-9]+ +([0-9Ee\\+-.]+) +([0-9Ee\\+-.]+) +([0-9Ee\\+-.]+)')


def parser(line):
    c = coords.findall(line)[0]
    return ' '.join((atom, *c)) + '\n'


if len(argv) < 2:
    print('Usage: toxyz coords.final [coords.xyz]')
    exit(-1)

out = argv[1]+'.xyz' if len(argv) < 3 else argv[2]

# first swap first and second rows
with open(argv[1]) as f:
    lines = f.readlines()


out_lines = list(
        map(parser, lines[2:]))

out_lines.insert(0, lines[0])
out_lines.insert(0, str(len(out_lines)-1)+'\n')

with open(out, 'w') as of:
    of.writelines(out_lines)
