#!/usr/bin/env python

import sys

import cairosvg

for argument in sys.argv[1:]:

    with open(argument) as f:
        svg = f.read()

    name = argument.split(".")[0]

    # save a PNG
    cairosvg.svg2png( bytestring=svg, scale=3.0, write_to=name+".png" )

