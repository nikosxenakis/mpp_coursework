#!/usr/bin/env python2
# test.py
import os
import filecmp

outputs = [
    ['./output/imagenew192x128.pgm', './test_output/imagenew192x128.pgm'],
    ['./output/imagenew256x192.pgm', './test_output/imagenew256x192.pgm'],
    ['./output/imagenew512x384.pgm', './test_output/imagenew512x384.pgm'],
    ['./output/imagenew768x768.pgm', './test_output/imagenew768x768.pgm'],
    ['./output/imagenew1600x1200.pgm', './test_output/imagenew1600x1200.pgm']
]

for output in outputs:
    if filecmp.cmp(output[0], output[1]) == False:
        print "ERROR: " + output[0]
    else:
        print "PASS: " + output[0]
