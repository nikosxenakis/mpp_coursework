#!/usr/bin/env python2
# test.py
import os
import filecmp

# for test_output in os.listdir('./test_output'):
#     if test_output.endswith('.pgm'):
#         for output in os.listdir('./output'):
#             if output.endswith('.pgm'):
#                 if test_output == output:
#                     if filecmp.cmp('./test_output/'+test_output, './output/'+output) == False:
#                         print "Error in: " + output + "\n"

# if filecmp.cmp('./test_output/image192x128.pgm', './output/imagenew192x128.pgm') == False:
#     print "Error is not the same with the case study"
# else:
#     print "It is the same with the case study"

outputs = [
    ['./output/imagenew192x128.pgm', './test_output/image192x128.pgm'],
    ['./output/imagenew256x192.pgm', './test_output/image256x192.pgm'],
    ['./output/imagenew512x384.pgm', './test_output/image512x384.pgm']
    # ['./output/imagenew768x768.pgm', './test_output/image768x768.pgm']
]

for output in outputs:
    if filecmp.cmp(output[0], output[1]) == False:
        print "ERROR: " + output[0]
    else:
        print "PASS: " + output[0]
