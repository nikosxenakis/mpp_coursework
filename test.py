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

if filecmp.cmp('./test_output/imagenew192x128.pgm', './output/imagenew192x128.pgm') == False:
    print "Error is not the same with previous runs"
else:
    print "It is the same with previous runs"
