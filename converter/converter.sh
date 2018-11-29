#!/bin/bash
#
#this script converts the given image.jpg 400 x 266 to a new edge image for reconstruction
convert image.jpg -colorspace gray -compress none -depth 8 image400x266.pgm
./pgm2edge image400x266.pgm edgenew400x266.pgm
