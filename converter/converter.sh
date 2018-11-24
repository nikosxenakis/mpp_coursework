#!/bin/bash
#
#this script converts the given image.jpg to a new edge image for reconstruction
convert archer.jpg -colorspace gray -compress none -depth 8 image400x266.pgm
./pgm2edge image400x266.pgm edgenew400x266.pgm
