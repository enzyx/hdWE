#!/usr/bin/env python
from __future__ import absolute_import, division, print_function
import sys
from PIL import Image
import numpy
 
 
def main():
    image = Image.open(sys.argv[1])
#    if image.mode == '1':
    width, height = image.size
#    for i in xrange(height):
#        print(
#            ''.join(
#                '0' if image.getpixel((j, i)) else '1'
#                for j in xrange(width)
#            )
#        )

    matrix = numpy.zeros(shape=(height,width))
    colors = []
    # defined area codes, 1 = helix, 2 = beta sheet:
    # helix at phi,psi = -90,-45
    colors.append(image.getpixel((89,255)))
#    colors.append(image.getpixel((225,89)))
    # beta sheet phi,psi = -135, 135
    colors.append(image.getpixel((44,45)))
#    colors.append(image.getpixel((45,44)))
    fraction = 10
    for i in xrange(height):
        for j in xrange(width):
            if i%fraction==0 and j%fraction==0:
                pixel = image.getpixel((j,i))
                if not pixel in colors:
                    colors.append(pixel)
                matrix[i,j] = colors.index(pixel)
#                sys.stdout.write("{0:3d} ".format(image.getpixel((j,i))[0]))
#        if i%fraction == 0:
#            sys.stdout.write("\n")

    sys.stderr.write("found {} areas\n".format(len(colors)))
    
    sys.stdout.write("{} = [".format(sys.argv[2]))
    for i in range(0,height,fraction):
        sys.stdout.write("[")
        for j in range(0,width,fraction):
            sys.stdout.write(" {0:1d}".format(int(matrix[i,j])+1))
            if j != width-1:
                sys.stdout.write(",")
        sys.stdout.write("]")
        if i != height-1:
            sys.stdout.write(",")
            sys.stdout.write("\n")
    sys.stdout.write("]")

 
if __name__ == '__main__':
    main()
