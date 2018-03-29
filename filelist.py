#!/usr/bin/python2

import os

outfile = open("filelists/defaultlist.txt", 'w')

for filename in os.listdir("data"):
    outfile.write("{}\n".format(filename))

outfile.close()
