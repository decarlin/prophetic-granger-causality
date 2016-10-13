#!/usr/bin/env	python

import sys

for line in sys.stdin:
	parts = line.rstrip().split("\t")
	if len(parts) == 3:
		print "\t".join([parts[0], parts[2], parts[1]])
		

