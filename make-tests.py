#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import subprocess

# "C:\MinGW\MinGW\bin\g++.exe" -o ash5tests-exe --std=c++11 ash5tests.cpp ash5.cpp -I C:\cpp-libs\boost_1_65_0\ -I C:\cpp-libs\Eigen3
if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Usage: {} pathToG++ pathToEigen pathToBoost".format(sys.argv[0]))
	else:
		gccPath = sys.argv[1]
		eigenPath = sys.argv[2]
		boostPath = sys.argv[3]
		subprocess.call('{} -o ash5tests-exe --std=c++11 ash5tests.cpp ash5.cpp -I {}/ - I {}'.format(gccPath, eigenPath), shell=True)
		subprocess.call('ash5tests-exe')
