#!/usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess

subprocess.call('g++ -o ash5tests-exe --std=c++11 ash5-tests.cpp ash5.cpp configuration.cpp', shell=True)
subprocess.call('./ash5tests-exe')
