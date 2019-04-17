#!/bin/bash
trash callgrind*
make
valgrind --tool=callgrind ./leath 10000
kcachegrind callgrind*
