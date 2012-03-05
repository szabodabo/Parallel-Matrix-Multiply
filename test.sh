#!/bin/bash

for NODES in 2 3 4 6 7 8 12 14 16
do
	echo "Running main with $NODES cores..."
	mpirun -np $NODES ./main
done
