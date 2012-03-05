#!/bin/bash

for i in {1..25}
do
	./test.sh
done

gnuplot plot
