#!/bin/bash
bsub -Is -q danuser_int_7d -R "span[hosts=1]" -n 8 "matlab -nodisplay -r 'runAlgorithm'"

