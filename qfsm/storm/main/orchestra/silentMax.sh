#!/bin/bash
bsub -q danuser_7d -R "span[hosts=1]" -n 1 "matlab -nodisplay -r 'runBatch;exit'"

