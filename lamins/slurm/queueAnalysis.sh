#!/bin/bash
mv job_* old
for l in `cat lists.txt`;
    do sbatch ~/slurm/analyzeLamins_20150610.sh $l;
done;
