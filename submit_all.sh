#!/bin/bash

mkdir "$TMPDIR"/FactorSim5/

cd "$HOME"/FactorSim5

sbatch -a 1-220 submit_jobs.sh


