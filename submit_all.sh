#!/bin/bash

mkdir "$TMPDIR"/FactorSim/

cd "$HOME"/FactorSim

sbatch -a 1-240 submit_jobs.sh


