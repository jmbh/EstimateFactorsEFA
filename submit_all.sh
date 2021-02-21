#!/bin/bash

mkdir "$TMPDIR"/FactorSi4/

cd "$HOME"/FactorSim4

sbatch -a 101-200 submit_jobs.sh
