#!/bin/bash
# Run driver_v3.py on runs/lhs_hf/config_1.yaml ... config_30.yaml in turn.
# Run from the folder holding driver_v3.py:  bash run_batch.sh

python3.13 make_yaml_file.py 

for i in $(seq 1 100); do
    python3.13 driver.py runs/lhs_hf/config_${i}.yaml
done
