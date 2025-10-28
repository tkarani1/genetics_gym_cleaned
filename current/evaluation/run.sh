#/bin/bash

json=gs://genetics-gym/ddd.json,gs://genetics-gym/asd.json

hailctl dataproc start hkf --network=broad-allow --max-idle 1h
hailctl dataproc submit hkf evaluation/1_evaluate.py --json $json --pyfiles evaluation,utils,analyses_for_manuscript
python -m evaluation.2_evaluate --json $json


# hailctl dataproc submit hkf evaluation/make_hail_tables.py --json gs://genetics-gym/json/gnomad_ungrouped.json --use-checkpoint gnomad_grouped --pyfiles evaluation,utils,analyses_for_manuscript
# python -m evaluation.process_hail_tables --json gs://genetics-gym/json/asd.json,gs://genetics-gym/json/chd.json,gs://genetics-gym/json/clinvar_grouped.json,gs://genetics-gym/json/ddd.json,gs://genetics-gym/json/fine_mapped.json,gs://genetics-gym/json/genebass.json,gs://genetics-gym/json/gnomad_grouped.json
# python -m evaluation.process_hail_tables --json gs://genetics-gym/json/gnomad_grouped.json,gs://genetics-gym/json/clinvar_grouped.json --do-subsets only

## If running on hail batch
# make sure environment is properly configured to have hail batch backend
# python evaluation/1_evaluate.py --json $json --pyfiles evaluation,utils,analyses_for_manuscript
# python evaluation/2_evaluate.py --json $json

## If running on dataproc. To limit dataproc use, only run evaluate_1.py on the cluster. Then run evaluate_2.py locally
# start dataproc cluster 
# hailctl dataproc start hkf --network=broad-allow --max-idle 1h
