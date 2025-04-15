import gsm_automated_plots
from hail.utils import hadoop_open
import json
# import hail as hl
import pandas as pd
import numpy as np
import argparse

REQUIRED_COMMAND_LIST = [
            'plot_title', 
            'eval_name', 
            'results_path', 
            'score_data_path', 
            'score_is_gsm', 
            'score_name_dict', 
            'score_join_on', 
            'score_delim', 
            'score_is_ht', 
            'eval_data_path', 
            'eval_join_on', 
            'eval_delim', 
            'eval_is_ht'
            ]
OPTIONAL_COMMAND_LIST = ['results_path', 'obs/exp', 'cutoffs']

def fix_dict(d): 
    if d['score_delim'] == None: 
            d['score_delim'] = '\t'
    if d['eval_delim'] == None: 
        d['eval_delim'] = '\t'
    for c in OPTIONAL_COMMAND_LIST:
        if c not in list(d.keys()): 
            print('getting deleted 1')
            d[c] = None
    if d['results_path'] == None: 
         d['results_path'] = './results/'
    return d

def find_key(gsm_labels, gle_labels): 
    key_priority = ['enst', 'ensg', 'gene_id', 'uniprot_id']
    keys_in_common = set(gsm_labels.keys()) & set(gle_labels.keys())
    print(keys_in_common)
    for k in key_priority: 
         if k in keys_in_common: 
              return k
    return None

def add_filters(all_d, names, type):
    names_copy = names.copy()
    for name in names_copy:
        info = all_d[name]
        print(info)
        if "filter" in info and info["filter"]:
            for  k, v in info["filter"].items():
                new_info = info.copy()
                new_info[f'{type}_filter'] = v
                all_d[f'{name}_{k}'] = new_info
                names.append(f'{name}_{k}')
    return all_d, names
                     
        

def main(scores_d, evals_d, GSM_list, GLE_list):
    all_results = pd.DataFrame()
    
    scores_d, GSM_list = add_filters(scores_d, GSM_list, 'score')
    evals_d, GLE_list = add_filters(evals_d, GLE_list, 'eval')

    print('looping')
    for gle in GLE_list:  
        for gsm in GSM_list: 
            d = {}
            gsm_d, gle_d = scores_d[gsm], evals_d[gle]
            
            print(d)
            d = gsm_d | gle_d
            join_on = find_key(gsm_d['labels'],gle_d['labels'])
            if not join_on: 
                raise Exception(f'No matching keys between {gsm} and {gle}')
            d['join_on'], d['score_join_on'], d['eval_join_on'] = join_on, gsm_d['labels'][join_on], gle_d['labels'][join_on]
            print(d)
            d = fix_dict(d)
            d['plot_title'] = gsm
            d['eval_name'] = gle
            d['cutoffs'] = [0.9, 0.95]

            results, stats = gsm_automated_plots.main(d)
            rename_dict = {s: f'{gle}_{s}' for s in stats}
            results = results.rename(columns=rename_dict)
            all_results = pd.concat([all_results, results], axis=1)

    path = 'results/.'
    all_results.to_csv(f'{path}multi_gene_eval.tsv', sep='\t', index=True)
    print(all_results)

    print('done')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--scores', type=str, required=True)
    parser.add_argument('--evals', type=str, required=True)
    parser.add_argument('--score-subset', type=str, help='To only do a subset of scores in the scores json file, add comma delimited names of scores')
    parser.add_argument('--eval-subset', type=str, help='To only do a subset of evaluation lists in the evals json file, add comma delimited names of evaluations')
    
    print(parser)
    args = parser.parse_args()


    scores_d = json.load(hadoop_open(args.scores))
    evals_d = json.load(hadoop_open(args.evals))

    if 'GSM' not in scores_d: 
        raise Exception('No gene-level scores provided')
    if 'GLE' not in evals_d: 
        raise Exception('No gene-level evaluation data provided')
    if args.score_subset: 
         GSM_list = args.score_subset.split(',')
         for s in GSM_list: 
              if s not in scores_d['GSM']: 
                   raise Exception('Score name not in file')  
    else: 
        GSM_list =  scores_d['GSM']
    if args.eval_subset: 
         GLE_list = args.eval_subset.split(',')
         for e in GLE_list: 
              if e not in evals_d['GLE']: 
                   raise Exception('Eval name not in file')
    else: 
         GLE_list = evals_d['GLE']
    print('ready')

    main(scores_d['GSM'], evals_d['GLE'], GSM_list, GLE_list)