import os
import sys
import gcsfs
fs = gcsfs.GCSFileSystem()

current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_dir, '..'))

import pandas as pd
from evaluation.compute_metrics import (
    enrichments_from_grouped_contingency_table,
    results_from_full_table_many_cutoffs
)
from analyses_for_manuscript.scores_to_use import PERCENTILE_SCORES, RAW_SCORES, RAW_SCORES_NO_PAI3D, PERCENTILE_SCORES_GENE_MEDIAN, SCORE_LIST_MAP, SUBSET_INFO
import argparse
import json
print('done imports')

def calc_write_grouped_enrichments(df, pval_cutoff, analysis_path, subset_df_path = None, subset_df_colname = None): 
    df_results_all_cutoffs, df_w_std_error, df_p = enrichments_from_grouped_contingency_table(df, pval_cutoff)
    df_p.to_csv(analysis_path + '_pvals.tsv', sep='\t')
    df_results_all_cutoffs.to_csv(analysis_path + '_results_all_cutoffs.tsv', sep='\t')
    df_w_std_error.to_csv(analysis_path + '_results_with_std_errors.tsv', sep='\t', index=False)


def main(
    analysis_name,
    is_pos_name,
    score_column_list,
    cutoff_list,
    pval_cutoff,
    hail_contingency,
    stat,
    ncase,
    ncontrol,
    do_subsets, #yes, no, or only
    gene_median,
    file_path
):
    analysis_path = file_path + analysis_name
    if hail_contingency:
        df_path = f'gs://genetics-gym/hail_evaluation_tables/{analysis_name}_contingency_tables.tsv'
        df = pd.read_csv(df_path, sep='\t')

        if do_subsets == 'yes' or do_subsets == 'no':  
            calc_write_grouped_enrichments(df, pval_cutoff, analysis_path)
        if do_subsets == 'yes' or do_subsets == 'only': 
            for subset_name, (subset_df_path, subset_df_colname) in SUBSET_INFO.items():
                print(subset_df_path, subset_df_colname)
                calc_write_grouped_enrichments(df, pval_cutoff, analysis_path + f"_{subset_name}", subset_df_path, subset_df_colname)
        
    else:
        df_path = f'gs://genetics-gym/hail_evaluation_tables/{analysis_name}.tsv'
        df = pd.read_csv(df_path, sep='\t')
        df_w_std_error, df_p = results_from_full_table_many_cutoffs(df, is_pos_name, score_column_list, cutoff_list, pval_cutoff, stat, ncase, ncontrol)
        df_p.to_csv(analysis_path + '_pvals.tsv', sep='\t')
        df_w_std_error.to_csv(analysis_path + '_results_with_std_errors.tsv', sep='\t', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--json', type=str)
    parser.add_argument('--do-subsets', type=str, default = 'no')
    parser.add_argument('--gene-median', action='store_true', default=False)
    args = parser.parse_args()

    for json_file in args.json.split(','):

        print('reading', json_file)
        d = json.load(fs.open(json_file))
        print('read')

        d['analysis_name'] = json_file.split('/')[-1][0:-5]
        if d['score_column_list'] in SCORE_LIST_MAP.keys(): 
            d['score_column_list'] = SCORE_LIST_MAP[d['score_column_list']]
        else:
            d['score_column_list'] = d['score_column_list'].split(',')

        for x in ['subset_df_path', 'subset_df_colname', 'stat', 'ncase', 'ncontrol']:
            if not x in d:
                d[x] = None
        
        d['file_path'] = "gs://trisha-tmp/recreating_workflow/"

        print('arguments:')
        print(json.dumps(d, indent=4))
        main(
            d['analysis_name'],
            d['is_pos_name'],
            d['score_column_list'],
            d['cutoff_list'],
            d['pval_cutoff'],
            d['hail_contingency'],
            d['stat'],
            d['ncase'],
            d['ncontrol'],
            args.do_subsets,
            args.gene_median,
        )

