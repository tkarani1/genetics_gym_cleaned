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
from analyses_for_manuscript.scores_to_use import PERCENTILE_SCORES, RAW_SCORES, RAW_SCORES_NO_PAI3D, PERCENTILE_SCORES_GENE_MEDIAN
import argparse
import json
print('done imports')


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
):
    if hail_contingency:
        if not do_subsets == 'only':
            df_path = f'gs://genetics-gym/hail_evaluation_tables/{analysis_name}_contingency_tables.tsv'
            df = pd.read_csv(df_path, sep='\t')
            df_results_all_cutoffs, df_w_std_error, df_p = enrichments_from_grouped_contingency_table(df, pval_cutoff)
            df_p.to_csv(f'gs://genetics-gym/results/{analysis_name}_pvals.tsv', sep='\t')
            df_results_all_cutoffs.to_csv(f'gs://genetics-gym/results/{analysis_name}_results_all_cutoffs.tsv', sep='\t')
            df_w_std_error.to_csv(f'gs://genetics-gym/results/{analysis_name}_results_with_std_errors.tsv', sep='\t', index=False)
        if not do_subsets == 'no':
            subset_info = (
                ('surfy', 'gs://missense-scoring/surfaceome_ids.txt', 'uniprot_name'),
                ('gpcrs', 'gs://missense-scoring/macarthur_gene_lists/gpcr_union.tsv', 'gene_symbol'),
                ('dominant', 'gs://missense-scoring/macarthur_gene_lists/all_ad.tsv',  'gene_symbol'),
                ('recessive', 'gs://missense-scoring/macarthur_gene_lists/all_ar.tsv', 'gene_symbol'),
                ('kinases', 'gs://missense-scoring/macarthur_gene_lists/kinases.tsv', 'gene_symbol'),
                ('brain', 'gs://genetics-gym/linkers-and-annotations/brain_list.txt', 'uniprot_id'),
                ('low_loeuf', 'gs://genetics-gym/linkers-and-annotations/low_loeuf_list.txt', 'gene_symbol'),
                ('high_loeuf', 'gs://genetics-gym/linkers-and-annotations/high_loeuf_list.txt', 'gene_symbol'),   
            )
            for (subset_name, subset_df_path, subset_df_colname) in subset_info:
                print(subset_df_path, subset_df_colname)
                df_path = f'gs://genetics-gym/hail_evaluation_tables/{analysis_name}_contingency_tables.tsv'
                df = pd.read_csv(df_path, sep='\t')
                df_results_all_cutoffs, df_w_std_error, df_p = enrichments_from_grouped_contingency_table(df, pval_cutoff, subset_df_path, subset_df_colname)
                df_p.to_csv(f'gs://genetics-gym/results/{analysis_name}_{subset_name}_pvals.tsv', sep='\t')
                df_results_all_cutoffs.to_csv(f'gs://genetics-gym/results/{analysis_name}_{subset_name}_results_all_cutoffs.tsv', sep='\t')
                df_w_std_error.to_csv(f'gs://genetics-gym/results/{analysis_name}_{subset_name}_results_with_std_errors.tsv', sep='\t', index=False)
        
    else:
        df_path = f'gs://genetics-gym/hail_evaluation_tables/{analysis_name}.tsv'
        df = pd.read_csv(df_path, sep='\t')
        df_w_std_error, df_p = results_from_full_table_many_cutoffs(df, is_pos_name, score_column_list, cutoff_list, pval_cutoff, stat, ncase, ncontrol)
        df_p.to_csv(f'gs://genetics-gym/results/{analysis_name}_pvals.tsv', sep='\t')
        df_w_std_error.to_csv(f'gs://genetics-gym/results/{analysis_name}_results_with_std_errors.tsv', sep='\t', index=False)

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
        if d['score_column_list'] == "PERCENTILE_SCORES":
            d['score_column_list'] = PERCENTILE_SCORES
        if d['score_column_list'] == "PERCENTILE_SCORES_GENE_MEDIAN":
            d['score_column_list'] = PERCENTILE_SCORES_GENE_MEDIAN
        if d['score_column_list'] == "RAW_SCORES_NO_PAI3D":
            d['score_column_list'] = RAW_SCORES_NO_PAI3D
        if d['score_column_list'] == "RAW_SCORES":
            d['score_column_list'] = RAW_SCORES

        for x in ['subset_df_path', 'subset_df_colname', 'stat', 'ncase', 'ncontrol']:
            if not x in d:
                d[x] = None

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

                    # elif d['evaluation'] == 'auc':
                    #         main_compute_auc(
                    #                 d['analysis_name'],
                    #                 d['ht_evaluation_path'],
                    #                 d['ht_score_path'],
                    #                 d['is_pos_name'],
                    #                 d['filter_table_list'],
                    #                 d['filter_column_list'],
                    #                 d['score_column_list'],
                    #                 d['group_by'],                                
                    #         )
                    # elif d['evaluation'] == 'abs_spearman':
                    #         main_abs_spearman_grouped(
                    #                 d['analysis_name'],
                    #                 d['ht_evaluation_path'],
                    #                 d['ht_score_path'],
                    #                 d['outcome_column'],
                    #                 d['filter_table_list'],
                    #                 d['filter_column_list'],
                    #                 d['score_column_list'],
                    #                 d['group_by'],                                       
                    #         )


    # def main_abs_spearman_grouped(
    #         analysis_name,
    #         ht_evaluation_path,
    #         ht_score_path,
    #         outcome_column,
    #         filter_table_list,
    #         filter_column_list,
    #         score_column_list,
    #         group_by,
    # ):
    #         ht = load_join_filter(
    #                 ht_evaluation_path,
    #                 ht_score_path,
    #                 outcome_column,
    #                 filter_table_list,
    #                 filter_column_list,
    #                 score_column_list,
    #         )
    #         ht = ht.checkpoint(f'gs://genetics-gym/checkpoints/{analysis_name}.ht', overwrite=True)
    #         print('checkpointed')
    #         results = compute_abs_spearman_grouped_many_scores(ht, score_column_list, outcome_column, group_by)
            
    #         df_results = pd.DataFrame({x:[results[x]] for x in results})
    #         df_results.rename(columns = score_name_mapper, inplace=True)
    #         df_results['method'] = [analysis_name]
    #         df_results.to_csv(f'gs://genetics-gym/results/{analysis_name}_results.tsv', sep='\t', index=False)
