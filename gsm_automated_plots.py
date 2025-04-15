import hail as hl
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve, average_precision_score, roc_curve, auc
from hail.utils import hadoop_open
import json
from scipy.stats import fisher_exact, binom


def results_from_full_table_one_cutoff(df, is_pos_name, score_column_list, cutoff, stat, ncase=None, ncontrol=None):
    # df is not filtered to all defined
    
    n_methods = len(score_column_list)
    oddsratios = np.ones((n_methods, n_methods))
    pvals = 0.5 * np.ones((n_methods, n_methods))

    for i, method1 in enumerate(score_column_list): # odds ratio of comparing one method to another
        for j, method2 in enumerate(score_column_list[0:i]):
            df_f = df[[is_pos_name, method1, method2]].dropna()
            df_f['method_1_cutoff'] = df_f[method1] >= np.quantile(df_f[method1], cutoff)   # using a different set so need to recompute the percentiles
            df_f['method_2_cutoff'] = df_f[method2] >= np.quantile(df_f[method2], cutoff)
            pos_1 = np.sum(df_f[is_pos_name] & df_f['method_1_cutoff'])
            pos_2 = np.sum(df_f[is_pos_name] & df_f['method_2_cutoff'])
            total_1 = np.sum(df_f['method_1_cutoff'])
            total_2 = np.sum(df_f['method_2_cutoff'])
            contingency_table = [
                [pos_1, pos_2],
                [total_1 - pos_1, total_2 - pos_2],
            ]
            oddsratio, pvalue = fisher_exact(contingency_table)
            oddsratios[i,j] = oddsratio
            oddsratios[j,i] = 1/oddsratio
            if oddsratio > 1:
                pvals[i,j] = pvalue/2
                pvals[j,i] = 1 - pvalue/2
            else:
                pvals[i,j] = 1 - pvalue/2
                pvals[j,i] = pvalue/2

    df_p = pd.DataFrame(pvals, index=score_column_list, columns=score_column_list)


    anchor_index =  np.argmax(np.sum(~df[score_column_list].isna(), axis=0))
    anchor_method = score_column_list[anchor_index]
    anchor_oddsratios = oddsratios[:,anchor_index]
    df_f = df[[is_pos_name, anchor_method]].dropna()
    df_f['anchor_cutoff'] = df_f[anchor_method] > np.quantile(df_f[anchor_method], cutoff)
 
    A = np.sum(df_f[is_pos_name] & df_f['anchor_cutoff'])
    B = np.sum((~df_f[is_pos_name]) & df_f['anchor_cutoff'])
    if stat=='enrichment':
        C = np.sum(df_f[is_pos_name] & ~df_f['anchor_cutoff'])
        D = np.sum(~df_f[is_pos_name] & ~df_f['anchor_cutoff'])
        anchor_results = (A / (A+C)) / (B / (B+D))
    elif stat=='rate_ratio':
        C = ncase
        D = ncontrol
        anchor_results = (A / C) / (B / D) 
    results = anchor_results * anchor_oddsratios
    df_results = pd.DataFrame(data=results.reshape((1, len(score_column_list))), columns=score_column_list)

    return df_results, df_p, anchor_method

def results_from_full_table_bootstrap(df, is_pos_name, score_column_list, cutoff, stat, ncase=None, ncontrol=None, n_bootstrap_samples=100 ):
    to_concat = []
    for i in range(n_bootstrap_samples):
        ii = np.random.choice(range(len(df)), len(df))
        df_sample = df.loc[ii].reset_index()
        df_results, _ , anchor_method = results_from_full_table_one_cutoff(df_sample, is_pos_name, score_column_list, cutoff, stat, ncase, ncontrol)    # is this getting called twice (in this function and in the full_table_many_cutoffs() function) # probably not...
        to_concat.append(df_results)
    
    df_results_all = pd.concat(to_concat)

    df_results_no_sampling, _, _ = results_from_full_table_one_cutoff(df, is_pos_name, score_column_list, cutoff, stat, ncase, ncontrol)

    df_bootstrap = pd.concat([df_results_no_sampling.mean(axis=0), df_results_all.std(axis=0)], axis=1)
    df_bootstrap.columns = [stat, 'std_error']
    df_bootstrap['cutoff'] = cutoff
    df_bootstrap['anchor_method'] = anchor_method
    df_bootstrap = df_bootstrap.reset_index(names=['score_column'])
    return df_bootstrap

def results_from_full_table_many_cutoffs(df, is_pos_name, score_column_list, cutoff_list, pval_cutoff, stat, ncase, ncontrol, n_bootstrap_samples=100):
    to_concat_bootstrap = []
    for cutoff in cutoff_list:
        df_bootstrap = results_from_full_table_bootstrap(df, is_pos_name, score_column_list, cutoff, stat, ncase, ncontrol, n_bootstrap_samples)
        to_concat_bootstrap.append(df_bootstrap)
        if cutoff == pval_cutoff:
            _, df_p = results_from_full_table_one_cutoff(df, is_pos_name, score_column_list, cutoff, stat, ncase, ncontrol)
    df_bootstrap_all_cutoffs = pd.concat(to_concat_bootstrap)
    return df_bootstrap_all_cutoffs, df_p



def calc_plot_pr_scipy(df, score_name_dict, is_pos_col, path, eval_name=None): 
    plt.figure(figsize=(7, 5))

    aps_dict = {}
    aps_scaled = {}
    # random_pr = np.sum(df[is_pos_col]) / len(df)
    aps_results = np.ones((len(score_name_dict)))

    for i, (score_name, score_col) in enumerate(score_name_dict.items()): 
        precision, recall, thresholds = precision_recall_curve(df[is_pos_col], df[score_col])
        aps = average_precision_score(df[is_pos_col], df[score_col])
        aps_dict[f'{score_name}_pr'] = [aps]
        aps_results[i] = aps
        plt.plot(recall, precision, label=f'{score_name} ({round(aps, 3)})')
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    if eval_name: 
        plt.title(f'PR for {eval_name[0]} on {eval_name[1]}')
    else: 
         plt.title(f'PR on {is_pos_col}')
    plt.legend()
    plt.savefig(f"{path}/pr_{eval_name[0]}_{eval_name[1]}.png")   
    #plt.show()
    df_results = pd.DataFrame(data=aps_results, index=list(score_name_dict.keys()), columns=['aps']) 

    plt.close()
    
    return df_results, ['aps']

def calc_plot_roc_scipy(df, score_name_dict, is_pos_col, path, eval_name=None): 
    plt.figure(figsize=(7, 5))

    roc_dict = {}
    roc_results = np.ones((len(score_name_dict)))
    for i, (score_name, score_col) in enumerate(score_name_dict.items()): 
        fpr, tpr, _ = roc_curve(df[is_pos_col], df[score_col])
        auc_val = auc(fpr, tpr)
        roc_dict[f'{score_name}_roc'] =  [auc_val]
        roc_results[i] = auc_val
        plt.plot(fpr, tpr, label=f'{score_name} ({round(auc_val, 3)})')
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    if eval_name: 
        plt.title(f'ROC for {eval_name[0]} on {eval_name[1]}')
    else: 
         plt.title(f'ROC on {is_pos_col}')
    plt.legend()
    plt.savefig(f"{path}/roc_{eval_name[0]}_{eval_name[1]}.png")   
    #plt.show()  
    plt.close()
    df_results = pd.DataFrame(data=roc_results, index=list(score_name_dict.keys()), columns=['auc']) 

    return df_results, ['auc']

def calc_obs_exp(df, d_info, score_name_dict, is_pos_col, cutoffs):
    df_copy = df.copy()
    df_copy = df_copy[df_copy[is_pos_col]]
    n_scores = len(score_name_dict)
    n_cutoffs = len(cutoffs)
    results = np.ones((n_scores, n_cutoffs))
    
    for i, c in enumerate(cutoffs):
        for j, (name, col) in enumerate(score_name_dict.items()): 
            df_copy['above_cutoff'] = (df_copy[col] >= np.quantile(df_copy[col], c))
            obs_sum = df_copy.loc[df_copy["above_cutoff"], d_info['observed']].sum()
            exp_sum = df_copy.loc[df_copy["above_cutoff"], d_info['expected']].sum()
            ratio = obs_sum / exp_sum
            results[j, i] = ratio
    cols = [f'oe_{c}' for c in cutoffs]
    df_results = pd.DataFrame(data=results, index=list(score_name_dict.keys()), columns=cols) 
    return df_results, cols

# def filter_scores(df, col_name, bounds): 
#     lower, upper = bounds
#     if lower is None and upper is None: 
#         return df
#     if lower: 
#         df = df[df[col_name] >=lower]
#     if upper: 
#         df = df[df[col_name] <= upper]
#     print(df)
#     return df
    


def main(d): 
    print('hi from automated')
    stat_cols = []

    # Load scores
    if d['score_is_ht']: 
        score_ht = hl.read_table(d['score_data_path'])
        if not d['score_is_gsm']: 
            score_ht = score_ht.key_by(d['score_join_on'])
            score_ht = score_ht.distinct()
            d['score_is_gsm'] = True
        score_df = score_ht.to_pandas()
    else: 
        compression = 'gzip' if d['score_data_path'].endswith('bgz') else None   # add in or condition
        score_df = pd.read_csv(d['score_data_path'], delimiter=d['score_delim'])
        if not d['score_is_gsm']:
            score_df = score_df.drop_duplicates(subset=d['score_join_on'], keep='first') 

    print('score df')
    print(len(score_df))
    print(score_df.isna().sum())
    score_df = score_df.dropna()
    # Load gene list
    if d['eval_is_ht']: 
        eval_ht = hl.read_table(d['eval_data_path'])
        eval_df = eval_ht.to_pandas()
    else: 
        compression = 'gzip' if d['eval_data_path'].endswith('bgz') else None   # add in or condition
        eval_df = pd.read_csv(d['eval_data_path'], delimiter=d['eval_delim'], compression=compression)
        eval_df['is_pos'] = True
    print('eval len')
    print(len(eval_df))

    if 'score_filter' in d and d['score_filter']: 
        score_df = score_df.query(d['score_filter'])
    if 'eval_filter' in d and d['eval_filter']: 
        eval_df = eval_df.query(d['eval_filter'])

    # Join tables
    score_eval_df = score_df.merge(eval_df, how='left', left_on=d['score_join_on'], right_on=d['eval_join_on'])
    score_eval_df['is_pos'] = score_eval_df['is_pos'].fillna(False)
    print('score eval len')
    print(len(score_eval_df))

    original_eval_genes = set(eval_df[d['eval_join_on']])
    #print('eval genes: ', len(original_eval_genes))

    score_eval_df_genes = set(score_eval_df[d['score_join_on']])

    #print('orig diff ', len(original_eval_genes - score_eval_df_genes))

    # GSM on GLE
    pr_results, c = calc_plot_pr_scipy(score_eval_df, d['score_name_dict'], 'is_pos', d['results_path'], [d['plot_title'], d['eval_name']])
    stat_cols = stat_cols + c
    roc_results, c = calc_plot_roc_scipy(score_eval_df, d['score_name_dict'], 'is_pos', d['results_path'], [d['plot_title'], d['eval_name']])
    stat_cols = stat_cols + c

    all_results = pd.concat([pr_results, roc_results], axis=1, join='inner')
    if d['obs/exp']: 
        oe_results, c = calc_obs_exp(score_eval_df, d['obs/exp'], d['score_name_dict'], 'is_pos', d['cutoffs'])
        all_results = pd.concat([all_results, oe_results], axis=1, join='inner')
        stat_cols = stat_cols + c
    x = d['plot_title']
    y = d['eval_name']
    path = d['results_path']
    all_results_sorted = all_results.sort_values(by='aps', ascending=False)
    all_results_sorted.to_csv(f'{path}gsm_eval_{x}_{y}.tsv', sep='\t', index=True)
    return all_results, stat_cols
    # GSM on GLE        # gene-level enrichment (is pos or not)
    # enr_dict = results_from_full_table_many_cutoffs(score_eval_df, 'is_pos', list(d['score_name_dict'].values()), [0.9, 0.95], 0.95, stat='enrichment', ncase=None, ncontrol=None, n_bootstrap_samples=100 )
    # rr_dict = results_from_full_table_many_cutoffs(score_eval_df, 'is_pos', list(d['score_name_dict'].values()), [0.9, 0.95], 0.95, stat='rate_ratio', ncase=None, ncontrol=None, n_bootstrap_samples=100 )
    # print(enr_dict)
    # print(rr_dict)

    # # Results df
    # scipy_df = pd.DataFrame(pr_dict)
    # scipy_df = scipy_df.assign(**roc_dict)
    # print(scipy_df)
    # x = d['plot_title']
    # y = d['eval_name']
    # path = d['results_path']
    # #scipy_df.to_csv(f'{path}gsm_pr_roc_{x}_{y}.tsv', sep='\t', index=False)
    # #enr_dict.to_csv(f'{path}gsm_enr_{x}_{y}.tsv', sep='\t', index=False)
    # #rr_dict.to_csv(f'{path}gsm_rr_{x}_{y}.tsv', sep='\t', index=False)




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--json', type=str)
    args = parser.parse_args()
    optional_commands = ['results_path', 'obs/exp']

    for json_file in args.json.split(','):
        print(json_file)
        print('reading', json_file)
        d = json.load(hadoop_open(json_file))
        print('read')

        d['analysis_name'] = json_file.split('/')[-1][0:-5]
        
        if d['score_delim'] == None: 
            d['score_delim'] = '\t'
        if d['eval_delim'] == None: 
            d['eval_delim'] = '\t'
        
        for c in optional_commands:
            if c not in list(d.keys()): 
                d[c] = None

        if d['results_path'] == None: 
            d['results_path'] = './results/'
        
        print(d)

        main(d)