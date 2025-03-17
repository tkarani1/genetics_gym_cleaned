import hail as hl
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve, average_precision_score, roc_curve, auc
from hail.utils import hadoop_open
import json

def calc_plot_pr_scipy(df, score_name_dict, is_pos_col, eval_name=None): 
    plt.figure(figsize=(7, 5))

    aps_dict = {}
    i = 0
    for score_name, score_col in score_name_dict.items(): 
        precision, recall, thresholds = precision_recall_curve(df[is_pos_col], df[score_col])
        aps = average_precision_score(df[is_pos_col], df[score_col])
        aps_dict[f'{score_name}_pr'] = [aps]
        plt.plot(recall, precision, label=f'{score_name} ({round(aps, 3)})')
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    if eval_name: 
        plt.title(f'PR for {eval_name[0]} on {eval_name[1]}')
    else: 
         plt.title(f'PR on {is_pos_col}')
    plt.legend()
    path = d['results_path']
    plt.savefig(f"{path}/pr_{eval_name[0]}_{eval_name[1]}.png")   
    plt.show()
    plt.close()
    return aps_dict

def calc_plot_roc_scipy(df, score_name_dict, is_pos_col, eval_name=None): 
    plt.figure(figsize=(7, 5))

    roc_dict = {}
    for score_name, score_col in score_name_dict.items(): 
        fpr, tpr, _ = roc_curve(df[is_pos_col], df[score_col])
        auc_val = auc(fpr, tpr)
        roc_dict[f'{score_name}_roc'] =  [auc_val]
        plt.plot(fpr, tpr, label=f'{score_name} ({round(auc_val, 3)})')
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    if eval_name: 
        plt.title(f'ROC for {eval_name[0]} on {eval_name[1]}')
    else: 
         plt.title(f'ROC on {is_pos_col}')
    plt.legend()
    path = d['results_path']
    plt.savefig(f"{path}roc_{eval_name[0]}_{eval_name[1]}.png")  
    plt.show()  
    plt.close()
    return roc_dict

def filter_scores(df, col_name, bounds): 
    lower, upper = bounds
    if lower is None and upper is None: 
        return df
    if lower: 
        df = df[df[col_name] >=lower]
    if upper: 
        df = df[df[col_name] <= upper]
    print(df)
    return df
    


def main(d): 

    # Load scores
    if d['score_is_ht']: 
        score_ht = hl.read_table(d['score_data_path'])
        if not d['is_gsm']: 
            score_ht = score_ht.key_by(d['score_join_on'])
            score_ht = score_ht.distinct()
            d['is_gsm'] = True
        score_df = score_ht.to_pandas()
    else: 
        score_df = pd.read_csv(d['score_data_path'], delimiter=d['score_delim'])
        if not d['is_gsm']:
            score_df = score_df.drop_duplicates(subset=d['score_join_on'], keep='first') 
    
    if d['score_filtering']: # this is kind of ugly
        score_df = filter_scores(score_df, d['score_filtering']['col_name'], (d['score_filtering']['lower_bound'], d['score_filtering']['upper_bound']))
    
    # Load gene list
    if d['eval_is_ht']: 
        eval_ht = hl.read_table(d['eval_data_path'])
        eval_df = eval_ht.to_pandas()
    else: 
        eval_df = pd.read_csv(d['eval_data_path'], delimiter=d['eval_delim'])
        eval_df['is_pos'] = True

    # Join tables
    score_eval_df = score_df.merge(eval_df, how='left', left_on=d['score_join_on'], right_on=d['eval_join_on'])
    score_eval_df['is_pos'] = score_eval_df['is_pos'].fillna(False)

    pr_dict = calc_plot_pr_scipy(score_eval_df, d['score_name_dict'], 'is_pos', [d['plot_title'], d['eval_name']])
    roc_dict = calc_plot_roc_scipy(score_eval_df, d['score_name_dict'], 'is_pos', [d['plot_title'], d['eval_name']])
    
    # Results df
    results_df = pd.DataFrame(pr_dict)
    results_df = results_df.assign(**roc_dict)
    x = d['plot_title']
    y = d['eval_name']
    path = d['results_path']
    #results_df.to_csv(f'{path}gsm_results_{x}_{y}.tsv', sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--json', type=str)
    args = parser.parse_args()
    optional_commands = ['results_path', 'score_filtering']

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