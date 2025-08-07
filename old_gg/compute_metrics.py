import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, binom
from sklearn.metrics import roc_auc_score, roc_curve

def enrichments_from_grouped_contingency_table(df, pval_cutoff, subset_df_path = None, subset_df_colname = None):
    # to do: generalize to group_by column instead of hard-coding uniprot_id

    assert(np.min(df.high_score_within_gene_n) > 0)
    df = df[df['high_score_within_gene_n'] - df['high_score_within_gene_pos'] > 0]
    df = df[(df.total_pos > 0) & (df.total_n - df.total_pos > 0)]
    df = df.copy()
    
    if subset_df_path is not None:
        df_mapper = pd.read_csv('gs://missense-scoring/uniprotkb_id_mapping.tsv.gz', sep='\t', compression='gzip')
        df_mapper['gene_symbol'] = [str(x).split() for x in df_mapper['Gene Names']]
        df_mapper = df_mapper.explode(column='gene_symbol')
        df_mapper = df_mapper.rename(columns = {'Entry':'uniprot_id', 'Entry Name':'uniprot_name'})
        subset_df = pd.read_csv(subset_df_path, header=None, names=[subset_df_colname])
        df_filter = pd.merge(df_mapper, subset_df, how='inner', on=subset_df_colname)
        df_filter = df_filter[['uniprot_id']].drop_duplicates()
        df = pd.merge(df, df_filter, how='inner', on='uniprot_id')
    
    n_uniprots = len(df['uniprot_id'].drop_duplicates())

    num1 = df['high_score_within_gene_pos']
    den1 = df['total_pos']
    num2 = df['high_score_within_gene_n'] - df['high_score_within_gene_pos']
    den2 = df['total_n'] - df['total_pos']
    df['enrichment'] = (num1 / den1) / (num2 / den2)

    df = df.pivot(index=['uniprot_id', 'cutoff'], columns='score_column', values='enrichment')    
    score_column_list = df.columns
    n_scores = len(score_column_list)

    df = df.reset_index()
    df = df.dropna()
    df = df.drop('uniprot_id', axis=1)

    df_results_all_cutoffs = df.groupby('cutoff').mean()

    df = df[df.cutoff==pval_cutoff]
    df = df.drop('cutoff', axis=1)

    df_results_cutoff = df.mean(axis=0)
    stderror_cutoff = df.std(axis=0)/np.sqrt(n_uniprots)
    df_w_std_error = pd.concat([df_results_cutoff, stderror_cutoff], axis=1)
    df_w_std_error.columns = ['enrichment', 'std_error']
    df_w_std_error['cutoff'] = pval_cutoff
    df_w_std_error = df_w_std_error.reset_index(names=['score_column'])

    pvals = np.ones((n_scores, n_scores))
    for i, method1 in enumerate(score_column_list):
        for j, method2 in enumerate(score_column_list):
            x = np.sum(df[method1] > df[method2])
            n = np.sum(~(df[method1]==df[method2]))
            null_p = 0.5
            if n == 0:
                pvals[i,j] = 1
            else:
                pvals[i,j] = binom.sf(x, n, null_p)
            print(method1, method2, x, n, null_p, pvals[i,j])
    df_p = pd.DataFrame(pvals, index=score_column_list, columns=score_column_list)
        
    return df_results_all_cutoffs, df_w_std_error, df_p

def pvals_from_contingency_table(df, cutoff, f_out):
    score_column_list = np.unique(df.score_column)
    df = df[df.cutoff == cutoff]
    
    n_methods = len(score_column_list)
    oddsratios = np.ones((n_methods, n_methods))
    pvals = 0.5 * np.ones((n_methods, n_methods))

    for i, method1 in enumerate(score_column_list):
        for j, method2 in enumerate(score_column_list[0:i]):
            pos_1 = df.loc[df.score_column == method1, 'high_score_pos'].values[0]
            pos_2 = df.loc[df.score_column == method2, 'high_score_pos'].values[0]
            total_1 = df.loc[df.score_column == method1, 'high_score_n'].values[0]
            total_2 = df.loc[df.score_column == method2, 'high_score_n'].values[0]
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
    df_p.to_csv(f_out, sep='\t', index=False)

def results_from_full_table_one_cutoff(df, is_pos_name, score_column_list, cutoff, stat, ncase=None, ncontrol=None):
    # df is not filtered to all defined
    
    n_methods = len(score_column_list)
    oddsratios = np.ones((n_methods, n_methods))
    pvals = 0.5 * np.ones((n_methods, n_methods))

    for i, method1 in enumerate(score_column_list):
        for j, method2 in enumerate(score_column_list[0:i]):
            df_f = df[[is_pos_name, method1, method2]].dropna()
            df_f['method_1_cutoff'] = df_f[method1] >= np.quantile(df_f[method1], cutoff)
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
    elif stat=='rate_ratio':
        C = ncase
        D = ncontrol
    anchor_results = (A / (A+C)) / (B / (B+D))
    results = anchor_results * anchor_oddsratios
    df_results = pd.DataFrame(data=results.reshape((1, len(score_column_list))), columns=score_column_list)

    return df_results, df_p

def results_from_full_table_bootstrap(df, is_pos_name, score_column_list, cutoff, stat, ncase, ncontrol, n_bootstrap_samples):
    to_concat = []
    for i in range(n_bootstrap_samples):
        print(f'bootstrap sample {i+1} out of {n_bootstrap_samples}')
        ii = np.random.choice(range(len(df)), len(df))
        df_sample = df.loc[ii].reset_index()
        df_results, _ = results_from_full_table_one_cutoff(df_sample, is_pos_name, score_column_list, cutoff, stat, ncase, ncontrol)
        to_concat.append(df_results)
    df_results_all = pd.concat(to_concat)
    df_results_no_sampling, _ = results_from_full_table_one_cutoff(df, is_pos_name, score_column_list, cutoff, stat, ncase, ncontrol)
    df_bootstrap = pd.concat([df_results_no_sampling.mean(axis=0), df_results_all.std(axis=0)], axis=1)
    df_bootstrap.columns = [stat, 'std_error']
    df_bootstrap['cutoff'] = cutoff
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

def compare_aucs(df, is_pos_name, score_column_list):
    df = df[[is_pos_name] + score_column_list].dropna()
    df[is_pos_name] = df[is_pos_name].astype(int)
    aucs = np.zeros(len(score_column_list))
    fprs = []
    tprs = []
    score_column = []
    for i,c in enumerate(score_column_list):
        aucs[i] = roc_auc_score(df[is_pos_name], df[c])
        fpr, tpr, _  = roc_curve(df[is_pos_name], df[c])
        fprs.extend(fpr)
        tprs.extend(tpr)
        score_column.extend([c]*len(fpr))
    df_auc = pd.DataFrame({'score_column':score_column_list, 'auc':aucs})
    df_roc = pd.DataFrame({'score_column':score_column, 'fpr':fprs, 'tpr':tprs})
    return df_auc, df_roc