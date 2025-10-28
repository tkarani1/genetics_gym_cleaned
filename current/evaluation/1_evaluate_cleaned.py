import hail as hl
from utils.hail_utils import join_tables
import argparse
from hail.utils import hadoop_open
from analyses_for_manuscript.scores_to_use import PERCENTILE_SCORES, RAW_SCORES, RAW_SCORES_NO_PAI3D, PERCENTILE_SCORES_GENE_MEDIAN, SCORE_LIST_MAP
import json

ht_linker_path = 'gs://genetics-gym/linkers-and-annotations/linker_no_mane_filter.ht'


def join_eval(
        ht_evaluation_path,
        ht_score_path,
        ht_linker_path, 
        is_pos_name,
        filter_table_list,  
        filter_column_list,
        score_column_list,
        intersect_scores = True,
):
    ht_linker = hl.read_table(ht_linker_path)
    ht_evaluation = hl.read_table(ht_evaluation_path, _n_partitions=200)
    ht_score = hl.read_table(ht_score_path, _n_partitions=200)
    ht = join_tables(ht_linker, ht_evaluation, [is_pos_name] + filter_column_list, filter=True)
    ht = ht.filter(hl.is_defined(ht[is_pos_name]))
    ht = join_tables(ht, ht_score, score_column_list, filter=intersect_scores) 
    
    for filter_table in filter_table_list:
        ht_filter = hl.read_table(filter_table)
        ht = join_tables(ht, ht_filter, filter_column_list, filter=True, broadcast=True)
            
    for c in filter_column_list:
        ht = ht.filter(ht[c])

    return ht

# Konrad wrote this

def compute_contingency_table(ht, score_columns, cutoffs, is_pos_name, group_by_column=None):
    if group_by_column:     # style - put all of the group_by together
        output = ht.aggregate(hl.struct(**{
            f'{score_column}_cutoffs': hl.agg.group_by(ht[group_by_column], hl.agg.approx_quantiles(ht[score_column], cutoffs, 1000))
            for score_column in score_columns   # approx quantiles 
        }), _localize=False)
        ht = ht.annotate_globals(**output).checkpoint(hl.utils.new_temp_file())
        grouped_ht = ht.group_by(ht[group_by_column])
    else:
        grouped_ht = ht.group_by(all_genes=True)    # WHAT IS THIS DOING??
    res = {
        'total_pos': hl.agg.count_where(ht[is_pos_name]),
        'total_n': hl.agg.count()
    } 
    # creating the function for calculating the contingency table
    if group_by_column:
        res['cutoff_high_score_pos_n'] = [(score_column,
                                           [(cutoff,
                                             hl.agg.count_where(ht[is_pos_name] & (ht[score_column] >= cutoff)),
                                             hl.agg.count_where(ht[score_column] >= cutoff),
                                             hl.agg.count_where(ht[is_pos_name] & (ht[score_column] >= ht[f'{score_column}_cutoffs'][ht[group_by_column]][i])),
                                             hl.agg.count_where(ht[score_column] >= ht[f'{score_column}_cutoffs'][ht[group_by_column]][i]))
                                            for i, cutoff in enumerate(cutoffs)])
                                          for score_column in score_columns]
    else:
        res['cutoff_high_score_pos_n'] = [(score_column,
                                           [(cutoff,
                                             hl.agg.count_where(ht[is_pos_name] & (ht[score_column] >= cutoff)),
                                             hl.agg.count_where(ht[score_column] >= cutoff))
                                            for cutoff in cutoffs])
                                          for score_column in score_columns]
    
    # apply the function by using aggregate on the ht and ** to unpack dict
    grouped_table = grouped_ht.aggregate(**res).checkpoint(hl.utils.new_temp_file())
    grouped_table = grouped_table.filter(grouped_table.total_n > 0)
    grouped_table = grouped_table.explode(grouped_table.cutoff_high_score_pos_n)
    grouped_table = grouped_table.annotate(
        score_column=grouped_table.cutoff_high_score_pos_n[0],
        cutoff_high_score_pos_n=grouped_table.cutoff_high_score_pos_n[1]
    )
    grouped_table = grouped_table.explode(grouped_table.cutoff_high_score_pos_n)
    res = {'cutoff': grouped_table.cutoff_high_score_pos_n[0],
           'high_score_pos': grouped_table.cutoff_high_score_pos_n[1],
           'high_score_n': grouped_table.cutoff_high_score_pos_n[2]}
    if group_by_column:
        res['high_score_within_gene_pos'] = grouped_table.cutoff_high_score_pos_n[3]
        res['high_score_within_gene_n'] = grouped_table.cutoff_high_score_pos_n[4]
    grouped_table = grouped_table.transmute(**res)
    return grouped_table

def main(
        analysis_name,
        ht_evaluation_path,
        ht_score_path,
        is_pos_name,
        filter_table_list,
        filter_column_list,
        score_column_list,
        cutoff_list,
        group_by,
        make_contingency,
        use_checkpoint,
        file_path
):
    # Step 1: JOIN scores, linker, eval
    if use_checkpoint: # load from checkpoint provided
        fname = file_path + f'hail_evaluation_tables/{analysis_name}.ht'
        ht = hl.read_table(fname, _n_partitions=2000)
    else: # compute joined table
        merged_ht = join_eval(ht_evaluation_path, 
                            ht_score_path, 
                            ht_linker_path,
                            is_pos_name,
                            filter_table_list,
                            filter_column_list, 
                            score_column_list, 
                            intersect_scores = make_contingency)
        ht_path = file_path + f'hail_evaluation_tables/{analysis_name}'
        ht.write(ht_path + '.ht', overwrite=True)
        ht.export(ht_path + '.tsv')

    # Step 2: compute contingency OR write joined table 
    if make_contingency: 
        results = compute_contingency_table(
                ht,
                score_column_list,
                cutoff_list,
                is_pos_name,
                group_by
            )
        df_path = file_path + f'hail_evaluation_tables/{analysis_name}_contingency_tables.tsv'
        results.export(df_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--json', type=str)
    parser.add_argument('--use-checkpoint', type=str)
    args = parser.parse_args()

    for json_file in args.json.split(','):
        print('reading', json_file)
        d = json.load(hadoop_open(json_file))
        print('read')

        d['analysis_name'] = json_file.split('/')[-1][0:-5]
        if d['score_column_list'] in SCORE_LIST_MAP.keys(): 
            t = d['score_column_list']
            d['score_column_list'] = SCORE_LIST_MAP[t]
        else:
            d['score_column_list'] = d['score_column_list'].split(',')
        if not 'group_by' in d:
            d['group_by'] = None
        if not 'cutoff_list' in d:
            d['cutoff_list'] = None
        
        d['file_path'] = "gs://trisha-tmp/recreating_workflow/"
        print('arguments:')
        print(json.dumps(d, indent=4))

        main(
            d['analysis_name'],
            d['ht_evaluation_path'],
            d['ht_score_path'],
            d['is_pos_name'],
            d['filter_table_list'],
            d['filter_column_list'],
            d['score_column_list'],
            d['cutoff_list'],
            d['group_by'],
            d['hail_contingency'],
            args.use_checkpoint,
            d['file_path']
        )