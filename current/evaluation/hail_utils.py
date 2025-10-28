import hail as hl

def join_tables(ht_to, ht_from, columns_to_join=None, filter=False, broadcast=False):
    # keys of ht_from must be columns in ht_to
    # columns_to_join are the columns from ht_from to add to ht_to
    if columns_to_join is None:
        columns_to_join = ht_from.row
    columns_to_join = [c for c in columns_to_join if c in ht_from.row and not c in ht_to.row]  
    if broadcast:
        ht_temp = create_broadcast_dict(
            ht_from.select(*columns_to_join)
        ).get(ht_to.row.select(*list(ht_from.key)))
    else:
        ht_temp = ht_from[ht_to.row.select(*list(ht_from.key))]
    ht_joined = ht_to.annotate(**{c: ht_temp[c] for c in columns_to_join})
    if filter:
        ht_joined = ht_joined.filter(hl.all([hl.is_defined(ht_joined[c]) for c in columns_to_join]))
    return ht_joined

def create_broadcast_dict(key, value = None):
    """
    Create broadcast join (local dictionary from key -> value)
    from a Hail Table.

    :param Expression key: Key Expression
    :param Expression value: Value Expression
    :return: Hail DictExpression (without an index)
    :rtype: DictExpression
    """
    if isinstance(key, hl.Table):
        key = key.key
    ht = key._indices.source
    if value is None:
        value = ht.row_value
    if ht != value._indices.source:
        from hail.expr.expressions import ExpressionException
        raise ExpressionException('create_broadcast_dict: key and value are from different sources')
    return hl.dict(ht.aggregate(hl.agg.collect((key, value)), _localize=False))


def checkpoint_randomfile(ht):
    randomfile = f'gs://missense-scoring/checkpoints/temp_{hl.eval(hl.rand_unif(0, 1))}.ht'
    ht = ht.checkpoint(randomfile)
    return ht

rg37 = hl.get_reference('GRCh37')  
rg38 = hl.get_reference('GRCh38')  
rg37.add_liftover('gs://hail-common/references/grch37_to_grch38.over.chain.gz', rg38)  
def liftover_37_38(ht):
    ht = ht.annotate(
        new_locus=hl.liftover(ht.locus, 'GRCh38', include_strand=True),
        old_locus=ht.locus
    )
    ht = ht.filter(hl.is_defined(ht.new_locus) & ~ht.new_locus.is_negative_strand)  
    ht = ht.key_by()
    ht = ht.annotate(locus=ht.new_locus.result)
    ht = ht.drop('new_locus', 'old_locus')
    return ht

def add_percentiles(ht, score_list):
    print('annotate globals')
    ht = ht.annotate_globals(
        cdfs=ht.aggregate(
            hl.struct(
                **{colname: hl.agg.approx_cdf(ht[colname], k=5000) for colname in score_list}
            ), _localize=False
    ))

    print('checkpointing')
    ht = ht.checkpoint('gs://missense-scoring/temp_checkpoint_2.ht', overwrite=True)
    print('done')

    num_rows = ht.count()
    ht = ht.annotate(**{
        colname+'_percentile':
            ht.cdfs[colname]['ranks'][
                hl.binary_search(ht.cdfs[colname]['values'], ht[colname])
            ] / num_rows for colname in score_list
        }
    )

    print('checkpointing')
    ht = ht.checkpoint('gs://missense-scoring/temp_checkpoint_3.ht', overwrite=True)
    print('done')
    return ht

# def add_quantile_column(ht, colname, reverse_order, k=4000):
#     percentile_name = colname+'_percentile'
#     if reverse_order:
#         colname_neg = colname+'_neg'
#         ht = ht.annotate(**{colname_neg: -ht[colname]})
#         ht = ht.annotate_globals(cdf=ht.aggregate(
#             hl.agg.approx_cdf(ht[colname_neg], k=4000), _localize=False
#         ))
#         ht = ht.annotate(
#             approx_rank=ht.cdf['ranks'][
#                 hl.binary_search(ht.cdf['values'], ht[colname_neg])
#             ]
#         )
#         ht = ht.drop(colname_neg)
#     else:
#         ht = ht.annotate_globals(cdf=ht.aggregate(
#             hl.agg.approx_cdf(ht[colname], k), _localize=False
#         ))
#         ht = ht.annotate(
#             approx_rank=ht.cdf['ranks'][
#                 hl.binary_search(ht.cdf['values'], ht[colname])
#             ]
#         )
#     ht = ht.transmute(**{percentile_name: ht.approx_rank / ht.count()})
#     return ht
