import hail as hl
hl.init(worker_memory="highmem", driver_memory='highmem') 

# all are keyed by locus, alleles
collected_paths = {
    'proteinmpnn': {'path': 'gs://trisha-tmp/new_VSM_temp/proteinmpnn_llr_collected.ht', 
             'score_cols': ['proteinmpnn_llr'], 
             'higher_is_less_deleterious' : {'proteinmpnn_llr': True}
            }, 
    'esm1b': {'path': 'gs://trisha-tmp/new_VSM_temp/esm1b_lj2_collected.ht', 
              'score_cols': ['esm_comb'],
              'higher_is_less_deleterious' : {'esm_comb': True}
             },
    'pai3d': {'path': 'gs://trisha-tmp/new_VSM_temp/score_PAI3D_collected.ht', 
              'score_cols': ['score_PAI3D'], 
              'higher_is_less_deleterious' : {'score_PAI3D': False}
            },
    'revel': {'path': 'gs://trisha-tmp/new_VSM_temp/revel_collected.ht',
            'score_cols': ['revel'], 
            'higher_is_less_deleterious' : {'revel': False}
            },
    'rasp': {'path': 'gs://trisha-tmp/new_VSM_temp/rasp_score_collected.ht', 
             'score_cols': ['rasp_score'], 
             'higher_is_less_deleterious' : {'rasp_score': False}
            },
    'AM': {'path': 'gs://trisha-tmp/new_VSM_temp/AM_comb_collected.ht', 
            'score_cols': ['AM_comb'],
            'higher_is_less_deleterious' : {'AM_comb': False}
            },
    'misfit_D': {'path': 'gs://trisha-tmp/new_VSM_temp/MisFit_D_collected.ht', 
            'score_cols': ['MisFit_D'],
            'higher_is_less_deleterious' : {'MisFit_D': False}
            },
    'misfit_S': {'path': 'gs://trisha-tmp/new_VSM_temp/MisFit_S_collected.ht', 
            'score_cols': ['MisFit_S'],
            'higher_is_less_deleterious' : {'MisFit_S': False}
            },
    'mpc': {'path': 'gs://missense-scoring/mpc_grch38_deduped_with_outliers_2024-04-30.ht',
            'score_cols': ['mpc'],
            'higher_is_less_deleterious' : {'mpc': False}
        }
}

joined_paths = {
    # 'proteinmpnn': {'path': 'gs://trisha-tmp/new_VSM_temp/vsm_proteinmpnn.ht', 
    #          'score_cols': ['proteinmpnn_llr'], 
    #          'higher_is_less_deleterious' : {'proteinmpnn_llr': True}
    #         }, 
    # 'esm1b': {'path': 'gs://trisha-tmp/new_VSM_temp/vsm_esm1b.ht', 
    #           'score_cols': ['esm1b'],
    #           'higher_is_less_deleterious' : {'esm1b': True}
    #          },
    # 'pai3d': {'path': 'gs://trisha-tmp/new_VSM_temp/vsm_pai3d.ht', 
    #           'score_cols': ['score_PAI3D'], 
    #           'higher_is_less_deleterious' : {'score_PAI3D': False}
    #         },
    # 'revel': {'path': 'gs://trisha-tmp/new_VSM_temp/vsm_revel.ht',
    #         'score_cols': ['revel'], 
    #         'higher_is_less_deleterious' : {'revel': False}
    #         },
    # 'rasp': {'path': 'gs://trisha-tmp/new_VSM_temp/vsm_rasp.ht', 
    #          'score_cols': ['rasp_score'], 
    #          'higher_is_less_deleterious' : {'rasp_score': False}
    #         },
    # 'AM': {'path': 'gs://trisha-tmp/new_VSM_temp/vsm_am.ht', 
    #         'score_cols': ['AM'],
    #         'higher_is_less_deleterious' : {'AM': False}
    #         },
    # 'misfit': {'path': 'gs://trisha-tmp/new_VSM_temp/vsm_misfit.ht', 
    #         'score_cols': ['MisFit_D', 'MisFit_S'],
    #         'higher_is_less_deleterious' : {'MisFit_D': False, 'MisFit_S': False}
    #         },
    'popeve': {
        'path': 'gs://trisha-tmp/new_VSM_temp/vsm_popeve.ht',
        'score_cols': ['popEVE', 'EVE', 'ESM_1v'],
        'higher_is_less_deleterious' : {'popEVE': True, 'EVE': False, 'ESM_1v': True}
    },
    # 'polyphen': {
    #     'path': 'gs://trisha-tmp/new_VSM_temp/vsm_polyphen.ht',
    #     'score_cols': ['polyphen_score'],
    #     'higher_is_less_deleterious' : {'polyphen_score': False}
    # },
    # 'cpt': {
    #     'path': 'gs://trisha-tmp/new_VSM_temp/vsm_cpt.ht',
    #     'score_cols': ['cpt1_score'],
    #     'higher_is_less_deleterious' : {'cpt1_score': False}
    # }
}

no_transcript_tables = { 
    'mpc': {'path': 'gs://missense-scoring/mpc_grch38_deduped_with_outliers_2024-04-30.ht',
            'score_cols': ['mpc'],
            'higher_is_less_deleterious' : {'mpc': False}
        },
    'cadd': {
        'path': 'gs://genetics-gym-not-public/Trisha/hail_VSM_tables_updated/cadd.ht',
        'score_cols': ['cadd_score'],
        'higher_is_less_deleterious' : {'cadd_score': False}
    },
    'gpn_msa': {
        'path': 'gs://genetics-gym-not-public/Trisha/hail_VSM_tables_updated/gpn_msa.ht', 
        'score_cols': ['gpn_msa_score'],
        'higher_is_less_deleterious' : {'gpn_msa_score': False}
    }
}

## 0. Make grouped Linker table
    # just need to run this the first time
# vsm_ht = hl.read_table('gs://trisha-tmp/new_VSM_temp/vsm_base_per_gene_SNP_2.ht')
# vsm_ht = vsm_ht.key_by('locus', 'alleles', 'ensg')
# vsm_ht = vsm_ht.checkpoint('gs://trisha-tmp/new_VSM_temp/vsm_base_em_gene_k_lag.ht', overwrite=True)
# vsm_ht = hl.read_table('gs://trisha-tmp/new_VSM_temp/vsm_base_em_gene_k_lag.ht')

# other_columns = ['aa_pos', 'aa_ref', 'aa_alt', 'gene_symbol',
#                 'ensp','uniprot_base', 'gene_id_em']
# columns_to_specify = ['enst', 'mane_NM', 'mane_select', 'canonical', 'uniprot_em_id',  'transcript_mane_select']
# vsm_ht = vsm_ht.collect_by_key()              
# for c in other_columns:
#     vsm_ht = vsm_ht.annotate(
#         **{f'{c}':  hl.find(lambda x: hl.is_defined(x[c]) == True, vsm_ht.values)[c]}
#     )
# vsm_ht = vsm_ht.annotate(
#     mane_enst = hl.find(lambda x: x.mane_select == True, vsm_ht.values)['enst'], 
#     canon_enst = hl.find(lambda x: x.canonical == True, vsm_ht.values)['enst'], 
#     **{'all_enst':  hl.filter(lambda x: hl.is_defined(x.enst) == True, vsm_ht.values).map(lambda x: x['enst'])},

#     mane_NM = hl.find(lambda x: x.mane_select == True, vsm_ht.values)['mane_NM'], 
#     mane_select = hl.if_else(hl.is_defined(hl.find(lambda x: x.mane_select == True, vsm_ht.values)['mane_select']), True, False),
#     canonical = hl.if_else(hl.is_defined(hl.find(lambda x: x.canonical == True, vsm_ht.values)['canonical']), True, False),

#     mane_uniprot_iso_em = hl.find(lambda x: x.mane_select == True, vsm_ht.values)['uniprot_em_id'],
#     canon_uniprot_iso_em = hl.find(lambda x: x.canonical == True, vsm_ht.values)['uniprot_em_id'],
#     **{'all_uniprot_iso_em':  hl.filter(lambda x: hl.is_defined(x.uniprot_em_id) == True, vsm_ht.values).map(lambda x: x['uniprot_em_id'])},

#     transcript_mane_select = hl.find(lambda x: x.transcript_mane_select == True, vsm_ht.values)['transcript_mane_select'],

# )
# vsm_ht = vsm_ht.checkpoint('gs://trisha-tmp/new_VSM_temp/vsm_base_per_gene_SNP.ht', overwrite=True)
# vsm_ht = vsm_ht.annotate(
#     mane_or_canon_enst = hl.coalesce(vsm_ht.mane_enst, vsm_ht.canon_enst, hl.is_missing(hl.tstr))
# )
# vsm_ht = vsm_ht.rename({'gene_id_em': 'gene_symbol_em'})
# vsm_ht = vsm_ht.checkpoint('gs://trisha-tmp/new_VSM_temp/vsm_base_per_gene_SNP_2.ht', overwrite=True)

# vsm_ht = hl.read_table('gs://trisha-tmp/new_VSM_temp/vsm_base_per_gene_SNP_2.ht')
# vsm_ht = vsm_ht.drop('values')
# vsm_ht = vsm_ht.checkpoint('gs://trisha-tmp/new_VSM_temp/vsm_base_per_gene_SNP_3.ht')

# 1.ENST based tables - join on locus, allele, enst
# x = 0
# final_str = ''
# for method, info in joined_paths.items(): 
#     print(method)
#     if method == 'popeve': 
#         continue
#     ht = hl.read_table(f'gs://trisha-tmp/new_VSM_temp/per_SNP_gene/max_per_SNP_gene_{method}.ht')
#     final_score_names = []
#     for s in info['score_cols']:
#         if info['higher_is_less_deleterious'][s]:
#             final_score_names.append(f'{s}_neg')
#         else:
#             final_score_names.append(s)
#     ht = ht.select(*final_score_names)
#     vsm_ht = vsm_ht.join(ht, how='left')
#     vsm_ht = vsm_ht.checkpoint(f'gs://trisha-tmp/new_VSM_temp/{x}_table_max_per_SNP_gene_{method}.ht', overwrite=True)
#     final_str = f'gs://trisha-tmp/new_VSM_temp/{x}_table_max_per_SNP_gene_{method}.ht'
#     x += 1

# print('done: ', x)
# print(final_str)
# 2. Refseq based tables - just popeve
x = 9
method = 'popeve'
info = joined_paths[method]
vsm_ht = hl.read_table('gs://trisha-tmp/new_VSM_temp/8_table_max_per_SNP_gene_cpt.ht')
# change keys to locus, alleles, gene_symbol to join with popeve table
vsm_ht = vsm_ht.key_by('locus', 'alleles', 'gene_symbol')
final_score_names = []
ht = hl.read_table(f'gs://trisha-tmp/new_VSM_temp/per_SNP_gene/max_per_SNP_gene_{method}.ht')
for s in info['score_cols']:
    if info['higher_is_less_deleterious'][s]:
        final_score_names.append(f'{s}_neg')
    else:
        final_score_names.append(s)
ht = ht.select(*final_score_names)
vsm_ht = vsm_ht.join(ht, how='left')
vsm_ht = vsm_ht.checkpoint(f'gs://trisha-tmp/new_VSM_temp/{x}_table_max_per_SNP_gene_{method}.ht', overwrite=True)
x += 1

print(x)

# 3. NO transcript tables - join on locus, alleles
# vsm_ht = hl.read_table(f'gs://trisha-tmp/new_VSM_temp/table_max_per_SNP_gene_{method}_{x}.ht')
vsm_ht = vsm_ht.key_by('locus', 'alleles')
for method, info in no_transcript_tables.items(): 
    print(method)
    final_score_names = []
    for s in info['score_cols']:
        if info['higher_is_less_deleterious'][s]:
            final_score_names.append(f'{s}_neg')
        else:
            final_score_names.append(s)
    ht = hl.read_table(info['path'])
    ht = ht.select(*final_score_names)
    vsm_ht = vsm_ht.join(ht, how='left')
    vsm_ht = vsm_ht.checkpoint(f'gs://trisha-tmp/new_VSM_temp/{x}_table_max_per_SNP_gene_{method}.ht', overwrite=True)
    x += 1
    vsm_ht.show()

vsm_ht = vsm_ht.checkpoint('gs://trisha-tmp/new_VSM_temp/vsm_all_tables/vsm_all_SNP_gene_la.ht', overwrite=True)
vsm_ht = vsm_ht.key_by('locus', 'alleles', 'ensg')
vsm_ht = vsm_ht.checkpoint('gs://trisha-tmp/new_VSM_temp/vsm_all_tables/vsm_all_SNP_gene_lag.ht', overwrite=True)

# vsm_ht.show()

# # Format
# vsm_ht = hl.read_table('gs://trisha-tmp/new_VSM_temp/table_max_per_SNP_gene_all_k_lag.ht')
# score_columns = ['proteinmpnn_llr_neg', 'esm1b_neg', 
#             'score_PAI3D', 'revel', 'rasp_score', 'AM', 
#             'MisFit_D', 'MisFit_S', 'polyphen_score', 'cpt1_score', 
#             'popEVE_neg', 'EVE', 'ESM_1v_neg', 'cadd_score', 'gpn_msa_score']
# other_columns = ['aa_pos', 'aa_ref', 'aa_alt', 'gene_symbol',
#                 'ensp','uniprot_base', 'em_gene']
# columns_to_specify = ['enst', 'mane_NM', 'mane_select', 'canonical', 'uniprot_em_id',  'transcript_mane_select']
# vsm_ht = vsm_ht.collect_by_key()              
# for s in score_columns:
#     vsm_ht = vsm_ht.annotate(
#         **{f'{s}':  hl.find(lambda x: hl.is_defined(x[s]) == True, vsm_ht.values)[s]}
#     )
# for c in other_columns:
#     vsm_ht = vsm_ht.annotate(
#         **{f'{c}':  hl.find(lambda x: hl.is_defined(x[c]) == True, vsm_ht.values)[c]}
#     )
# vsm_ht = vsm_ht.annotate(
#     mane_enst = hl.find(lambda x: x.mane_select == True, vsm_ht.values)['enst'], 
#     canon_enst = hl.find(lambda x: x.canonical == True, vsm_ht.values)['enst'], 
#     **{'all_enst':  hl.filter(lambda x: hl.is_defined(x.enst) == True, vsm_ht.values).map(lambda x: x['enst'])},

#     mane_NM = hl.find(lambda x: x.mane_select == True, vsm_ht.values)['mane_NM'], 
#     mane_select = hl.if_else(hl.is_defined(hl.find(lambda x: x.mane_select == True, vsm_ht.values)['mane_select']), True, False),
#     canonical = hl.if_else(hl.is_defined(hl.find(lambda x: x.canonical == True, vsm_ht.values)['canonical']), True, False),

#     mane_uniprot_em_id = hl.find(lambda x: x.mane_select == True, vsm_ht.values)['uniprot_em_id'],
#     canon_uniprot_em_id = hl.find(lambda x: x.canonical == True, vsm_ht.values)['uniprot_em_id'],
#     **{'all_uniprot_em_id':  hl.filter(lambda x: hl.is_defined(x.uniprot_em_id) == True, vsm_ht.values).map(lambda x: x['uniprot_em_id'])},

#     transcript_mane_select = hl.find(lambda x: x.transcript_mane_select == True, vsm_ht.values)['transcript_mane_select'],

# )
# vsm_ht = vsm_ht.checkpoint('gs://genetics-gym-not-public/Trisha/vsm_all_tables/vsm_all_SNP_gene.ht', overwrite=True)

# vsm_ht = vsm_ht.drop('values')
# vsm_ht.write('gs://genetics-gym-not-public/Trisha/vsm_all_tables/vsm_all_SNP_gene_cleaned.ht', overwrite=True)
