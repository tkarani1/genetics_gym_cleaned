import hail as hl
hl.init(worker_memory="highmem")

all_vsm = hl.read_table('gs://genetics-gym-not-public/Trisha/updated_data_tables/new_vsm_all_outer.ht')
all_vsm = all_vsm.drop('aa_pos_1', 'aa_ref_1', 'aa_alt_1', 'gene_symbol_1', 
                       'uniprot_isoform_1', 'mane_select_1', 'canonical_1', 
                      'uniprot_mapped_1', 'uniprot_merged_1')

#all_vsm = all_vsm.head(100)
score_list = ['esm1b', 'rasp', 'MisFit_S', 'MisFit_D', 'AM', 'mpc', 'esm1b_neg', 'popEVE','EVE', 'ESM_1v', 'esm1b_rs', 'esm1b_rs_neg']
#score_list = ['esm1b', 'AM']

# # 1: add negatives
# higher_is_less_deleterious = {
#     # 'esm_score':True,
#     'esm1b': True,
#     'proteinmpnn_llr':True,
#     'am_pathogenicity':False,
#     'score_PAI3D':False,
#     'rasp':False,
#     'MisFit_S':False,
#     'MisFit_D':False,
#     'popEVE':True,
#     'EVE':False,
#     'ESM_1v':True,
#     'mpc':False,
# }

# for s in score_list: 
#     if s in list(higher_is_less_deleterious.keys()):
#         if higher_is_less_deleterious[s]:
#             all_vsm = all_vsm.annotate(
#                 **{f'{s}_neg': 1-all_vsm[s]}
#             )
#             score_list.append(f'{s}_neg')

# all_vsm = all_vsm.checkpoint(f'gs://genetics-gym-not-public/Trisha/updated_data_tables/new_vsm_all_outer_negs.ht', overwrite=True)
# # print(0)
all_vsm = hl.read_table(f'gs://genetics-gym-not-public/Trisha/updated_data_tables/new_vsm_all_outer_negs.ht')


all_vsm = all_vsm.key_by('locus', 'alleles')
all_vsm = all_vsm.checkpoint(f'gs://genetics-gym-not-public/Trisha/updated_data_tables/cp{0}.ht', overwrite=True)
# all_vsm = hl.read_table(f'gs://genetics-gym-not-public/Trisha/updated_data_tables/cp{0}.ht')
all_vsm = all_vsm.collect_by_key()
# print(1)
annotate_dict = {}

score_list = ['esm1b', 'rasp', 'MisFit_S', 'MisFit_D', 'AM', 'mpc', 'esm1b_neg', 'popEVE','EVE', 'ESM_1v', 'esm1b_rs', 'esm1b_rs_neg', 'popEVE_neg', 'ESM_1v_neg']
for s in score_list: 
    annotate_dict[f'{s}_mane'] = hl.find(lambda x: x.mane_select == True, all_vsm.values)[s]
    annotate_dict[f'{s}_canonical'] = hl.find(lambda x: x.canonical == True, all_vsm.values)[s]
    annotate_dict[f'{s}_other'] = all_vsm.values.filter(lambda x: hl.is_defined(x[s])).map(lambda x: x[s])
print(2)
all_vsm = all_vsm.annotate(
    **annotate_dict
)
print(3)
all_vsm = all_vsm.checkpoint(f'gs://genetics-gym-not-public/Trisha/updated_data_tables/cp{1}.ht', overwrite=True)
# print(all_vsm.describe())
print(4)

# all_vsm = hl.read_table(f'gs://genetics-gym-not-public/Trisha/updated_data_tables/cp{1}.ht')

annotate_dict = {}
for s in score_list:
    annotate_dict[f'{s}'] = hl.coalesce(all_vsm[f'{s}_mane'], all_vsm[f'{s}_canonical'], hl.find(lambda x: hl.is_defined(x), all_vsm[f'{s}_other']))

    
print(5)
all_vsm = all_vsm.annotate(
    **annotate_dict
)
print(6)

for s in score_list: 
    all_vsm = all_vsm.drop(f'{s}_mane', f'{s}_canonical', f'{s}_other')
print(7)
all_vsm = all_vsm.checkpoint(f'gs://genetics-gym-not-public/Trisha/updated_data_tables/cp{2}.ht', overwrite=True)
print(8)

    
all_vsm = all_vsm.annotate(
    mane_ENST = hl.find(lambda x: hl.is_defined(x), all_vsm.values.mane_ENST),
    aa_pos = hl.find(lambda x: hl.is_defined(x), all_vsm.values.aa_pos) ,
    aa_ref = hl.find(lambda x: hl.is_defined(x), all_vsm.values.aa_ref),
    aa_alt = hl.find(lambda x: hl.is_defined(x), all_vsm.values.aa_alt), 
    gene_symbol = hl.find(lambda x: hl.is_defined(x), all_vsm.values.aa_ref), 
    uniprot_mapped = hl.find(lambda x: hl.is_defined(x), all_vsm.values.uniprot_mapped)

)
print(9)

all_vsm = all_vsm.drop('values')
print(10)
all_vsm = all_vsm.checkpoint('gs://genetics-gym-not-public/Trisha/updated_data_tables/new_vsm_outer_grouped.ht', overwrite=True)
print(11)
# all_vsm.show()