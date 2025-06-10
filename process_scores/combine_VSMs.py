import hail as hl
hl.init()

vsm_ens = {
    'esm1b': {'path': 'gs://genetics-gym-not-public/Trisha/updated_data_tables/esm1b_scores.ht', 
              'score_cols': ['esm1b']
             },
    'misfit': {'path': 'gs://genetics-gym-not-public/Trisha/updated_data_tables/misfit_scores.ht', 
               'score_cols': ['MisFit_D', 'MisFit_S']
              },
    'AM': {'path': 'gs://genetics-gym-not-public/Trisha/hail_VSM_tables_updated/AM_fixed.ht', 
           'score_cols': ['AM']
           }, 
    'RASP': {'path': 'gs://missense-scoring/mutation/rasp_scores.ht', 
             'score_cols': ['rasp_score']
            },
    'MPC': {'path': 'gs://missense-scoring/mpc_grch38_deduped_with_outliers_2024-04-30.ht', 
             'score_cols': ['mpc']
            },
    
}

vsm_rs = {
    'popeve': {'path': 'gs://genetics-gym-not-public/Trisha/updated_data_tables/popeve_2025.ht',  
               'score_cols': ['popEVE', 'EVE', 'ESM_1v']}
}

higher_is_less_deleterious = {
    'esm_score':True,
    'esm1b': True,
    'proteinmpnn_llr':True,
    'AM':False,
    'score_PAI3D':False,
    'rasp_score':False,
    'rasp': False,
    'misfit': False,
    'MisFit_S':False,
    'MisFit_D':False,
    'popeve':True,
    'eve':False,
    'esm1_v':True,
    'mpc':False,
}


print('setup')

### Step 1: ENST based VSM Table
vsm_ht = hl.read_table('gs://genetics-gym-not-public/Trisha/linker_ens_iso_all_mapped.ht')
score_list = []

## Uniprot keys
# esm1b
ht = hl.read_table(vsm_ens['esm1b']['path'])
ht = ht.key_by('uniprot_isoform', 'aa_pos', 'aa_ref', 'aa_alt')
vsm_ht = vsm_ht.annotate(
    esm1b = ht[vsm_ht.uniprot_merged, vsm_ht.aa_pos, vsm_ht.aa_ref, vsm_ht.aa_alt].esm_score
)
score_list.extend(vsm_ens['esm1b']['score_cols'])
print(score_list)

# RASP
ht = hl.read_table(vsm_ens['RASP']['path'])
ht = ht.rename({'uniprot_aa_pos': 'aa_pos', 'uniprot_aa_ref': 'aa_ref', 'uniprot_aa_alt': 'aa_alt'})
ht = ht.key_by('uniprot_id', 'aa_pos', 'aa_ref', 'aa_alt')
joined = ht[vsm_ht.uniprot_merged, vsm_ht.aa_pos, vsm_ht.aa_ref, vsm_ht.aa_alt]
vsm_ht = vsm_ht.annotate(
    rasp = joined.rasp_score
)
score_list.extend(vsm_ens['RASP']['score_cols'])
print(score_list)


## ENST keys
# misfit
ht = hl.read_table(vsm_ens['misfit']['path'])
ht = ht.key_by('enst', 'aa_pos', 'aa_alt')
joined = ht[vsm_ht.enst, vsm_ht.aa_pos, vsm_ht.aa_alt]
vsm_ht = vsm_ht.annotate(
    MisFit_S = joined.MisFit_S, 
    MisFit_D = joined.MisFit_D
)
score_list.extend(vsm_ens['misfit']['score_cols'])
print(score_list)



## Locus, allele keys
# AM
ht = hl.read_table(vsm_ens['AM']['path'])
ht = ht.key_by('locus', 'alleles')
vsm_ht = vsm_ht.annotate(
    AM = ht[vsm_ht.locus, vsm_ht.alleles].AM_score, 
)
score_list.extend(vsm_ens['AM']['score_cols'])
print(score_list)



# MPC
ht = hl.read_table(vsm_ens['MPC']['path'])
ht = ht.key_by('locus', 'alleles')
vsm_ht = vsm_ht.annotate(
    mpc = ht[vsm_ht.locus, vsm_ht.alleles].mpc, 
)
score_list.extend(vsm_ens['MPC']['score_cols'])
print(score_list)
print(vsm_ht.describe())

vsm_ht.write('gs://genetics-gym-not-public/Trisha/updated_data_tables/new_vsm_ens_keys.ht')

print('done')

### Step 2: RefSeq based VSM table
vsm_ht = hl.read_table('gs://genetics-gym-not-public/Trisha/linker_alt_refseq_iso_all_mapped.ht')
ht = hl.read_table(vsm_rs['popeve']['path'])
ht = ht.key_by('NP', 'aa_pos', 'aa_ref', 'aa_alt')
joined = ht[vsm_ht.NP, vsm_ht.aa_pos, vsm_ht.aa_ref, vsm_ht.aa_alt]
vsm_ht = vsm_ht.annotate(
        popEVE = joined.popEVE, 
        EVE = joined.EVE, 
        ESM_1v = joined.ESM_1v
    )
vsm_ht.write('gs://genetics-gym-not-public/Trisha/updated_data_tables/new_vsm_2_v2.ht')


### Step 3: Join to create all locus, alleles, transcripts VSM table

# fix neg
for s in score_list:
    print(s)
    if higher_is_less_deleterious[s]:
        vsm_ht = vsm_ht.annotate(
            **{f'{s}_neg': 1-vsm_ht[s]}
        )
print('finished table')

rs_ht = hl.read_table('gs://genetics-gym-not-public/Trisha/updated_data_tables/new_vsm_ens_keys.ht')
ens_table = hl.read_table('gs://genetics-gym-not-public/Trisha/updated_data_tables/new_vsm_2_v2.ht')
rs_ht = rs_ht.annotate(
   mane_ENST_short = rs_ht.mane_ENST.split("\.")[0]
)
ens_table = ens_table.key_by('locus', 'alleles', 'enst')
rs_ht = rs_ht.key_by('locus', 'alleles', 'mane_ENST_short')
joined = ens_table.join(rs_ht, how='inner')

joined.write('gs://genetics-gym-not-public/Trisha/updated_data_tables/new_vsm_all_inner.ht')

joined = ens_table.join(rs_ht, how='outer')
joined.write('gs://genetics-gym-not-public/Trisha/updated_data_tables/new_vsm_all_outer.ht')

