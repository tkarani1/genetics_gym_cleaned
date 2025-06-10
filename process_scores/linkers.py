# This file was originally from an ipynb 

import hail as hl
hl.init()

def convert_aa_three_to_one(ht):
    # Mapping dictionary from three-letter to one-letter amino acid codes
    aa_dict = {
        'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D',
        'Cys': 'C', 'Gln': 'Q', 'Glu': 'E', 'Gly': 'G',
        'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
        'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S',
        'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
        'Sec': 'U', 'Pyl': 'O', 'Asx': 'B', 'Glx': 'Z',
        'Xaa': 'X', 'Ter': '*'
    }

    # Create a Hail dictionary expression
    aa_mapping = hl.dict(aa_dict)

    # Annotate the table by mapping the three-letter codes to one-letter codes
    ht = ht.annotate(
        aa_ref = aa_mapping.get(ht.aa_ref, ht.aa_ref),
        aa_alt = aa_mapping.get(ht.aa_alt, ht.aa_alt)
    )

    return ht

uni_mapping = hl.import_table('gs://genetics-gym-not-public/Trisha/uniprot_mapping_cleaned.tsv')
uni_mapping_2 = uni_mapping.key_by(uni_mapping.Gene_Name)

ht = hl.read_table('gs://missense-scoring/mutation/everything_raw.ht')
ht_2 = ht.explode(ht.vep.transcript_consequences)


## Missense 
ht_3 = ht_2.annotate(uniprot_isoform = hl.if_else(          # fill NA to avoid drop when exploding
    hl.is_defined(ht_2.vep.transcript_consequences.uniprot_isoform), 
    ht_2.vep.transcript_consequences.uniprot_isoform, ['x']))
ht_4 = ht_3.explode(ht_3.uniprot_isoform)
ht_5 = ht_4.filter(ht_4.vep.transcript_consequences.most_severe_consequence == "missense_variant")

pattern = r".*:p.(\D+)(\d+)(\D+)"
array_match = ht_5.vep.transcript_consequences.hgvsp.first_match_in(pattern)
ht_6 = ht_5.annotate(aa_ref = array_match[0], aa_pos = hl.int(array_match[1]), aa_alt = array_match[2])

# ENST based 
ht_6_ens = ht_6.filter(ht_6.vep.transcript_consequences.transcript_id.startswith('ENST') )
ht_38_mis_ens = ht_6_ens.select(aa_pos = ht_6_ens.aa_pos,
                                       aa_ref = ht_6_ens.aa_ref,
                                       aa_alt = ht_6_ens.aa_alt,
                                        gene_symbol = ht_6_ens.vep.transcript_consequences.gene_symbol,
                                    uniprot_isoform = ht_6_ens.uniprot_isoform,
                                       enst = ht_6_ens.vep.transcript_consequences.transcript_id,
                                       ensp = ht_6_ens.vep.transcript_consequences.protein_id,
                                       ensg = ht_6_ens.vep.transcript_consequences.gene_id, 
                                    mane_NM = ht_6_ens.vep.transcript_consequences.mane_select,
                                   mane_select = hl.is_defined(ht_6_ens.vep.transcript_consequences.mane_select),
                                canonical = hl.is_defined(ht_6_ens.vep.transcript_consequences.canonical))
ht_38_mis_ens = convert_aa_three_to_one(ht_38_mis_ens)


    # fill in missing uniprot name
ht_38_mis_ens_2 = ht_38_mis_ens.annotate(uniprot_mapped = uni_mapping_2[ht_38_mis_ens.gene_symbol].uniprot_id)
ht_38_mis_ens_3 = ht_38_mis_ens_2.annotate(
                        uniprot_merged = hl.if_else(~(ht_38_mis_ens_2.uniprot_isoform == 'x'), 
                                                    ht_38_mis_ens_2.uniprot_isoform, 
                                                   ht_38_mis_ens_2.uniprot_mapped)
                        )

ht_38_mis_ens_3.write('gs://genetics-gym-not-public/Trisha/linker_alt_ens_iso_all_mapped.ht')

# RefSeq based
ht_6_rs = ht_6.filter(ht_6.vep.transcript_consequences.transcript_id.startswith('NM'))
ht_38_mis_rs = ht_6_rs.select(aa_pos = ht_6_rs.aa_pos,
                                       aa_ref = ht_6_rs.aa_ref,
                                       aa_alt = ht_6_rs.aa_alt,
                                        gene_symbol = ht_6_rs.vep.transcript_consequences.gene_symbol,
                                    uniprot_isoform = ht_6_rs.uniprot_isoform,
                                       NM = ht_6_rs.vep.transcript_consequences.transcript_id,
                                       NP = ht_6_rs.vep.transcript_consequences.protein_id,
                                       ncbi = ht_6_rs.vep.transcript_consequences.gene_id,
                                    mane_ENST = ht_6_rs.vep.transcript_consequences.mane_select,
                                   mane_select = hl.is_defined(ht_6_rs.vep.transcript_consequences.mane_select),
                                canonical = hl.is_defined(ht_6_rs.vep.transcript_consequences.canonical))
ht_alt_38 = convert_aa_three_to_one(ht_38_mis_rs)

    # fill in missing uniprot name
ht_alt_38_2 = ht_alt_38.annotate(uniprot_mapped = uni_mapping_2[ht_alt_38.gene_symbol].uniprot_id)
ht_alt_38_3 = ht_alt_38_2.annotate(
                        uniprot_merged = hl.if_else(~(ht_alt_38_2.uniprot_isoform == 'x'), 
                                                    ht_alt_38_2.uniprot_isoform, 
                                                   ht_alt_38_2.uniprot_mapped)
                        )

ht_alt_38_3.write('gs://genetics-gym-not-public/Trisha/linker_alt_refseq_iso_all_mapped.ht')

