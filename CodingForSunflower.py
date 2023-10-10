import pandas as pd
import numpy as np

def DataProcessing(phenotype = "Amylose content"):
    genotype = pd.read_csv("sunflower/Annuus.ann_gwas.tranche90_snps_bi_AN50_beagle_AF99_missing20_maf5_id_assigned_rand30k.csv",
                           header= 0,
                           index_col= None)
    genotype = genotype.set_index("ID")
    genotype = genotype.drop(['CHROM', "POS", "REF", "ALT", "QUAL", "FILTER"], axis= 1)

    #Coding
    genotype = genotype.T
    genotype = genotype.replace('0|0', 0)
    genotype = genotype.replace('0|1', 1)
    genotype = genotype.replace('1|0', 1)
    genotype = genotype.replace('1|1', 2)
    genotype = genotype.replace('.|.', -1)
    # print(genotype)

    # Fusion genotype data and table type data
    phenotypetxt = pd.read_csv("sunflower/wild.csv",
                               header= 0,
                               index_col= 0,
                               sep=",")


    phenotypetxt = phenotypetxt[[phenotype]]

    phenotypetxt = phenotypetxt.dropna(axis=0, how='any')  # drop all rows that have any NaN values
    #
    PhenoGenotype = pd.merge(phenotypetxt, genotype, how="inner", left_index=True, right_index=True)



    PhenoGenotype.to_csv("sunflower/PhenoGenotype"+ "_" + phenotype +".csv", header= True, sep= ",", index= True) #index = 0不保存行索引
    # print(GenoPhenotype)

#Amylose content    Alkali spreading value     Protein content
# DataProcessing(phenotype = "Leaf area")
# DataProcessing(phenotype = "Flower head diameter")
DataProcessing(phenotype = "Seed rectangular")