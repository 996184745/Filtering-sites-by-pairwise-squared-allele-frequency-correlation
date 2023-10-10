import pandas as pd
import numpy as np

def DataProcessing(phenotype = "Amylose content"):
    genotype = pd.read_csv("maize282/SNP55K_maize282_AGP3_20190419.csv", header= 0, index_col= None)
    genotype['CHROM'] = genotype['CHROM'].astype('str')
    genotype['POS'] = genotype['POS'].astype('str')
    genotype["POSITION"] = genotype['CHROM'] + "_" + genotype['POS']
    genotype = genotype.set_index("POSITION")
    genotype = genotype.drop(['CHROM', "POS", "ID", "REF", "ALT", "QUAL", "FILTER"], axis= 1)

    #Coding
    genotype = genotype.T
    genotype = genotype.replace('0/0', 0)
    genotype = genotype.replace('0/1', 1)
    genotype = genotype.replace('1/0', 1)
    genotype = genotype.replace('1/1', 2)
    genotype = genotype.replace('./.', -1)
    # print(genotype)

    # Fusion genotype data and table type data
    phenotypetxt = pd.read_csv("maize282/traitMatrix_maize282NAM_v15-130212.txt",
                               header= None,
                               index_col= 0,
                               sep="\t")

    header_series = phenotypetxt.apply(lambda x: x[0] + '_' + x[1], axis=0)
    phenotypetxt.columns = header_series
    phenotypetxt = phenotypetxt.drop(["<Trait>", "<Header name=env>"])

    phenotypetxt = phenotypetxt[[phenotype]]

    phenotypetxt = phenotypetxt.dropna(axis=0, how='any')  # drop all rows that have any NaN values
    #
    PhenoGenotype = pd.merge(phenotypetxt, genotype, how="inner", left_index=True, right_index=True)



    PhenoGenotype.to_csv("maize282/PhenoGenotype"+ "_" + phenotype +".csv", header= True, sep= ",", index= True) #index = 0不保存行索引
    # print(GenoPhenotype)

#DaystoSilk PlantHeight
#065 06CL1 06PR 07A 07CL1 26M3
DataProcessing(phenotype = "PlantHeight_26M3")