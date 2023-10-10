import pandas as pd
import numpy as np

def DataProcessing(phenotype = "Amylose content", NumberOfSNP = 16000):
    # sativas413 is a genotype table
    sativas413 = pd.read_csv("rice413/sativas413.csv", header= 0)
    sativas413['CHROM'] = sativas413['CHROM'].astype('str')
    sativas413['POS'] = sativas413['POS'].astype('str')
    sativas413["POSITION"] = sativas413['CHROM'] + "_" + sativas413['POS']
    sativas413 = sativas413.set_index("POSITION")
    sativas413 = sativas413.drop(['CHROM', "POS", "REF", "ALT", "QUAL", "FILTER"], axis= 1)

    #get p-value
    # pval_all = pd.read_csv("rice413/MixedModel_Pval_all.txt", sep= "\t")
    # pval_all = pval_all.drop(['CHR', "position.IRGSP.v4", "position.MSU.v6.corrected"], axis= 1)
    #
    # #The names betwenn sativas413 and pval_all is the same.
    # pval_all.rename(columns= {'SNPID': "ID"}, inplace= True)
    #
    # pval_pheno = pval_all[["ID", phenotype]]
    # pval_pheno[phenotype] = -np.log10(pval_pheno[phenotype])
    #
    # # genoPval = pd.concat([sativas413, pval_pheno], axis= 1)
    genoPval = sativas413



    # #Sort by P-Value
    # genoPval.sort_values(axis= 0, by = phenotype, ascending= False, inplace= True)
    # # print(genoPval)
    #
    # #Delete some SNP sites
    # genoPval = genoPval[:NumberOfSNP]
    #
    genotype = genoPval.drop(["ID"], axis=1)
    # genotype.set_index("POSITION")
    #Coding
    genotype = genotype.T
    genotype = genotype.replace('0/0', 0)
    genotype = genotype.replace('0/1', 1)
    genotype = genotype.replace('1/1', 2)
    genotype = genotype.replace('./.', -1)
    # print(genotype)

    # Fusion genotype data and table type data
    phenotypetxt = pd.read_csv("rice413/RiceDiversity_44K_Phenotypes_34traits_PLINK.txt", header=0, sep="\t")
    phenotypetxt.index = phenotypetxt['HybID'] + "_" + phenotypetxt["NSFTVID"].astype(np.str_)
    phenotypetxt = phenotypetxt.drop(['HybID', "NSFTVID", ], axis=1)
    phenotypetxt = phenotypetxt[[phenotype]]

    phenotypetxt = phenotypetxt.dropna(axis=0, how='any')  # drop all rows that have any NaN values
    #
    PhenoGenotype = pd.merge(phenotypetxt, genotype, how="left", left_index=True, right_index=True)



    PhenoGenotype.to_csv("rice413/PhenoGenotype"+ "_" + phenotype + str(NumberOfSNP) +".csv", header= True, sep= ",", index= True) #index = 0不保存行索引
    # print(GenoPhenotype)

#Amylose content    Alkali spreading value     Protein content




# DataProcessing(phenotype = "Plant height", NumberOfSNP = "all")
DataProcessing(phenotype = "Seed number per panicle", NumberOfSNP = "all")
# DataProcessing(phenotype = "Seed surface area", NumberOfSNP = "all")
