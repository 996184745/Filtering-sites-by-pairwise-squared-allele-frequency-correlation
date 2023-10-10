import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
import matplotlib.pyplot as plt
import matplotlib
from sklearn.inspection import permutation_importance
import allel
from scipy.spatial.distance import squareform
import time
import multiprocessing


def LDHeatmap(pheno, data_file):
    # When r2 = 1, it means that the linkage is completely unbalanced and there is no recombination.
    #
    # When r2 = 0, it means that the linkage is completely balanced and the random group

    # pheno = "Plant height"
    # data_file = "rice413/PhenoGenotype_" + pheno + "all.csv"
    # data_file = "rice413/PhenoGenotype_" + pheno + "_filter0.9.csv"
    # pheno = "Seed number per panicle"
    # pheno = "Seed surface area"
    # data_file = "rice413/PhenoGenotype_Seed surface areaall.csv"


    # pheno = "EarDia"
    # data_file = "maize/PhenoGenotype_" + pheno + ".csv"


    # pheno = "tkw"
    # data_file = "wheat2000/PhenoGenotype_1.csv"

    # pheno = "DaystoSilk_07CL1"
    # data_file = "maize282/PhenoGenotype_DaystoSilk_07CL1.csv"
    # data_file = "maize282/PhenoGenotype_maize_filter0.5.csv"

    # pheno = "Leaf area"
    # data_file = "sunflower/PhenoGenotype_sunflower_filter0.1.csv"
    # pheno = "Seed rectangular"
    # data_file = "sunflower/PhenoGenotype_Seed rectangular.csv"


    data_file = pd.read_csv(data_file, header=0, index_col=0)
    print("data_file.columns[-30:]: ", data_file.columns[-30:])

    X = data_file.drop([pheno], axis=1)
    X = X.T


    r = allel.rogers_huff_r(X)
    r2 = squareform(r ** 2)
    r2 = pd.DataFrame(r2, columns= X.index, index= X.index)
    # print(r2)
    # print(X.index)

    # r2.to_csv("result/r2" + "_" + "rice" + "_filter0.9.csv", header=True, sep=",",
    #                      index=True)  # index = 0 does not save row index
    # r2.to_csv("result/r2" + "_" + "sunflower" + "_filter0.1.csv", header=True, sep=",",
    #           index=True)  # index = 0 does not save row index
    # r2.to_csv("result/r2" + "_" + "maize" + "_filter0.5.csv", header=True, sep=",",
    #           index=True)  # index = 0 does not save row index
    #
    r2.to_csv("result/r2_maize.csv", header=True, sep=",", index=True)
    return r2
    # fig, ax = plt.subplots()

    # im, cbar = sns.heatmap(data = r2)
    # sns.heatmap(data= r2, xticklabels = X.index)
    # plt.imsave('result/LDHeatmapChr1.png', r2)
    # # fig.tight_layout()
    # plt.show()

def LDHeatmapChromosome(pheno, data_file, Chr = '1'):
    data_file = pd.read_csv(data_file, header=0, index_col=0)

    X = data_file.drop([pheno], axis=1)

    selected_pos = []
    for item in X.columns:
        if item.split('_')[0] == Chr:
            selected_pos.append(item)
    x_selected = X.loc[:, selected_pos]
    x_selected = x_selected.T
    r = allel.rogers_huff_r(x_selected)
    r2 = squareform(r ** 2)
    r2 = pd.DataFrame(r2, columns=x_selected.index, index=x_selected.index)

    r2.to_csv("result/r2_maize" + "_Chr" + Chr + ".csv", header=True, sep=",", index=True)
    return r2

def filterSite(pheno, data_file, r2, threshold = 0.8):
    ColNames = r2.columns.tolist()

    delCol = []  #Deleted columns
    NanCol = []

    for index, row in r2.iterrows():
        #If a null value appears, delete this column
        if pd.isnull(r2[index][ColNames[0]]):
            delCol.append(index)
            NanCol.append(index)
            continue

        if index not in delCol:  #This line is no longer in delCol, so do the following processing
            diagonal = ColNames.index(index)
            for ColName in ColNames[(diagonal+1):]:
                if r2[index][ColName] > threshold and ColName not in delCol:
                    #If r2 is greater than the threshold and this ColName is not in delCol,
                    # add the current ColName to delCol.
                    delCol.append(ColName)

                # print(r2[index][ColName])

        # print("index: ", index)
        # print(row)
        # print(len(row))

    # print("delCol: ", delCol)
    print("len(delCol): ", len(delCol))

    # print("NanCol: ", NanCol)
    print("len(NanCol): ", len(NanCol))

    # data_file = "rice413/PhenoGenotype_Plant heightall.csv"
    # data_file = "maize/PhenoGenotype_EarDia.csv"
    # pheno = "Seed rectangular"
    # data_file = "sunflower/PhenoGenotype_Seed rectangular.csv"
    # data_file = "sunflower/PhenoGenotype_sunflower.csv"

    data_file = pd.read_csv(data_file, header=0, index_col=0)


    data_file = data_file.drop(columns= delCol)

    data_file.to_csv("maize282/PhenoGenotype_" + pheno + "_filter" + str(threshold) + ".csv", header=True, sep=",",
                     index=True)
    # data_file.to_csv("maize282/PhenoGenotype_maize_filter" + str(threshold) + ".csv", header= True, sep= ",", index= True)
    # data_file.to_csv("sunflower/PhenoGenotype_sunflower_filter" + str(threshold) + ".csv", header=True, sep=",", index=True)

    # print(data_file.describe())

def filterSiteChromosome(r2, threshold = 0.8):
    ColNames = r2.columns.tolist()

    delCol = []  #Deleted columns
    NanCol = []

    for index, row in r2.iterrows():
        #If a null value appears, delete this column
        if pd.isnull(r2[index][ColNames[0]]):
            delCol.append(index)
            NanCol.append(index)
            continue

        if index not in delCol:  #This line is no longer in delCol, so do the following processing
            diagonal = ColNames.index(index)
            for ColName in ColNames[(diagonal+1):]:
                if r2[index][ColName] > threshold and ColName not in delCol:
                    # If r2 is greater than the threshold and this ColName is not in delCol,
                    # add the current ColName to delCol.
                    delCol.append(ColName)

                # print(r2[index][ColName])

        # print("index: ", index)
        # print(row)
        # print(len(row))

    # print("delCol: ", delCol)
    print("len(delCol): ", len(delCol))

    # print("NanCol: ", NanCol)
    # print("len(NanCol): ", len(NanCol))
    #
    # # data_file = "rice413/PhenoGenotype_Plant heightall.csv"
    # # data_file = "maize/PhenoGenotype_EarDia.csv"
    # # pheno = "Seed rectangular"
    # # data_file = "sunflower/PhenoGenotype_Seed rectangular.csv"
    # # data_file = "sunflower/PhenoGenotype_sunflower.csv"
    #
    # data_file = pd.read_csv(data_file, header=0, index_col=0)
    #
    #
    # data_file = data_file.drop(columns= delCol)
    #
    # data_file.to_csv("maize282/PhenoGenotype_" + pheno + "_filter" + str(threshold) + ".csv", header=True, sep=",",
    #                  index=True)


def GetMaxR2Position(R2File, length = 13):
    '''
    :param R2File: pairwise squared allele-frequency correlation file
    :param length: Returned site length
    :return: the position list
    '''
    R2File = pd.read_csv(R2File, header=0, index_col=0)

    # 选用代码
    # print("18612:18624")
    # print(R2File.iloc[18612:18624, 18612:18624])
    # print("18624:18636")
    # print(R2File.iloc[18624:18636, 18624:18636])
    # print("18636:18648")
    # print(R2File.iloc[18636:18648, 18636:18648])

    max_by_columns = R2File.max()
    max_column = max_by_columns.idxmax()
    max_index = R2File[max_column].idxmax()

    ColNames = R2File.columns.tolist()
    # print("ColNames[18612:18648]: ", ColNames[18612:18648])

    firstSitePos = ColNames.index(max_column)
    secondSitePos = ColNames.index(max_index)
    middlePos = (firstSitePos + secondSitePos) // 2
    selectedCol = ColNames[middlePos - length // 2 : middlePos + length //2]

    R2File_max = max(max_by_columns)

    print("R2File_max: ", R2File_max)
    print("max_index: ", max_index)
    print("max_column: ", max_column)
    print("R2File.loc[max_index, max_column]): ", R2File.loc[max_index, max_column])
    print("selectedCol: ", selectedCol)

    R2selected = R2File.iloc[middlePos - length // 2 : middlePos + length //2,
                 middlePos - length // 2 : middlePos + length //2] # 指定连续列，用数字
    print("R2selected: ", R2selected)

    return selectedCol


def GetMaxAreaR2Position(R2File, length = 12):
    R2File = pd.read_csv(R2File, header=0, index_col=0)
    ColNames = R2File.columns.tolist()
    selectedCol = []
    Maximum = 0

    for i in range(len(ColNames) - length):
        if R2File.iloc[i:(i+length), i:(i+length)].sum().sum() > Maximum:
            Maximum = R2File.iloc[i:(i+length), i:(i+length)].sum().sum()
            selectedCol = ColNames[i:(i+length)]
            print(R2File.iloc[i:(i+length), i:(i+length)])
            print("Maximum: ", Maximum)
            print()

    print("selectedCol: ", selectedCol)

    return selectedCol


if __name__ == '__main__':
    # DaystoSilk PlantHeight
    # 065 06CL1 06PR 07A 07CL1 26M3
    # pheno = "PlantHeight_26M3"
    # data_file = "maize282/PhenoGenotype_PlantHeight_26M3.csv"
    #
    # r2 = LDHeatmap(pheno, data_file)
    # filterSite(pheno, data_file, r2, threshold = 0.5)

    # R2File = "result/r2_sunflower.csv"
    # GetMaxAreaR2Position(R2File)

    # pheno = "Plant height"
    # data_file = "rice413/PhenoGenotype_" + pheno + "all.csv"
    # pheno = "Leaf area"
    # data_file = "sunflower/PhenoGenotype_sunflower.csv"
    # pheno = "DaystoSilk_07CL1"
    # data_file = "maize282/PhenoGenotype_DaystoSilk_07CL1.csv"

    # t1 = time.time()
    # r2 = LDHeatmap(pheno, data_file)

    # r2 = LDHeatmapChromosome(pheno, data_file, Chr= '1')
    # r2 = LDHeatmapChromosome(pheno, data_file, Chr='2')
    # r2 = LDHeatmapChromosome(pheno, data_file, Chr='3')
    # r2 = LDHeatmapChromosome(pheno, data_file, Chr='4')
    # r2 = LDHeatmapChromosome(pheno, data_file, Chr='5')
    # r2 = LDHeatmapChromosome(pheno, data_file, Chr='6')
    # r2 = LDHeatmapChromosome(pheno, data_file, Chr='7')
    # r2 = LDHeatmapChromosome(pheno, data_file, Chr='8')
    # r2 = LDHeatmapChromosome(pheno, data_file, Chr='9')
    # r2 = LDHeatmapChromosome(pheno, data_file, Chr='10')
    # r2 = LDHeatmapChromosome(pheno, data_file, Chr='11')
    # r2 = LDHeatmapChromosome(pheno, data_file, Chr='12')
    # r2 = LDHeatmapChromosome(pheno, data_file, Chr='13')
    # r2 = LDHeatmapChromosome(pheno, data_file, Chr='14')
    # r2 = LDHeatmapChromosome(pheno, data_file, Chr='15')
    # r2 = LDHeatmapChromosome(pheno, data_file, Chr='16')
    # r2 = LDHeatmapChromosome(pheno, data_file, Chr='17')



    # p1 = multiprocessing.Process(target= LDHeatmapChromosome, args= (pheno, data_file, '1'))
    # p2 = multiprocessing.Process(target=LDHeatmapChromosome, args=(pheno, data_file, '2'))
    # p3 = multiprocessing.Process(target=LDHeatmapChromosome, args=(pheno, data_file, '3'))
    # p4 = multiprocessing.Process(target=LDHeatmapChromosome, args=(pheno, data_file, '4'))
    # p5 = multiprocessing.Process(target=LDHeatmapChromosome, args=(pheno, data_file, '5'))
    # p6 = multiprocessing.Process(target=LDHeatmapChromosome, args=(pheno, data_file, '6'))
    # p7 = multiprocessing.Process(target=LDHeatmapChromosome, args=(pheno, data_file, '7'))
    # p8 = multiprocessing.Process(target=LDHeatmapChromosome, args=(pheno, data_file, '8'))
    # p9 = multiprocessing.Process(target=LDHeatmapChromosome, args=(pheno, data_file, '9'))
    # p10 = multiprocessing.Process(target=LDHeatmapChromosome, args=(pheno, data_file, '10'))
    # p11 = multiprocessing.Process(target=LDHeatmapChromosome, args=(pheno, data_file, '11'))
    # p12 = multiprocessing.Process(target=LDHeatmapChromosome, args=(pheno, data_file, '12'))
    # p13 = multiprocessing.Process(target=LDHeatmapChromosome, args=(pheno, data_file, '13'))
    # p14 = multiprocessing.Process(target=LDHeatmapChromosome, args=(pheno, data_file, '14'))
    # p15 = multiprocessing.Process(target=LDHeatmapChromosome, args=(pheno, data_file, '15'))
    # p16 = multiprocessing.Process(target=LDHeatmapChromosome, args=(pheno, data_file, '16'))
    # p17 = multiprocessing.Process(target=LDHeatmapChromosome, args=(pheno, data_file, '17'))
    #
    # p1.start()
    # p2.start()
    # p3.start()
    # p4.start()
    # p5.start()
    # p6.start()
    # p7.start()
    # p8.start()
    # p9.start()
    # p10.start()
    # p11.start()
    # p12.start()
    # p13.start()
    # p14.start()
    # p15.start()
    # p16.start()
    # p17.start()
    #
    # p1.join()
    # p2.join()
    # p3.join()
    # p4.join()
    # p5.join()
    # p6.join()
    # p7.join()
    # p8.join()
    # p9.join()
    # p10.join()
    # p11.join()
    # p12.join()
    # p13.join()
    # p14.join()
    # p15.join()
    # p16.join()
    # p17.join()


    # t2 = time.time()
    # print("t2 - t1: ", t2 - t1)

    for i in range(10):
        print("i+1:", i+1)
        r2 = "result/r2_maize_Chr" + str(i+1) + ".csv"
        r2 = pd.read_csv(r2, header=0, index_col=0)
        filterSiteChromosome(r2, threshold= 0.5)