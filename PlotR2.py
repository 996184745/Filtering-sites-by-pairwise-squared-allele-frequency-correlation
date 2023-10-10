import pandas as pd
import numpy as np


def AggregateR2(dataFile, RatioOfAggregation):
    '''
    :param datafile: source file of r2
    :param RatioOfAggregation: like 4 in 1, 5 in 1, is an integer
    :return: a new matrix file of r2, but much smaller than datafile
    '''
    dataFile = pd.read_csv(dataFile, header=0, index_col=0)
    dataFile.fillna(value= 0, inplace= True)
    # dataFile.dropna(axis= 1, how= 'any', inplace=True)
    # dataFile.dropna(axis= 0, how= 'any', inplace= True)

    r2Matrix = np.array(dataFile)

    print("r2Matrix.shape: ", r2Matrix.shape)

    AggregationMatrixLength = len(r2Matrix) // RatioOfAggregation
    r2Aggregation = np.zeros((AggregationMatrixLength, AggregationMatrixLength))
    for i in range(AggregationMatrixLength):
        for j in range(AggregationMatrixLength):
            r2Aggregation[i][j] = r2Matrix[i*RatioOfAggregation : (i+1)* RatioOfAggregation, j*RatioOfAggregation: (j+1)*RatioOfAggregation].sum()

    np.savetxt("result/r2_sunflower_12_Aggregation.csv", r2Aggregation, delimiter= ",")

    m = np.argmax(r2Aggregation)
    r, c = divmod(m, r2Aggregation.shape[1])
    print("np.max(r2Aggregation): ", np.max(r2Aggregation))
    print("r, c: ", r, c)

    return r2Aggregation




# dataFile = "result/r2_rice_filter0.9.csv"
# RatioOfAggregation = 1000
dataFile = "result/r2_sunflower.csv"
RatioOfAggregation = 12
AggregateR2(dataFile, RatioOfAggregation)