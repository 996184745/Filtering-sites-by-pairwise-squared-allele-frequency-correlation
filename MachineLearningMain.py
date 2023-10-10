import time
import numpy as np
import pandas as pd
import scipy.stats as stats
from sklearn.linear_model import Ridge
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold, StratifiedShuffleSplit
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt

import allel
from scipy.spatial.distance import squareform

def main():
    # When r2 = 1, it means that the linkage is completely unbalanced and there is no recombination.
    #
    # When r2 = 0, it means that the linkage is completely balanced and the random group




    model = Ridge(random_state=0)
    # model = RandomForestRegressor(random_state=0)


    # pheno = "Plant height"
    # pheno = "Seed number per panicle"
    # pheno = "Seed surface area"
    # data_file = "rice413/PhenoGenotype_Plant height_filter0.9.csv"
    # data_file = "rice413/PhenoGenotype_Seed surface areaall.csv"
    # data_file = "rice413/PhenoGenotype_Seed surface area_filter0.9.csv"


    # pheno = "Leaf area"
    # data_file = "sunflower/PhenoGenotype_sunflower.csv"
    # data_file = "sunflower/PhenoGenotype_sunflower_filter0.9.csv"
    # pheno = "Flower head diameter"
    # data_file = "sunflower/PhenoGenotype_Flower head diameter.csv"
    # data_file = "sunflower/PhenoGenotype_Flower head diameter_filter0.9.csv"
    # pheno = "Seed rectangular"
    # data_file = "sunflower/PhenoGenotype_Seed rectangular.csv"
    # data_file = "sunflower/PhenoGenotype_Seed rectangular_filter0.9.csv"


    # DaystoSilk PlantHeight
    # 065 06CL1 06PR 07A 07CL1 26M3
    pheno = "PlantHeight_26M3"
    data_file = "maize282/PhenoGenotype_PlantHeight_26M3_filter0.5.csv"
    # data_file = "maize282/PhenoGenotype_maize_filter0.9.csv"





    data_file = pd.read_csv(data_file, header= 0, index_col= 0)

    y = data_file[pheno]
    X = data_file.drop([pheno], axis=1)

    X = np.array(X)

    pearns = []
    kf = KFold(n_splits=10, shuffle= True, random_state= 0)
    t1 = time.time()
    for train_index, test_index in kf.split(X):
        # print("TRAIN:", train_index, "TEST:", test_index)
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]

        model.fit(X_train, y_train)


        y_pre = model.predict(X_test)

        pearn, p = stats.pearsonr(y_pre, y_test)
        pearns.append(pearn)
    t2 = time.time()
    print("model training time: ", t2 - t1)

    t3 = time.time()
    model.predict(X_test[[0], :])
    t4 = time.time()
    model.predict(X_test[[1], :])
    t4 = time.time()
    model.predict(X_test[[2], :])
    t4 = time.time()
    model.predict(X_test[[3], :])
    t4 = time.time()
    model.predict(X_test[[4], :])
    t4 = time.time()
    model.predict(X_test[[5], :])
    t4 = time.time()
    model.predict(X_test[[6], :])
    t4 = time.time()
    model.predict(X_test[[7], :])
    t4 = time.time()
    model.predict(X_test[[8], :])
    t4 = time.time()
    model.predict(X_test[[9], :])
    t4 = time.time()
    print("model inference time: ", (t4 - t3)/10)


    print(pearns)
    print(np.mean(pearns))







    # FI = pd.Series(model.feature_importances_, index = X.columns.values)  # pySpark
    # positionList = list(FI.index)
    # FI = FI.sort_values(ascending=False)
    # featureImportances = pd.DataFrame({'POSITION': FI.index, "featureImportances": FI.values})
    #
    # #================linkage disequilibrium==========================
    # r2List = [None] * len(positionList)
    # LDArray = pd.DataFrame()
    # wide = 100  #number of neighbors
    # for i in range(50):
    #     centerPoint = featureImportances.at[i, 'POSITION']
    #     centerIndex = positionList.index(centerPoint)
    #
    #
    #     for j in range(centerIndex - wide, centerIndex + wide + 1, 1):
    #         LDArray[positionList[j]] = X[positionList[j]]
    #     LDArray = LDArray.T
    #
    #     r = allel.rogers_huff_r(LDArray)
    #     r2 = squareform(r ** 2)
    #
    #     left = sum(r2[:wide, wide]) / wide
    #
    #     right = sum(r2[(wide+1):, wide]) / wide
    #
    #     if np.isnan(left):
    #         r2List[i] = right
    #     elif np.isnan(right):
    #         r2List[i] = left
    #     else:
    #         r2List[i] = max(left, right)
    #     print("i:", i)
    #     print("left:", len(r2[:wide, wide]))
    #     print("right:", len(r2[(wide + 1):, wide]))
    #     print("left:", left)
    #     print("right:", right)
    #     print("r2List[i]:", r2List[i])
    #     #
    #     # r2List[i] = sum(r2[:,wide]) / wide / 2
    #     # print("r2List[i]:", r2List[i])
    #
    #     LDArray = pd.DataFrame()  #clear again
    #
    # featureImportances["LD"] = r2List
    # # ================linkage disequilibrium==========================
    #
    #
    # featureImportances.to_csv("result/"+ pheno + type(model).__name__ + ".csv", header= True)

    #=====================shap=======================================
    # explainer = shap.TreeExplainer(model)
    # shap_values = explainer.shap_values(X)
    # print(shap_values)
    # =====================shap=======================================


    #==============permutation_importance============================
    # perm_importance = permutation_importance(model, X, y)
    # sorted_importances_idx = perm_importance.importances_mean.argsort()
    # importances = pd.DataFrame(
    #     perm_importance.importances[sorted_importances_idx].T,
    #     columns=X.columns[sorted_importances_idx],
    # )
    # ax = importances.plot.box(vert=False, whis=10)
    # ax.set_title("Permutation Importances (test set)")
    # ax.axvline(x=0, color="k", linestyle="--")
    # ax.set_xlabel("Decrease in accuracy score")
    # ax.figure.tight_layout()
    # plt.show()
    # ==============permutation_importance============================


    # fig = plt.figure(figsize=(12, 5))
    # plt.bar(FI.index, FI.values, color="blue")
    # plt.xlabel('features')
    # plt.ylabel('importances')
    # plt.show()


if __name__ == '__main__':
    main()