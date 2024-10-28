import ase.db
import numpy as np
from tpot import TPOTRegressor
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import math
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
import pandas as pd
import matplotlib.pyplot as plt
import sklearn.gaussian_process as gp
from sklearn.pipeline import Pipeline
import pandas as pd
import pickle

db = ase.db.connect('cdpbs.db')


# Load data, select data by 3 sigma
feature = np.loadtxt('feature.csv', delimiter=',',encoding="UTF-8-sig")
target = np.loadtxt('energy.csv', delimiter=',',encoding="UTF-8-sig")
ratio = np.loadtxt('ratio.csv', delimiter=',',encoding="UTF-8-sig")
pre_ratio = np.loadtxt('pre_ratio.csv', delimiter=',',encoding="UTF-8-sig")

X_train=feature[0:92,:]
y_train=target[0:92]
print(X_train.shape)





tpot = TPOTRegressor(generations = 5, population_size = 32, offspring_size = 32, verbosity = 2, scoring = 'neg_median_absolute_error', n_jobs = -1, random_state = 42)

pre_feature = np.loadtxt('pre_feature.csv', delimiter=',',encoding="UTF-8-sig")
X_pre=pre_feature[0:1001,:]




tpot.fit(X_train, y_train)
y_predict = tpot.predict(X_pre)

np.savetxt('pre_energy.csv',y_predict,delimiter=',')

best_tpot_model = tpot.fitted_pipeline_
fw = open('best_model','wb')
pickle.dump(best_tpot_model,fw)
fw.close()


# target_predict, sigma_predict = gpr_model.predict(pre_feature, return_std=True)
# Id=[]
# for row in db.select():
#     Id.append(row.id)
#
# Id = np.array(Id)
# print(type(Id),type(pre_ratio),type(sigma_predict))
#
# print(len(Id),len(pre_ratio),len(sigma_predict))
#
# Feature = np.vstack((Id,pre_ratio,sigma_predict)).T
#
#
# np.savetxt('select.csv',Feature,delimiter=',')


