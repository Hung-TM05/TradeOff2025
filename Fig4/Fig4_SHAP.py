import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import OneHotEncoder
import shap
from matplotlib import pyplot as plt

# date = 'test_250110'
# set initial seed (=0) which is the anchor of seed after random forest loop
org_seed = my_seed = 0
# set main seed
main_seed = 99
np.random.seed(main_seed)
location_list = ['All']
coor_list = ['CTmin', 'TTrange'] #'CTmax', 
df_all = pd.read_csv('/Users/tmhung/Library/CloudStorage/OneDrive-個人/ESSLAB/1. Manuscript/1. Trade-off/0. Nature EE format/Final codes/Empirical_moth_data.csv')

# In[] 3 Locs predictors
for loc in location_list:
    my_seed = org_seed
    for coor in coor_list:
        if loc == 'All':
            df = df_all
        else:
            df = df_all[df_all.Location == loc]
                
        y_columns = coor
                        
        df.dropna(axis=0, how='any', inplace=True)
        
        y_ = df[y_columns]
        if coor == 'CTmax':
            x_ = df[['Family', 'B_length','CTmin',  
                     'STmax', 'STmin', 'STmean' ]]
        elif coor == 'CTmin':
            x_ = df[['Family', 'B_length','CTmax',  
                     'STmax', 'STmin', 'STmean' ]]
        elif coor == 'TTrange':
            x_ = df[['Family', 'B_length',
                     'STmax', 'STmin', 'STR' ]]
        
        target = coor
        # for target in targets:
        y = y_.values
        rfr_df = None
        for i in range(10):
            print(f'Processing... {i+1}', end='\r')
            
            rfr = RandomForestRegressor(n_estimators=500, max_samples=.7, oob_score=True, random_state=my_seed)
            ohe = OneHotEncoder(sparse=False)
            fam_ohe = ohe.fit_transform(df.Family.values.reshape(-1,1))
            
            x_cols = list(ohe.categories_[0]) + list(x_.columns[1:])
            
            x = np.concatenate([fam_ohe, x_.iloc[:,1:].values], axis=1)
            
            
            rfr.fit(x, y)
            feature_importances_ = rfr.feature_importances_ * 100
            rfr_dict = dict(zip(x_cols, feature_importances_))
            rfr_dict['r2_oob'] = rfr.oob_score_
            rfr_dict['r2_x'] = rfr.score(x, y)
            
            explainer = shap.TreeExplainer(rfr)
            # explainer = shap.Explainer(rfr)
            
            if rfr_df is None:
                rfr_df = pd.DataFrame(rfr_dict, index=[i])
            else:
                rfr_df = pd.concat([rfr_df, pd.DataFrame(rfr_dict, index=[i])])
            if i < 9:
                shap_values = explainer.shap_values(x)
                df_ = pd.DataFrame(shap_values)
                df_.to_csv(f'/Users/tmhung/Library/CloudStorage/OneDrive-個人/ESSLAB/1. Manuscript/1. Trade-off/0. Nature EE format/Final codes/SHAP results/SHAP_value_csv/10_rfr_{target}_SHAP_value_{i}_{loc}_seed_{my_seed}.csv', index=False, sep=',')
            else:
                shap_values2 = explainer(x)
                shap_values22 = shap_values2.values
                df_ = pd.DataFrame(shap_values22)
                df_.to_csv(f'/Users/tmhung/Library/CloudStorage/OneDrive-個人/ESSLAB/1. Manuscript/1. Trade-off/0. Nature EE format/Final codes/SHAP results/SHAP_value_csv/10_rfr_{target}_SHAP_value_{i}_{loc}_seed_{my_seed}.csv', index=False, sep=',')
           
            print(f'Done. Seed: {my_seed}')
            feature_importances = rfr_df.iloc[:,:-2] 
            # rfr_df.to_csv(f'rfr_{target}.csv', index=False, sep=',')
            pd.concat([rfr_df[feature_importances.mean().sort_values().index.values], rfr_df.iloc[:,-2:]], axis=1).to_csv(f'/Users/tmhung/Library/CloudStorage/OneDrive-個人/ESSLAB/1. Manuscript/1. Trade-off/0. Nature EE format/Final codes/SHAP results/10_rfr_All_{target}_{loc}.csv', index=False, sep=',')           
            my_seed = my_seed+1
        
          
        colname = list(ohe.categories_[0])
        # for n in range(1,21):
        #     fam_ = x_cols[n]
        #     colname.append(fam_)
        if coor == 'CTmax':
            ls_col = ['Body length', 'CTmin', 'STmax','STmin','STmean']
        elif coor == 'CTmin':
            ls_col = ['Body length', 'CTmax', 'STmax','STmin','STmean']
        elif coor == 'TTrange':
            ls_col = ['Body length', 'STmax','STmin','STR']
        for m in ls_col:
            colname.append(m)
    
        li = []
        seed = my_seed-10
        for k in range(9):
            df = pd.read_csv(f'/Users/tmhung/Library/CloudStorage/OneDrive-個人/ESSLAB/1. Manuscript/1. Trade-off/0. Nature EE format/Final codes/SHAP results/SHAP_value_csv/10_rfr_{target}_SHAP_value_{k}_{loc}_seed_{seed}.csv')
            df.columns = colname
            li.append(df)
            seed = seed+1
    
        df_10 = pd.DataFrame(shap_values2.values)
        df_10.columns = colname
        li.append(df_10)
    
        from functools import reduce
        df_sum = reduce(lambda x, y: pd.DataFrame.add(x, y, fill_value=0), li)
        df_sum = df_sum/10
        colname_sum = []
        for n in range(0,len(colname)):
            colname_sum.append(n)
        df_sum.columns = colname_sum
            
        shap_values2.values = df_sum
        shap_values2.feature_names = colname
    
        ## workable plotting method
        # shap.plots.bar(shap_values2)
        shap.plots.beeswarm(shap_values2, show=False)
        plt.tight_layout()
        plt.savefig(f'/Users/tmhung/Library/CloudStorage/OneDrive-個人/ESSLAB/1. Manuscript/1. Trade-off/0. Nature EE format/Final codes/SHAP results/SHAP_plots/{loc}_{coor}_beeswarm.pdf',format="pdf", bbox_inches="tight")
        plt.close()
        clust = shap.utils.hclust(x, y, linkage="single")
        shap.plots.bar(shap_values2, show=False)
        plt.tight_layout()
        plt.savefig(f'/Users/tmhung/Library/CloudStorage/OneDrive-個人/ESSLAB/1. Manuscript/1. Trade-off/0. Nature EE format/Final codes/SHAP results/SHAP_plots/{loc}_{coor}_bar.pdf',format="pdf", bbox_inches="tight")
        plt.close()
        