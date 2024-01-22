import sklearn

#import dataframes and merge
import pandas as pd
df1 = pd.read_csv('SRR3401407.csv')
df2 = pd.read_csv('SRR3401415.csv')
df3 = pd.read_csv('SRR3401416.csv')
df4 = pd.read_csv('SRR3401417.csv')
df5 = pd.read_csv('SRR3401418.csv')   
df = pd.concat([df1, df2, df3, df4, df5], axis=0, ignore_index=True)


# Keep SNPs only - regular expression pattern to match SNP variants
snp_pattern = r'^chr\d+:\d+:[ACGT]>[ACGT]$'
# boolean mask to filter SNPs
snp_mask = df['Variant'].str.match(snp_pattern)
# Apply  mask to the dataframe to keep only SNP rows
snp_df = df[snp_mask]
# Reset the index to make the resulting dataframe clean
snp_df.reset_index(drop=True, inplace=True)

import pandas as pd
from sklearn.model_selection import train_test_split

#drop rows with missing values and 'Read Depth' 1 or lower
snp_df = snp_df.dropna()
snp_df = snp_df[snp_df['Read Depth'] > 1]

#scale **after** dropping READ DEPTH 1 or lower
from sklearn.preprocessing import StandardScaler
column_to_scale = 'Read Depth'
scaler = StandardScaler()
# Reshape Read Depth column & fit-transform the data
snp_df.loc[:, [column_to_scale]] = scaler.fit_transform(snp_df[[column_to_scale]].values.reshape(-1, 1))

# Split the data into features and targets
X = snp_df.drop(['High Confidence', 'Variant','Reads Supporting Reference','Reads Supporting Alternate',"Sample"], axis=1)  # features
y = snp_df['High Confidence']  # targets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=0, shuffle=True)


# imbalance handling and training model
from imblearn.under_sampling import RandomUnderSampler
from sklearn.ensemble import RandomForestClassifier
rcf_opt = {'bootstrap': True, 'criterion': 'gini', 'max_depth': 233, 'max_features': 'log2', 'min_samples_leaf': 2, 'min_samples_split': 2, 'n_estimators': 817}
rfc_model = RandomForestClassifier(random_state=0, **rcf_opt)
rus = RandomUnderSampler(random_state=0, sampling_strategy=0.4) 
X_rus, y_rus = rus.fit_resample(X_train, y_train)
rfc_model.fit(X_rus, y_rus)


class_probabilities = rfc_model.predict_proba(X_test)
predicted_labels = (class_probabilities[:, 1] >= 0.8).astype(int)
X_test['Class1_Probabilities'] = class_probabilities[:, 1]
X_test['Actual_Labels'] = y_test
X_test['Predicted_Labels'] = predicted_labels

from sklearn.metrics import average_precision_score
y_scores = X_test['Class1_Probabilities']
y_true = X_test['Actual_Labels']
average_precision = average_precision_score(y_true, y_scores)
print('Average precision-recall score: {0:0.2f}'.format(average_precision))