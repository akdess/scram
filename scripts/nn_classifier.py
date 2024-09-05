# ******
# readme
#
# [03/06/2023] adata.h5ad is accepted as reference or test dataset
#
#              **** Important assumptions: ****
#
#              -----------------------------
#              Reference and test file paths
#              -----------------------------
#              4. ref_directory directly contains adata.h5ad OR { genes.tsv, barcodes.tsv, matrix.mtx }, same for test_directory
#
#              ------------------
#              Reference metadata
#              ------------------                 
#              5. If adata.h5ad is loaded as reference, metadata has been included thus no additional metadata file is needed
#
# [03/06/2023] scanpy preprocessing workflow added
#              df_ref replaced by adata_ref
#              df_test replaced by adata_test
#
# [03/06/2023] Updated to be compatible with both sparse and dense input
#              Optimized to speed up for sparse input
#
#
#
# [02/16/2023] Updated with user input arguments and automatic input check
#
#              **** Important assumptions: ****
#
#              ----------
#              User input
#              ----------
#              1. Arguments are provided by user in order (see examples below): 
#                 $ python nn_classifier.py [ref_name] [ref_directory] [ref_class_col] [test_name] [test_directory] --marker [filename] --shap
#              
#              2. Positional arguments [ref_name] [ref_directory] [ref_class_col] [test_name] [test_directory] are all provided by user
#              3. Each positional or optional argument string doesn't contain space 
#
#              -----------------------------
#              Reference and test file paths
#              -----------------------------
#              4. ref_directory directly contains genes.tsv, barcodes.tsv, and matrix.mtx, same for test_directory
#
#              ------------------
#              Reference metadata
#              ------------------
#              5. ref_directory directly contains metadata file, filename is meta.tsv, meta.csv, or meta.xlsx
#              6. If filename is meta.tsv, sep='\t', else if filename is meta.csv, sep=',' 
#              7. Cell order in metadata is the same as reference dataset
#
#              ---------------------------------
#              Reference marker genes (Optional)
#              ---------------------------------
#              8. If the file is provided by user, ref_directory directly contains it, filename ends with .tsv, .csv, or .xlsx
#              9. If filename ends with .tsv, sep='\t', else if filename ends with .csv, sep=',' 
#              10. In this file, the gene symbol column name is Gene 
#
#
# [02/16/2023] Commented out previous data loading blocks
# 
#
#
#
# This python file loads one reference dataset (df_ref) and one test dataset (df_test).
# It selects overlapped genes and trains a neural network model as a cell type classifier based on the reference dataset.
# Then it predicts the cell type probabilities and extracts features for the test dataset.
# Finally, it evaluates the model using samples from test dataset by SHAP.
# 
# To save memory usage, intermediate data are deleted explicitly.
# 
# 
# The code is developed in Python 3.7.3. Install important packages with specific versions in command line:
# $ pip install tensorflow==1.13.1
# $ pip install keras==2.2.4
# 
# Install other packages in command line:
# $ pip install pandas
# $ pip install scipy
# $ pip install scikit-learn
# $ pip install matplotlib
# $ pip install joblib
# $ pip install shap
# $ pip install openpyxl
# $ pip install scanpy
# 
# 
# [02/16/2023] 
# Run this file in command line with arguments: 
# $ python nn_classifier.py [ref_name] [ref_directory] [ref_class_col] [test_name] [test_directory] --marker [filename] --shap
#
# Example 1: Use cerebellum as reference and meningioma as test, with reference marker genes provided, opt to run SHAP
# $ python nn_classifier.py cerebellum ./developing_human_cerebellum/ Cluster meningioma ./meningioma/ --marker CellTypeMarker_DevelopingHumanData.xlsx --shap
#
# Simpler command for Example 1:
# $ python nn_classifier.py cerebellum developing_human_cerebellum Cluster meningioma meningioma -m CellTypeMarker_DevelopingHumanData.xlsx -s
# 
# Example 2: Use codex as reference and glioma as test, with reference marker genes provided, not to run SHAP
# $ python nn_classifier.py codex codex Cluster glioma glioma -m codex_cluster_markers.xlsx
#
# Example 3: Use codex as reference and DIPG as test, reference marker genes not provided, opt to run SHAP
# $ python nn_classifier.py codex codex Cluster DIPG DIPG -s
#
# See details in help messages in code below.
#
# ******


# Import packages
print('Importing packages...')

import os
import random
import numpy as np
import pandas as pd
import tensorflow as tf
from scipy.io import mmread
from sklearn.preprocessing import MaxAbsScaler, LabelEncoder
from tensorflow.keras.layers import Dense, Dropout, BatchNormalization, Activation
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.utils import to_categorical
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras import backend as K
from tensorflow.keras import regularizers
from joblib import dump
import shap
import matplotlib.pyplot as plt
import argparse
import scanpy as sc
from sklearn.preprocessing import MinMaxScaler, LabelEncoder
from sklearn.utils import class_weight
from scipy.sparse import isspmatrix
from tensorflow.keras.optimizers import SGD


def obs_key_wise_subsampling(adata, obs_key, N):
    adatas = [adata[adata.obs[obs_key].isin([clust])] for clust in adata.obs[obs_key].astype('category').cat.categories]
    for dat in adatas:
        if dat.n_obs > N:
            sc.pp.subsample(dat, n_obs=N, random_state=0)
    return adatas[0].concatenate(*adatas[1:]).copy()

print('Packages imported.')
print()


# Set random seeds for reproducibility
os.environ['PYTHONHASHSEED'] = '1000'
random.seed(1000)
np.random.seed(1000)
tf.random.set_seed(1000)
#session_conf = tf.ConfigProto(intra_op_parallelism_threads=1, inter_op_parallelism_threads=1)
#sess = tf.Session(graph=tf.get_default_graph(), config=session_conf)
#K.set_session(sess)


# **** User command line input ****
# **** Example inputs are provided in help messages, with cerebellum used as reference and meningioma as test ****
parser = argparse.ArgumentParser(description='neural network-based cell type classifier')

# Positional inputs for reference dataset
parser.add_argument('ref_name', help='Input name of reference dataset for model training. e.g., cerebellum')
parser.add_argument('ref_directory', help='Input reference directory that contains genes.tsv, barcodes.tsv, matrix.mtx, and metadata. e.g., ./developing_human_cerebellum/')
parser.add_argument('ref_class_col', help='Input the column name in metadata that contains cell type information. e.g., Cluster')

# Positional inputs for test dataset
parser.add_argument('test_name', help='Input name of test dataset for cell type classification. e.g., meningioma')
parser.add_argument('test_directory', help='Input test directory that contains genes.tsv, barcodes.tsv, and matrix.mtx. e.g., meningioma')

# Optional input for reference marker genes, assume the file is in reference directory, filename ends with .tsv, .csv, or .xlsx
parser.add_argument('-m', '--marker', help='Input -m [filename] or --marker [filename] that is in reference directory and contains marker gene information. e.g., -m CellTypeMarker_DevelopingHumanData.xlsx')

# Optional on/off flag for SHAP evaluation
parser.add_argument('-s', '--shap', action='store_true', help='Add -s or --shap at the end of command line if you opt to run SHAP, otherwise it will be skipped. e.g., -s')
# Optional on/off flag to normalize reference dataset
parser.add_argument('-r', '--normalize_ref', action='store_true', help='Add -r or --normalize_ref if you opt to normalize reference dataset, otherwise it will not be normalized. e.g., -r')

# Optional on/off flag to normalize testing dataset
parser.add_argument('-t', '--normalize_test', action='store_true', help='Add -t or --normalize_test if you opt to normalize testing dataset, otherwise it will not be normalized. e.g., -t')

parser.add_argument('-v', '--var_genes', action='store_true', help='Add -t or --var_genes e.g., -v')


# **** User input finished ****


# **** Parse input arguments ****
args = parser.parse_args()

# Arguments for reference dataset
ref_name = args.ref_name
ref_class_col = args.ref_class_col
ref_directory = args.ref_directory

# Arguments for test dataset
test_name = args.test_name
test_directory = args.test_directory 

# Argument for reference marker genes 
marker = args.marker

# Argument for using SHAP
run_shap = args.shap 
# Argument to normalize reference dataset
normalize_ref = args.normalize_ref

# Argument to normalize testing dataset
normalize_test = args.normalize_test

var_genes = args.var_genes

# **** Input arguments parsed ****


# **** Load all datasets according to user input ****
print('**** Loading all datasets according to user input... ****')
print()

# ----------------------
# Load reference dataset
# ----------------------
print(f'Loading {ref_name} as reference dataset...')

# Assume reference directory directly contains adata.h5ad OR { genes.tsv, barcodes.tsv, matrix.mtx }
#N=2000
adata_ref = sc.read_h5ad(ref_directory)
#adata_ref = obs_key_wise_subsampling(adata_ref,  'CellType', N).copy()
N=2000
adata_ref = obs_key_wise_subsampling(adata_ref,  ref_class_col, N)
df_ref_meta = adata_ref.obs.copy()


print(f'{ref_name} loaded as reference dataset.')
print('--------')

# In metadata cell type strings, replace each '/' with '.' as they will appear in output filenames
assert ref_class_col in df_ref_meta.columns, \
    f'IndexNotFoundError: {ref_class_col} column not found in {ref_name} metadata'
df_ref_meta[ref_class_col] = df_ref_meta[ref_class_col].map(lambda x: x.replace('/', '.'))


# -----------------
# Load test dataset
# -----------------
print(f'Loading {test_name} as test dataset...')

# Assume test directory directly contains adata.h5ad OR { genes.tsv, barcodes.tsv, matrix.mtx }
adata_test = sc.read_h5ad(test_directory)
adata_test_raw = adata_test.copy()


print(f'{test_name} loaded as test dataset.')
print('--------')


# --------------------------------------
# Load reference marker genes (Optional)
# --------------------------------------
if marker is None:
    # Not provided by user
    df_ref_markers = None
    print('Reference marker genes not provided by user, thus skipped.')
    print()

else:
    # Assume filename ends with .tsv, .csv, or .xlsx
    assert marker.endswith('.tsv') or marker.endswith('.csv') or marker.endswith('.xlsx'), \
        f'FileTypeError: {marker} provided by user is not a .tsv, .csv, or .xlsx file.'
    
    print(f'Loading {marker} for retrieving reference marker genes...')
    
    # Assume the file is in reference directory
    assert os.path.isfile(marker), \
        f'FileNotFoundError: {marker} not found in {marker}'
    
    # If filename ends with .tsv, assume sep='\t', else if filename ends with .csv, assume sep=',' 
    if marker.endswith('.tsv'):
        df_ref_markers = pd.read_csv(marker, sep='\t')
    elif marker.endswith('.csv'):
        df_ref_markers = pd.read_csv( marker, sep=',')
    else:
        df_ref_markers = pd.read_excel( marker)

    # In df_ref_markers, assume column name of gene symbols is Gene 
    gene_col = 'Gene'
    assert gene_col in df_ref_markers.columns, \
        f'IndexNotFoundError: {gene_col} column not found in {marker}'

    df_ref_markers = df_ref_markers.set_index(gene_col)

    print(f'{marker} loaded and {gene_col} column retrieved as reference marker genes.')
    print()

# **** Datasets all loaded according to user input ****
print('**** Datasets all loaded according to user input. ****')
print()


# **** Neural network model training and prediction starts ****
# Create a directory for output
output_dir = f'{test_name}_{ref_name}_neuralNetwork/'
os.makedirs(output_dir, exist_ok=True)


# scanpy preprocessing 
# Scale each gene into interval (min=-1, max=1) for neuralNetwork 
# Use MaxAbsScaler instead of MinMaxScaler for compatibility with sparse input
#sc.pp.filter_genes(adata_ref, min_counts=1) 
#target_sum=1e6
if normalize_ref:
    sc.pp.normalize_total(adata_ref,target_sum=1e4)
    sc.pp.log1p(adata_ref)
sc.pp.scale(adata_ref, zero_center=True, max_value=10)
if var_genes:
    sc.pp.highly_variable_genes(adata_ref, n_top_genes=2000)

if normalize_test:
    sc.pp.normalize_total(adata_test,target_sum=1e4)
    sc.pp.log1p(adata_test)
sc.pp.scale(adata_test, zero_center=True, max_value=10)

# Select common genes that overlap among training data, testing data, and (optional) training marker genes
if df_ref_markers is None:
    if var_genes:
        common_genes = set(adata_ref.var_names) & set(adata_test.var_names) & set(adata_ref.var_names[adata_ref.var.highly_variable])
    common_genes = set(adata_ref.var_names) & set(adata_test.var_names)
else:
    common_genes = set(adata_ref.var_names) & set(adata_test.var_names) & set(df_ref_markers.index)
    if var_genes:
        common_genes = common_genes | set(adata_ref.var_names[adata_ref.var.highly_variable])    
    common_genes= set(adata_test.var_names) & common_genes
    
common_genes = sorted(list(common_genes))

#print("Highly variable genes: %d"%sum(adata_ref.var.highly_variable))
#print("Used genes: %d"%sum(common_genes))
with open(output_dir + 'common_genes.tsv', 'w') as f:
    f.write('\n'.join(common_genes) + '\n')
    print('common_genes saved.')

adata_ref = adata_ref[:, common_genes]
adata_test = adata_test[:, common_genes]

# Split adata_ref into 80% training and 20% validation
print('Split adata_ref into 80% training and 20% validation...')
cells_train, cells_valid = [], []
for celltype, df in adata_ref.obs.groupby(ref_class_col):
    if df.shape[0] < 10:
        print(f'**** Note: {celltype} is excluded for too few samples. ****')
        continue
    i = int(len(df.index) * 0.8)
    cells_train += df.index[:i].tolist()
    cells_valid += df.index[i:].tolist()

adata_train = adata_ref[cells_train, :]
adata_valid = adata_ref[cells_valid, :]
del adata_ref

# One-hot encoding of celltypes in training data
label_encoder = LabelEncoder()
y_train = label_encoder.fit_transform(adata_train.obs[ref_class_col].values)
dump(label_encoder, output_dir + 'label_encoder.joblib')

num_classes = len(np.unique(y_train))
y_train = to_categorical(y_train, num_classes=num_classes)

# One-hot encoding of celltypes in validation data
y_valid = label_encoder.transform(adata_valid.obs[ref_class_col].values)
y_valid = to_categorical(y_valid, num_classes=num_classes)

# Scale each gene in training data
scaler = MinMaxScaler()
if isspmatrix(adata_train.X):
    adata_train.X = adata_train.X.toarray()
X_train = scaler.fit_transform(adata_train.X)
#X_train =adata_train.X
dump(scaler, output_dir + 'scaler.joblib')
del adata_train

# Scale each gene in validation data
if isspmatrix(adata_valid.X):
    adata_valid.X = adata_valid.X.toarray() 
X_valid = scaler.transform(adata_valid.X)
#X_valid =adata_valid.X
del adata_valid

# Scale each gene in testing data
# scaler = load(output_dir + 'scaler.joblib')
scaler2 = MinMaxScaler()
if isspmatrix(adata_test.X):
    adata_test.X = adata_test.X.toarray() 
X_test = scaler.transform(adata_test.X)
#X_test=adata_test.X
cells_test = adata_test.obs_names
del adata_test

"""model = Sequential()
model.add(Dense(256, activation='relu', input_dim=X_train.shape[1]))
model.add(Dropout(0.2))
model.add(Dense(64, activation='relu'))
model.add(Dropout(0.2))
model.add(Dense(32, activation='relu', name='feature_extraction'))
model.add(Dropout(0.2))
model.add(Dense(num_classes, activation='softmax'))"""
#### without drouput it gets very bad
model = Sequential()
model.add(Dense(256, input_dim=X_train.shape[1]))
model.add(BatchNormalization())
model.add(Activation('relu'))
model.add(Dropout(0.1))
model.add(Dense(64))
model.add(BatchNormalization())
model.add(Activation('relu'))
model.add(Dropout(0.1))
model.add(Dense(32,  name='feature_extraction'))
model.add(BatchNormalization())
model.add(Activation('relu'))
model.add(Dropout(0.1))
model.add(Dense(num_classes, activation='softmax'))

final_learning_rate = 0.001 
initial_learning_rate = 0.1
epochs=100
learning_rate_decay_factor = (final_learning_rate / initial_learning_rate)**(1/epochs)
steps_per_epoch = int(X_train.shape[1]/32)
#steps_per_epoch = 1500

lr_schedule = tf.keras.optimizers.schedules.ExponentialDecay(
                initial_learning_rate=initial_learning_rate,
                decay_steps=steps_per_epoch,
                decay_rate=learning_rate_decay_factor,
                staircase=True)
optimizer = tf.keras.optimizers.Adam(learning_rate=lr_schedule)
model.compile(loss='categorical_crossentropy',optimizer=optimizer, metrics=['accuracy'])

print(f'{np.isnan(X_train.any())} number of X_train NAs.')
print(f'{np.isnan(X_valid.any())} number of X_valid NAs.')
print(f'{np.isnan(X_test.any())} number of X_test NAs.')

print(f'{np.all(np.isfinite(X_train))} number of X_train finite.')
print(f'{np.all(np.isfinite(X_valid))} number of X_valid finite.')
print(f'{np.all(np.isfinite(X_test))} number of X_test finite.')

print(f'{len(np.isfinite(X_train))/X_train.shape[0]} number of X_train finite.')
print(f'{len(np.isfinite(X_valid))/X_valid.shape[0]} number of X_valid finite.')
print(f'{len(np.isfinite(X_test))/X_test.shape[0]} number of X_test finite.')

print(f'{len(np.isnan(X_train))/X_train.shape[0]} number of X_train NAs.')
print(f'{len(np.isnan(X_valid))/X_valid.shape[0]} number of X_valid NAs.')
print(f'{len(np.isnan(X_test))/X_test.shape[0]} number of X_test NAs.')



history = model.fit(X_train, y_train, validation_data=(X_valid, y_valid), epochs=100, batch_size=32, callbacks=[EarlyStopping('val_accuracy', patience=20)], verbose=2)


# Model prediction
y_predict_prob = model.predict(X_test, batch_size=32, verbose=0)
y_predict_int = np.argmax(y_predict_prob, axis=1)
y_predict_class = label_encoder.inverse_transform(y_predict_int)
y_predict_score = np.amax(y_predict_prob, axis=1)

# Save results
col_cluster = 'neuralNetwork_' + ref_name + '_cluster'
col_score = 'neuralNetwork_' + ref_name + '_clusterScore'

df_map_prob = pd.DataFrame(y_predict_prob, index=cells_test)
df_map_prob.index.name = 'cells'
df_map_prob.columns = ['neuralNetwork_' + ref_name + '_' + label_encoder.inverse_transform([x])[0] for x in df_map_prob.columns]
df_map_prob[col_cluster] = y_predict_class
df_map_prob[col_score] = y_predict_score
df_map_prob.to_csv(output_dir + f'{test_name}_{ref_name}_neuralNetwork_scores.csv', sep=',')
adata_test_raw.obs = pd.concat([adata_test_raw.obs, df_map_prob], axis=1)
print('y_predict_prob saved.')

model.save(f'{output_dir}{test_name}_{ref_name}_neuralNetwork.h5')
print('model saved.')

del y_predict_prob, df_map_prob



# Save model
model.save(output_dir + f'{ref_name}_neuralNetwork.h5')
print('model saved.')

# Plot history
for metric in ['loss', 'accuracy']:
    plt.figure()
    plt.plot(history.history[metric])
    plt.plot(history.history['val_' + metric])
    plt.title(f'{ref_name} neural network {metric} history')
    plt.xlabel('epoch')
    plt.ylabel(metric)
    plt.legend(['train', 'valid'])
    figure = plt.gcf()
    figure.patch.set_facecolor('white')
    figure.savefig(output_dir + f'{ref_name}_neuralNetwork_{metric}_history.png', bbox_inches='tight', dpi=300)
    plt.close('all')
    print(f'{metric} history saved.')

# # Load pretrained model
# label_encoder = load(output_dir + 'label_encoder.joblib')
# model = load_model(output_dir + f'{ref_name}_neuralNetwork.h5')

# **** Neural network model training and prediction finished ****


# **** SHAP evaluation starts if user added -s or --shap at the end of command line ****
if run_shap:
    print('**** SHAP evaluation starts... ****')

    # Select X_test samples with high scores for SHAP, ideally top 50 cells in each class
    df_test_scaled = pd.DataFrame(X_test, index=cells_test, columns=common_genes)
    df_test_scaled.index.name = 'cells'
    df_test_scaled[col_cluster] = y_predict_class
    df_test_scaled[col_score] = y_predict_score

    df_test = []
    for name, group in df_test_scaled.groupby(col_cluster):
        group = group.sort_values(by=col_score, ascending=False).iloc[:50, :]
        if group.shape[0] < 50:
            print('sample size less than 50:', name, group.shape[0])
        df_test.append(group)

    df_test = pd.concat(df_test)
    df_test[[col_cluster, col_score]].to_csv(output_dir + f'{test_name}_cellsForShap.csv', sep=',')
    df_test = df_test.drop(columns=[col_cluster, col_score])
    print('df_test.shape', df_test.shape)

    del X_test, df_test_scaled

    # Summarize X_train that will be used as SHAP background
    df_train = pd.DataFrame(X_train, columns=common_genes)
    df_train = shap.kmeans(df_train, 50)

    del X_train

    # SHAP evaluation
    explainer = shap.KernelExplainer(model.predict, df_train)
    shap_values = explainer.shap_values(df_test)

    # Save SHAP values
    df_shap_summary = []
    for i in range(len(shap_values)):
        label_name = label_encoder.inverse_transform([i])[0]

        df_shap_values = pd.DataFrame(shap_values[i], index=df_test.index.values, columns=common_genes)
        df_shap_values.index.name = 'cells'
        df_shap_values.to_csv(output_dir + f'{test_name}_{ref_name}_neuralNetwork_{label_name}_shapValues.csv', sep=',')

        s_shap_mean = df_shap_values.abs().mean()
        s_shap_mean.index.name = 'Gene'
        s_shap_mean.name = label_name
        df_shap_summary.append(s_shap_mean)

        s_shap_mean.sort_values(ascending=False).to_csv(output_dir + f'{test_name}_{ref_name}_neuralNetwork_{label_name}_shapValues_absMeanAcrossCells.csv', sep=',')

        # SHAP summary plot
        plt.figure()
        shap.initjs()
        shap.summary_plot(shap_values[i], df_test, show=False)
        figure = plt.gcf()
        figure.patch.set_facecolor('white')
        figure.savefig(output_dir + f'{test_name}_{ref_name}_neuralNetwork_{label_name}_shapSummaryPlot.png', bbox_inches='tight', dpi=300)
        plt.close('all')

    df_shap_summary = pd.concat(df_shap_summary, axis=1)
    df_shap_summary.to_csv(output_dir + f'{test_name}_{ref_name}_neuralNetwork_allClusters_shapValues_absMeanAcrossCells.csv', sep=',')
    
    print('**** SHAP evaluation finished. ****')
    print()

else:
    print('**** SHAP evaluation skipped. ****')
    print()

# **** SHAP evaluation finished or skipped ****

# done
print('done')

