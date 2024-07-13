# This script loads one testing dataset and multiple pre-trained models for prediction.
# The output adata_scores_concat.h5ad has concatenated scores predicted by multiple models.


import os
import random
import numpy as np
import pandas as pd
import scanpy as sc
import tensorflow as tf
from scipy.io import mmread
from scipy.sparse import issparse
from joblib import load
from keras.models import load_model
from keras import backend as K
from tqdm import tqdm


os.environ['PYTHONHASHSEED'] = '0'
random.seed(0)
np.random.seed(0)
tf.random.set_seed(0)


# --------------------------
# Assign key parameters here
# --------------------------
test_name = '041423_glioma_withscanpyscale'
ref_names = ['suva_idh_a_o_v2', 'hpa_brain_simple', 'allen_class_label_main', 'allen_neurons_only', 'TissueImmune', 'aldinger', 'codex', 'suva', 'bhaduri','dirks_primary_gbm_combined']
normalize_test = True


# -------------------
# Assume folder names
# -------------------
# One folder contains a set of files pre-defined using one reference dataset:
# common_genes.tsv, label_encoder.joblib, scaler.joblib, {ref_name}_neuralNetwork.h5
model_dirs = {ref_name: f'{test_name}_{ref_name}_neuralNetwork/' for ref_name in ref_names}


# -----------------
# Load test dataset
# -----------------
print(f'Loading {test_name} as test dataset...')

# Assume test directory directly contains adata.h5ad OR { genes.tsv, barcodes.tsv, matrix.mtx, meta.csv }
test_directory = f'{test_name}/'
filename = 'glioma.h5ad'

if os.path.isfile(test_directory + filename):
    print(f'Loading {filename}...')
    adata_test = sc.read_h5ad(test_directory + filename)
    print(f'{filename} loaded.')

else:
    print('Loading genes.tsv, barcodes.tsv, and matrix.mtx...')
    for filename in ['genes.tsv', 'barcodes.tsv', 'matrix.mtx']:
        assert os.path.isfile(test_directory + filename), \
            f'FileNotFoundError: {filename} not found in {test_directory}'

    with open(test_directory + 'genes.tsv') as f:
        genes = f.read().rstrip().split('\n')

    with open(test_directory + 'barcodes.tsv') as f:
        barcodes = f.read().rstrip().split('\n')

    mat = mmread(test_directory + 'matrix.mtx')
    df = pd.DataFrame.sparse.from_spmatrix(mat, index=genes, columns=barcodes).fillna(0)
    adata_test = sc.AnnData(df.T)

    df_test_meta = pd.read_csv(test_directory + 'meta.csv', sep=',', index_col=0)
    adata_test.obs = df_test_meta.loc[barcodes, :]
    
    del genes, barcodes, mat, df, df_test_meta
    print('genes.tsv, barcodes.tsv, matrix.mtx, meta.csv loaded.')

adata_test_raw = adata_test.copy()
print(f'{test_name} loaded as test dataset.')
print('--------')


# **** Neural network model loading and prediction starts ****

# Create a directory for output
output_dir = f'{test_name}_{"_".join(ref_names)}_multipleNeuralNetworks/'
os.mkdir(output_dir)

# Preprocessing 
print('Preprocessing adata_test...')
if normalize_test:
    sc.pp.normalize_total(adata_test)
    sc.pp.log1p(adata_test)
sc.pp.scale(adata_test, zero_center=True, max_value=10)

# Prediction using each model
record = []
for ref_name in tqdm(ref_names):
    model_dir = model_dirs[ref_name]

    # Load common genes
    with open(model_dir + 'common_genes.tsv') as f:
        common_genes = f.read().rstrip().split('\n')

    adata = adata_test[:, common_genes]

    # Load pre-trained models
    label_encoder = load(model_dir + 'label_encoder.joblib')
    scaler = load(model_dir + 'scaler.joblib')
    best_model = load_model(model_dir + f'{ref_name}_neuralNetwork.h5')

    # Scale each gene in testing data
    X_test = adata.X
    if issparse(X_test):
        X_test = X_test.toarray()
    X_test = scaler.transform(X_test)
    cells_test = adata.obs_names
    del adata
    
    # Model prediction
    y_predict_prob = best_model.predict(X_test, batch_size=32, verbose=0)
    y_predict_int = np.argmax(y_predict_prob, axis=1)
    y_predict_class = label_encoder.inverse_transform(y_predict_int)
    y_predict_score = np.amax(y_predict_prob, axis=1)

    # Save results
    col_cluster = 'neuralNetwork_' + ref_name + '_cluster'
    col_score = 'neuralNetwork_' + ref_name + '_clusterScore'

    df_predict_prob = pd.DataFrame(y_predict_prob, index=cells_test)
    df_predict_prob.index.name = 'cells'
    df_predict_prob.columns = ['neuralNetwork_' + ref_name + '_' + label_encoder.inverse_transform([x])[0] for x in df_predict_prob.columns]
    df_predict_prob[col_cluster] = y_predict_class
    df_predict_prob[col_score] = y_predict_score
    df_predict_prob.to_csv(output_dir + f'{test_name}_{ref_name}_neuralNetwork_scores.csv', sep=',')

    # Add to record
    record.append(df_predict_prob)

    # Lower memory usage
    del best_model
    K.clear_session()
    

# Output concatenated df
df_predict_prob = pd.concat(record, axis=1)
df_predict_prob.to_csv(output_dir + f'{test_name}_{"_".join(ref_names)}_multipleNeuralNetworks_scores.csv', sep=',')

# Output concatenated adata
adata_test_raw.obs = pd.concat([adata_test_raw.obs, df_predict_prob], axis=1)
adata_test_raw.obs = adata_test_raw.obs.loc[:, ~adata_test_raw.obs.columns.duplicated()]
del adata_test_raw.raw # necessary to avoid bugs, otherwise adata can't be saved
adata_test_raw.write_h5ad(output_dir + 'adata_scores_concat.h5ad')
