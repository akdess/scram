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
from sklearn.preprocessing import MinMaxScaler, LabelEncoder
from tensorflow.keras.layers import Input, Dense, Dropout, BatchNormalization, Activation
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.models import Model
from tensorflow.keras.utils import to_categorical
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.wrappers.scikit_learn import KerasClassifier
from joblib import dump
import shap
import matplotlib.pyplot as plt
import argparse
import scanpy as sc
from scipy.sparse import issparse
from sklearn.model_selection import GridSearchCV, StratifiedKFold
import time

# from keras.models import load_model # Only needed when loading pretrained neural network model 
# from joblib import load # Only needed when loading fitted labelEncoder 

print('Packages imported.')
print()


# Set logging
# os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
# tf.logging.set_verbosity(tf.logging.ERROR)


# Set random seeds for reproducibility
os.environ['PYTHONHASHSEED'] = '0'
random.seed(0)
np.random.seed(0)
tf.random.set_seed(0)
# session_conf = tf.ConfigProto(intra_op_parallelism_threads=1, inter_op_parallelism_threads=1)
# sess = tf.Session(graph=tf.get_default_graph(), config=session_conf)
# K.set_session(sess)


# Define nn_classifier
def create_model(activation, hidden_layers, dropout_rate, learning_rate):
    # Get data dimension
    input_shape = (X_train.shape[1],)
    output_size = num_classes

    # Input layer
    inputs = Input(shape=input_shape)
    x = inputs
    
    # Hidden layers
    for layer_size in hidden_layers[:-1]:
        x = Dense(layer_size)(x)
        x = BatchNormalization()(x)
        x = Activation(activation)(x)
        x = Dropout(dropout_rate)(x)
    
    # Last hidden layer for feature extraction
    x = Dense(hidden_layers[-1], name='feature_extraction')(x)
    x = BatchNormalization()(x)
    x = Activation(activation)(x)
    x = Dropout(dropout_rate)(x)
    
    # Output layer
    outputs = Dense(output_size, activation='softmax')(x)
    
    # Build model
    model = Model(inputs=inputs, outputs=outputs)
    
    # Compile model
    optimizer = Adam(learning_rate=learning_rate)
    model.compile(loss='categorical_crossentropy', optimizer=optimizer, metrics=['accuracy'])
    
    return model


def decay_learning_rate(initial_learning_rate, final_learning_rate, epochs, batch_size):
    decay_steps = X_train.shape[0] // batch_size
    decay_rate = (final_learning_rate / initial_learning_rate) ** (1 / epochs)
    lr_schedule = tf.keras.optimizers.schedules.ExponentialDecay(
                    initial_learning_rate=initial_learning_rate,
                    decay_steps=decay_steps,
                    decay_rate=decay_rate,
                    staircase=True)
    return lr_schedule


def recall_learning_rate(lr_schedule, epochs):
    if isinstance(lr_schedule, tf.keras.optimizers.schedules.ExponentialDecay):
        initial_learning_rate = lr_schedule.initial_learning_rate
        decay_rate = lr_schedule.decay_rate
        final_learning_rate = (decay_rate ** epochs) * initial_learning_rate
    else:
        initial_learning_rate = lr_schedule
        final_learning_rate = lr_schedule
    return initial_learning_rate, final_learning_rate


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

# Optional on/off flag to run SHAP evaluation
parser.add_argument('-s', '--shap', action='store_true', help='Add -s or --shap if you opt to run SHAP, otherwise SHAP will not run. e.g., -s')

# Optional on/off flag to normalize reference dataset
parser.add_argument('-r', '--normalize_ref', action='store_true', help='Add -r or --normalize_ref if you opt to normalize reference dataset, otherwise it will not be normalized. e.g., -r')

# Optional on/off flag to normalize testing dataset
parser.add_argument('-t', '--normalize_test', action='store_true', help='Add -t or --normalize_test if you opt to normalize testing dataset, otherwise it will not be normalized. e.g., -t')

# **** User input finished ****


# **** Parse input arguments ****
args = parser.parse_args()

# Arguments for reference dataset
ref_name = args.ref_name
ref_directory = args.ref_directory if args.ref_directory.endswith('/') else args.ref_directory + '/'
ref_class_col = args.ref_class_col

# Arguments for test dataset
test_name = args.test_name
test_directory = args.test_directory if args.test_directory.endswith('/') else args.test_directory + '/'

# Argument for reference marker genes 
marker = args.marker

# Argument to run SHAP
run_shap = args.shap 

# Argument to normalize reference dataset
normalize_ref = args.normalize_ref

# Argument to normalize testing dataset
normalize_test = args.normalize_test

# **** Input arguments parsed ****


# **** Load all datasets according to user input ****
print('**** Loading all datasets according to user input... ****')
print()

# ----------------------
# Load reference dataset
# ----------------------
print(f'Loading {ref_name} as reference dataset...')

# Assume reference directory directly contains adata.h5ad OR { genes.tsv, barcodes.tsv, matrix.mtx }
filename = 'adata.h5ad'
if os.path.isfile(ref_directory + filename):
    print(f'Loading {filename}...')
    adata_ref = sc.read_h5ad(ref_directory + filename)
    print(f'{filename} loaded.')

else:
    print('Loading genes.tsv, barcodes.tsv, and matrix.mtx...')
    for filename in ['genes.tsv', 'barcodes.tsv', 'matrix.mtx']:
        assert os.path.isfile(ref_directory + filename), \
            f'FileNotFoundError: {filename} not found in {ref_directory}'

    with open(ref_directory + 'genes.tsv') as f:
        genes = f.read().rstrip().split('\n')

    with open(ref_directory + 'barcodes.tsv') as f:
        barcodes = f.read().rstrip().split('\n')

    mat = mmread(ref_directory + 'matrix.mtx')
    df = pd.DataFrame.sparse.from_spmatrix(mat, index=genes, columns=barcodes).fillna(0)
    adata_ref = sc.AnnData(df.T)
    del genes, barcodes, mat, df
    print('genes.tsv, barcodes.tsv, and matrix.mtx loaded.')

    # Load reference metadata
    # **** Assume cell ID column name is CellID, otherwise please re-assign it here ****
    index_col = 'CellID'

    # Assume the file is in reference directory, filename is expected as meta.tsv, meta.csv, or meta.xlsx
    # If filename is meta.tsv, assume sep='\t', else if filename is meta.csv, assume sep=',' 
    if os.path.isfile(ref_directory + 'meta.tsv'):
        df_ref_meta = pd.read_csv(ref_directory + 'meta.tsv', sep='\t', index_col=index_col) 
    elif os.path.isfile(ref_directory + 'meta.csv'): 
        df_ref_meta = pd.read_csv(ref_directory + 'meta.csv', sep=',', index_col=index_col)
    elif os.path.isfile(ref_directory + 'meta.xlsx'):
        df_ref_meta = pd.read_excel(ref_directory + 'meta.xlsx', index_col=index_col)
    else:
        df_ref_meta = None
    assert df_ref_meta, \
        f'FileNotFoundError: None of meta.tsv, meta.csv, or meta.xlsx found in {ref_directory}'
    
    # Attach metadata to adata
    df_ref_meta = df_ref_meta.loc[adata_ref.obs_names, :]
    adata_ref.obs = df_ref_meta
    del df_ref_meta

print(f'{ref_name} loaded as reference dataset.')
print('--------')

# In metadata cell type strings, replace each '/' with '.' as they will appear in output filenames
assert ref_class_col in adata_ref.obs_keys(), \
    f'IndexNotFoundError: {ref_class_col} column not found in {ref_name} metadata'
adata_ref.obs[ref_class_col] = adata_ref.obs[ref_class_col].map(lambda x: x.replace('/', '.'))


# -----------------
# Load test dataset
# -----------------
print(f'Loading {test_name} as test dataset...')

# Assume test directory directly contains adata.h5ad OR { genes.tsv, barcodes.tsv, matrix.mtx }
filename = 'adata.h5ad'
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
    print('genes.tsv, barcodes.tsv, and matrix.mtx loaded.')

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
    assert os.path.isfile(ref_directory + marker), \
        f'FileNotFoundError: {marker} not found in {ref_directory}'
    
    # If filename ends with .tsv, assume sep='\t', else if filename ends with .csv, assume sep=',' 
    if marker.endswith('.tsv'):
        df_ref_markers = pd.read_csv(ref_directory + marker, sep='\t')
    elif marker.endswith('.csv'):
        df_ref_markers = pd.read_csv(ref_directory + marker, sep=',')
    else:
        df_ref_markers = pd.read_excel(ref_directory + marker)

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
print('ref_n_samples', adata_ref.n_obs, 'ref_n_genes',  adata_ref.n_vars, 'ref_minValue', adata_ref.X.min(), 'ref_maxValue', adata_ref.X.max())
print('test_n_samples', adata_test.n_obs, 'test_n_genes', adata_test.n_vars, 'test_minValue', adata_test.X.min(), 'test_maxValue', adata_test.X.max())

# Create a directory for output
output_dir = f'{test_name}_{ref_name}_neuralNetwork/'
os.mkdir(output_dir)

# Preprocessing 
print('Preprocessing adata_ref...')
if normalize_ref:
    sc.pp.normalize_total(adata_ref)
    sc.pp.log1p(adata_ref)
sc.pp.highly_variable_genes(adata_ref, n_top_genes=2000)
sc.pp.scale(adata_ref, zero_center=True, max_value=10)

print('Preprocessing adata_test...')
if normalize_test:
    sc.pp.normalize_total(adata_test)
    sc.pp.log1p(adata_test)
sc.pp.scale(adata_test, zero_center=True, max_value=10)

# Select common genes that overlap among training data, testing data, and (optional) training marker genes
markers = set(adata_ref.var_names[adata_ref.var.highly_variable])
if df_ref_markers is not None:
    markers = markers | set(df_ref_markers.index)
    del df_ref_markers

common_genes = markers & set(adata_ref.var_names) & set(adata_test.var_names)
common_genes = sorted(list(common_genes))
with open(output_dir + 'common_genes.tsv', 'w') as f:
    f.write('\n'.join(common_genes) + '\n')
    print('common_genes saved.')

adata_ref = adata_ref[:, common_genes]
adata_test = adata_test[:, common_genes]

# Split adata_ref into 80% training and 20% validation
print('Split adata_ref into 80% training and 20% validation...')
cells_train, cells_valid = [], []
counts = adata_ref.obs[ref_class_col].value_counts()
num_subsample = counts[counts >= 20].quantile(0.25, interpolation='lower')
for celltype, df in adata_ref.obs.groupby(ref_class_col):
    if df.shape[0] < 20:
        print(f'**** Note: cell type {celltype} is excluded for too few samples. ****')
        time.sleep(1)
        continue
    
    cells_this_type = df.index.tolist()
    np.random.shuffle(cells_this_type)

    # Subsample to same level
    cells_this_type = cells_this_type[:num_subsample]

    i = int(len(cells_this_type) * 0.8)
    cells_train += cells_this_type[:i]
    cells_valid += cells_this_type[i:]

# Shuffle the sample order to mix up cell types
np.random.shuffle(cells_train)
np.random.shuffle(cells_valid)

adata_train = adata_ref[cells_train, :]
adata_valid = adata_ref[cells_valid, :]
del adata_ref

# One-hot encoding of celltypes in training data
label_encoder = LabelEncoder()
y_train_int = label_encoder.fit_transform(adata_train.obs[ref_class_col].values)
dump(label_encoder, output_dir + 'label_encoder.joblib')
num_classes = len(np.unique(y_train_int))
y_train_one_hot = to_categorical(y_train_int, num_classes=num_classes)

# One-hot encoding of celltypes in validation data
y_valid_int = label_encoder.transform(adata_valid.obs[ref_class_col].values)
y_valid_one_hot = to_categorical(y_valid_int, num_classes=num_classes)

# Scale each gene in training data
scaler = MinMaxScaler()
X_train = adata_train.X
if issparse(X_train):
    X_train = X_train.toarray()
X_train = scaler.fit_transform(X_train)
dump(scaler, output_dir + 'scaler.joblib')
del adata_train

# Scale each gene in validation data
X_valid = adata_valid.X
if issparse(X_valid):
    X_valid = X_valid.toarray()
X_valid = scaler.transform(X_valid)
del adata_valid

# Scale each gene in testing data
# scaler = load(output_dir + 'scaler.joblib')
X_test = adata_test.X
if issparse(X_test):
    X_test = X_test.toarray()
X_test = scaler.transform(X_test)
cells_test = adata_test.obs_names
del adata_test

# Define decayed learning rate
epochs = 100
batch_size = 32
lr_schedule1 = decay_learning_rate(initial_learning_rate=0.1, final_learning_rate=0.01, epochs=epochs, batch_size=batch_size)
lr_schedule2 = decay_learning_rate(initial_learning_rate=0.1, final_learning_rate=0.001, epochs=epochs, batch_size=batch_size)
lr_schedule3 = decay_learning_rate(initial_learning_rate=0.01, final_learning_rate=0.001, epochs=epochs, batch_size=batch_size)

# Define hyperparameters to tune
activation = ['relu']
# hidden_layers = [
#     [256], [128], [64], [32], 
#     [256, 128], [256, 64], [256, 32], [128, 64], [128, 32], [64, 32], 
#     [256, 128, 64], [256, 128, 32], [256, 64, 32], [128, 64, 32], 
#     [256, 128, 64, 32]
#     ]
hidden_layers = [[256, 64, 32]]
# dropout_rate = [0.1, 0.2]
dropout_rate = [0.1]
# learning_rate = [0.1, 0.01, 0.001, lr_schedule1, lr_schedule2, lr_schedule3]
learning_rate = [lr_schedule2]

# Create a dictionary of hyperparameters
param_grid = dict(
    activation=activation,
    hidden_layers=hidden_layers,
    dropout_rate=dropout_rate,
    learning_rate=learning_rate
)

# Search for best parameters
model = KerasClassifier(build_fn=create_model, epochs=epochs, batch_size=batch_size, verbose=2)
grid = GridSearchCV(estimator=model, param_grid=param_grid, scoring='f1_macro', n_jobs=-1, cv=StratifiedKFold(n_splits=2), error_score='raise')

# Run the parameter tuning loop
grid.fit(X_train, y_train_int)

# Save the best parameters
initial_learning_rate, final_learning_rate = recall_learning_rate(grid.best_params_['learning_rate'], epochs=epochs)

with open(output_dir + 'best model parameters.txt', 'w') as f:
    message = f'Best: {grid.best_score_:.4f} using {grid.best_params_}' + '\n'
    message += f'initial_learning_rate: {initial_learning_rate}, final_learning_rate: {final_learning_rate}' + '\n'
    f.write(message)
    print(message)

# Use the best parameters to create the best model
best_model = create_model(
    activation=grid.best_params_['activation'],
    hidden_layers=grid.best_params_['hidden_layers'],
    dropout_rate=grid.best_params_['dropout_rate'],
    learning_rate=grid.best_params_['learning_rate']
)

# Train the best model
history = best_model.fit(X_train, y_train_one_hot, 
                         validation_data=(X_valid, y_valid_one_hot), 
                         epochs=epochs, 
                         batch_size=batch_size, 
                         callbacks=[EarlyStopping('val_accuracy', patience=10)], 
                         verbose=2)

# Save the best model
best_model.save(output_dir + f'{ref_name}_neuralNetwork.h5')
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
    for ex in ['png', 'pdf']:
        figure.savefig(output_dir + f'{ref_name}_neuralNetwork_{metric}_history.{ex}', bbox_inches='tight', dpi=300)
    plt.close('all')
    print(f'{metric} history saved.')

# # Load pretrained model
# label_encoder = load(output_dir + 'label_encoder.joblib')
# model = load_model(output_dir + f'{ref_name}_neuralNetwork.h5')

# Model prediction
y_predict_prob = best_model.predict(X_test, batch_size=batch_size, verbose=0)
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
print('y_predict_prob saved.')

# Feature extraction
X_test_feature = Model(inputs=best_model.input, outputs=best_model.get_layer('feature_extraction').output).predict(X_test, batch_size=batch_size, verbose=0)

df_test_feature = pd.DataFrame(X_test_feature, index=cells_test)
df_test_feature.index.name = 'cells'
df_test_feature.columns = [f'feature_{x+1}' for x in df_test_feature.columns]
df_test_feature[col_cluster] = y_predict_class
df_test_feature[col_score] = y_predict_score
df_test_feature.to_csv(output_dir + f'{test_name}_{ref_name}_neuralNetwork_features.csv', sep=',')
print('X_test_feature saved.')

# Add to test metadata
# adata_test_raw.obs = pd.concat([adata_test_raw.obs, df_predict_prob, df_test_feature], axis=1)
adata_test_raw.obs = pd.concat([adata_test_raw.obs, df_predict_prob], axis=1)
adata_test_raw.obs = adata_test_raw.obs.loc[:, ~adata_test_raw.obs.columns.duplicated()]
del adata_test_raw.raw # necessary to avoid bugs, otherwise adata can't be saved
# adata_test_raw.write_h5ad(output_dir + 'adata_scores_features.h5ad')
adata_test_raw.write_h5ad(output_dir + 'adata_scores.h5ad')

del y_predict_prob, df_predict_prob, X_test_feature, df_test_feature, adata_test_raw

# **** Neural network model training and prediction finished ****


# **** SHAP evaluation starts if user added -s or --shap at the end of command line ****
if not run_shap:
    print('**** SHAP evaluation skipped. ****')
    print()

else:
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
    explainer = shap.KernelExplainer(best_model.predict, df_train)
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
        for ex in ['png', 'pdf']:
            figure.savefig(output_dir + f'{test_name}_{ref_name}_neuralNetwork_{label_name}_shapSummaryPlot.{ex}', bbox_inches='tight', dpi=300)
        plt.close('all')

    df_shap_summary = pd.concat(df_shap_summary, axis=1)
    df_shap_summary.to_csv(output_dir + f'{test_name}_{ref_name}_neuralNetwork_allClusters_shapValues_absMeanAcrossCells.csv', sep=',')
    
    print('**** SHAP evaluation finished. ****')
    print()

# **** SHAP evaluation finished or skipped ****

# done
print('done')
print()
