# This file is the instruction for running the HE2RNA code on Great Lakes
# Note that you probably need to specify the correct path

# Gene expression prediction

Predict gene expression from WSIs taken from TCGA with HE2RNA [1]. The model takes as inputs arrays of size n_tiles * 2048, where n_tiles = 100 when super-tile preprocessing is used, and n_tiles = 8,000 when all tiles are treated separately. The model is implemented as a succession of 1D convolution (equivalent to an MLP shared among all tiles).
Additionally, Model interpretability can be explored at: https://owkin.com/he2rna-result-visualization/.

## Installation

Create a virtual environment and install the required packages (the variable CUDA_TOOLKIT_ROOT_DIR is needed to install libKMcuda):

```bash

# Create environment
conda create -n env 
conda activate env

# Install required packages
cd <turbo drive>/<your name>/HE2RNA_code
pip install -r requirements1.txt

# Find CUDA path
env | grep cuda

# Load module
ml Bioinformatics openslide

# Manually build LibKMCUDA package (see https://github.com/src-d/kmcuda)
git clone https://github.com/src-d/kmcuda
cd kmcuda/src
cmake -D CUDA_ARCH=70 -D CUDA_TOOLKIT_ROOT_DIR=/sw/arcts/centos7/cuda/11.0.2/bin -DCMAKE_BUILD_TYPE=Release . && make

```
Note that the libKMCUDA is in "~/.local/lib/python3.7/site-packages/"

## Data collection and preprocessing

To ensure reproducibility of the results, coordinates of the tiles used in the paper (necessary to extract tile images and features from whole-slide images) are provided in the archive tile_coordinates.gz.

EDIT: due to an issue related to data quota, file tile_coordinates.gz should be downloaded instead from https://drive.google.com/file/d/1PJsUv1SQieJs7hqtWOqW68v1K9c-mIF6/view?usp=sharing.

To uncompress it, run
```bash
tar -xzvf tile_coordinates.gz
```
Splits used in the paper are also provided in patient_splits.pkl.

### TCGA

We originally downloaded the whole-slide images from the TCGA data portal https://portal.gdc.cancer.gov/ via the gdc-client tool. To access all of TCGA data used in this work, follow the steps described below.

Paths to folders containing slides, tile features and RNAseq data should be consistent with the contant of file constant.py. If data is saved in a different location, constant.py has to be modified accordingly, as well as the example config files.

#### Download TCGA slides (if needed)
Go to the dedicated folder to store files (all FFPE slides from TCGA is approx. 10To)
```bash
cd <turbo drive>/shared_data/TCGA-LGG_SVS_Raw
```
Download images using the corresponding manifest:
```bash
gdc-client download -m ../<your name>/HE2RNA_code/gdc_manifests/gdc_manifest.2018-06-26_TCGA-LGG.txt

```

#### Tile feature extraction

The code in extract_tile_features_from_slides.py is designed to extract resnet features of tile images directly from whole-slide images, using the coordinates of the tiles in Openslide format. To extract tile features from WSIs from a given TCGA project, e.g. LGG, run:
```bash
cd <turbo drive>/<your name>
mkdir TCGA_tiles/

python extract_tile_features_from_slides.py --path_to_slides <turbo drive>/shared_data/TCGA-LGG_SVS_Raw --tile_coordinates <turbo drive>/<your name>/HE2RNA_code/tile_coordinates/tile_coordinates_TCGA_LGG.pkl --path_to_save_features <turbo drive>/<your name>/TCGA_tiles/TCGA_LGG

python extract_tile_features_from_slides_ori.py --path_to_slides <turbo drive>/shared_data/TCGA-GBM_SVS_Raw --tile_coordinates <turbo drive>/<your name>/HE2RNA_code/tile_coordinates/tile_coordinates_TCGA_GBM.pkl --path_to_save_features <turbo drive>/<your name>/TCGA_tiles/TCGA_GBM

```
Note that if you download the image using the manifest data in the previous step, please refer to the original github code for tile feature extraction!

#### Download and preprocess RNAseq data
Create a folder to store rnaseq data and download transcriptomes:
```bash
cd <turbo drive>/<your name>
mkdir TCGA_transcriptome
cd TCGA_transcriptome
gdc-client download -m ../HE2RNA_code/gdc_manifests/gdc_manifest.2018-03-13_alltranscriptome.txt
```
At this stage, there should be one folder per sample, containing a .gz archive. Extract the archives, using for instance gunzip
```bash
gunzip */*.txt.gz
```
To make things more convenient, we already save a file containing transcriptomes matched to whole-slide images, using
```bash
cd <turbo drive>/<your name>/HE2RNA_code
python transcriptome_data.py
```

#### Supertile preprocessing (optional)
Finally, once all previous steps have been performed, supertile preprocessing can be performed using the following command (the csv file containing transcriptome is used here to ensure consistency between preprocessed image samples and RNAseq data),
```bash
python supertile_preprocessing.py --path_to_transcriptome <turbo drive>/<your name>/TCGA_transcriptome/all_transcriptomes.csv --path_to_save_processed_data <turbo drive>/<your name>/TCGA_100_supertiles.h5 --n_tiles 100

```

## Gene expression prediction experiment

To run an experiment, write first a config file or use one of the examples available in folder condigs. 

* config_CD3_LGG.ini: prediction of CD3 genes on LGG, using supertiles, and starting training from scratch!

Launch experiment with a single train-test split:
```bash
python main.py --configs/config <config_file> --run single_run --logdir ./exp
```
Launch cross-validation:
```bash
python main.py --config configs/config_CD3_LGG.ini --run cross_validation --n_folds 5 --logdir ./exp_cd3_gbm_lgg
python main.py --config configs/config_EGFR_LGG_GBM.ini --run cross_validation --n_folds 5 --logdir ./exp_egfr_gbm_lgg
```

Binary classification (example):
```bash
python main_cat.py --config configs/config_CD3_LGG_cat.ini --run cross_validation --n_folds 5 --logdir ./exp_cd3d_01
```
Launch TensorboardX for visualizing training curves
```bash
tensorboard --logdir=./exp --port=6006
```

Results will be saved in the specified path as follows:
* for a single train/valid/test split, the model will be saved as model.pt and the correlation per gene and cancer type will be saved as results_single_split.csv
* for a cross-validation, each model will be saved in a dedicated folder model_i/model.pt, the correlation per gene, cancer type and fold will be saved as results_per_fold.csv.

### Config file options

* [main]
	* path: Path to the directory where model's weights will be saved.
	* use_saved_model (optional): Path to previous experiment to reload saved models
	* splits (optional): Path to Pickle file containing saved patient splits for cross-validation, useful in particular when finetuning a model on a subset of the data, to ensure consistency of the train and test set with those used for pretraining.
    * single_split (optional): Path to Pickle file containing saved patient split for single run

* [data]
	* genes (optional): List of coma-separated Ensembl IDs, or path to a pickle file containing such a list. If None, all available genes with nonzero median expression are used.
	* path_to_transcriptome (optional): If None, build targets from projectname and list of genes. Otherwise, load transcriptome data from a saved csv file.
	* path_to_data (optional): Path to the data, saved either in a pickle file (for aggregated data) or in an hdf5 file. If None, build the dataset from .npy files.

* [architecture]
	* layers: Integers defining the number of feature maps of the model's 1D convolutional layers
	* dropout: Float between 0 and 1.
	* ks: List of ks to sample from
	* nonlin: 'relu', 'sigmoid' or 'tanh'.
	* device: 'cpu' or 'cuda'.

* [training]
	* max_epochs: Integer, defaults to 200.
    * patience: Integer, defaults to 20.
    * batch_size: Integer, defaults to 16.
    * num_workers: number of workers used for loading batches, defaults to 0 (value should be 0 when working with hdf5-stored data)

* [optimization]
	* algo: 'sgd' or 'adam'.
	* lr: Float.
	* momentum: Float, optional
    
## Spatialization of gene expression

### Spatialization of lymphocyte genes in LGG

Once a model has been trained to predict the expression of genes specifically expressed by lymphocytes (for instance CD3), the following script can be used to visualize the mRNA expression (take one HE slide for example).

```bash
#LGG
python spatialization.py --path_to_tile_features <turbo drive>/<your name>/TCGA_tiles/TCGA_LGG/0.50_mpp/TCGA-HT-8113-01Z-00-DX1.npy -–path_to_slide <turbo drive>/shared_data/TCGA-LGG_SVS_Raw/TCGA-HT-8113-01Z-00-DX1.641F5405-47CF-41C6-8EA4-4ABD6C677A46.svs 

#GBM
python spatialization.py --path_to_tile_features <turbo drive>/<your name>/TCGA_tiles/TCGA_GBM/0.50_mpp/TCGA-06-0152-01Z-00-DX9.npy -–path_to_slide <turbo drive>/shared_data/TCGA-GBM_SVS_Raw/c8cf5f17-f0ee-4587-a4e8-d3461f374cdc/TCGA-06-0152-01Z-00-DX9.8620b410-9d16-418e-88e1-b75c020da50a.svs
```

## References

[1] Schmauch, B., Romagnoni, A., Pronier, E., Saillard, C., Maillé, P., Calderaro, J., ... & Courtiol, P. (2019). Transcriptomic learning for digital pathology. bioRxiv, 760173.

[2] Kather, J. N et al. 100,000 histological images of human colorectal cancer and healthy tissue (Version v0.1). Zenodo. http://doi.org/10.5281/zenodo.1214456 (2018).

[3] Bulten, W., et al. PESO: Prostate Epithelium Segmentation on H&E-stained prostatectomy whole slide images (Version 1). Zenodo. http://doi.org/10.5281/zenodo.1485967 (2018).

# License

GPL v3.0
