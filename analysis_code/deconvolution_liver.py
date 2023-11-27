#!/Python-3.9.1/bin/python3.9
#SBATCH -J MASH_deconvolution
#SBATCH -o MASH_deconvolution.log
#SBATCH -e MASH_deconvolution.err
#SBATCH --mem=70G
#SBATCH --partition=gpu  --gpus=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=120:00:00
# -*- coding: utf-8 -*-

import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib as mpl

import diopy
import cell2location
import os
os.environ["THEANO_FLAGS"] = 'device=cuda,floatX=float32,force_device=True'
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs


results_folder = path_results #plots
data_folder = path_data

# create paths and names to results folders for reference regression and cell2location models
ref_run_name = f'{data_folder}/reference_signatures'
run_name = f'{data_folder}/cell2location_map'


try: os.mkdir(data_folder)
except: print("Directory already exists")
os.chdir(data_folder)

get_new_model = False
if get_new_model:
    #read REF DATA after have used dior in R to convert it from seurat
    ref_path = path_sn_liver

    adata_ref = diopy.input.read_h5(file = ref_path)

    adata_ref.obs["cluster_annotation"] #this would be the annotation
    adata_ref.var['SYMBOL'] = adata_ref.var_names #set gene as symbol
    # adata_ref.X = adata_ref.layers["counts"] #set non-normalized data

    # adata_ref.var.set_index('index', drop=True, inplace=True)
    
    selected = cell2location.utils.filtering.filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

    # filter the object
    adata_ref = adata_ref[:, selected].copy()

    cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                                # 10X reaction / sample / batch
                                batch_key='library_id', 
                                # cell type, covariate used for constructing signatures
                                labels_key='cluster_annotation',
                                # multiplicative technical effects (platform, 3' vs 5', donor effect)
                                categorical_covariate_keys=None
                            )


    # create the regression model
    mod = cell2location.models.RegressionModel(adata_ref)

    mod.train(max_epochs=500, batch_size = 2048, use_gpu=True)
    # mod.view_anndata_setup() #see summary of the ncells, annotation and so on
    plt.figure()
    mod.plot_history(20)
    plt.legend(labels=['Reference']);
    plt.savefig(f"{results_folder}/loss_ref.png")
    
    adata_ref = mod.export_posterior( #summarise posterior distribution, export cell abundance to anndata object
        adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2048, 'use_gpu': True}
    )
    
    plt.figure()
    mod.plot_QC()
    plt.title("QC ref")
    plt.savefig(f"{results_folder}/QC_ref.png")
    
    # Save model
    mod.save(f"{ref_run_name}", overwrite=True)

    # Save anndata object with results
    adata_file = f"{ref_run_name}/sc_ref_MASH.h5ad"
    adata_ref.write(adata_file)
    print("SAVING MODEL REF")
    
else:
    print(f"loading model ref")
    adata_file = f"{ref_run_name}/sc_ref_MASH.h5ad"
    adata_ref = sc.read_h5ad(adata_file)
    mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", adata_ref)


# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']
print(f"Estimated expression per cluster{inf_aver.iloc[0:5,:]}")

#READ ST
ST_path = path_ST_liver

ST_adata = diopy.input.read_h5(file = ST_path, assay_name="SCT")

ST_adata.var['SYMBOL'] = ST_adata.var_names #set gene as symbol
ST_adata.X = ST_adata.layers["counts"] 

# find mitochondria-encoded (MT) genes
ST_adata.var['MT_gene'] = [gene.startswith('MT') for gene in ST_adata.var['SYMBOL']]
# remove MT genes for spatial mapping (keeping their counts in the object)
ST_adata.obsm['MT'] = ST_adata[:, ST_adata.var['MT_gene'].values].X.toarray()
ST_adata = ST_adata[:, ~ST_adata.var['MT_gene'].values]

dup = ST_adata.var.index.duplicated(keep='first')
ST_adata.var['dup'] = dup
ST_adata = ST_adata[:, ~ST_adata.var['dup'].values]

# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(ST_adata.var_names, inf_aver.index)
print('Intersect with sc:', len(intersect), "genes")
ST_adata = ST_adata[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=ST_adata, batch_key="orig.ident")

mod = cell2location.models.Cell2location(
    ST_adata, cell_state_df=inf_aver, 
    # the expected average cell abundance: tissue-dependent 
    N_cells_per_location=10,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=20
) 
mod.view_anndata_setup()

mod.train(max_epochs=20000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1,
          use_gpu=True,
         )

print("Training done")

# plot ELBO loss history during training, removing first 100 epochs from the plot
plt.figure()
mod.plot_history(1000)
plt.legend(labels=['full ST data training']);
plt.savefig(f"{results_folder}/loss_ST.png")

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
ST_adata = mod.export_posterior(
    ST_adata, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
)

# Save model
mod.save(f"{run_name}", overwrite=True)
print("model ST saved")

# Save anndata object with results
adata_file = f"{run_name}/sp_ST_MASH.h5ad"
ST_adata.write(adata_file)
print("ST object saved")

#load them again
# adata_file = f"{run_name}/sp_ST_CKD.h5ad"
# ST_adata = sc.read_h5ad(adata_file)
# mod = cell2location.models.Cell2location.load(f"{run_name}", ST_adata)

plt.figure()
mod.plot_QC()
plt.title("Quality control ST")
plt.savefig(f"{results_folder}/QC_ST.png")


import pandas as pd
rint("\ncreating csv\n")
# add 5% quantile, representing confident cell abundance, 'at least this amount is present', to adata.obs with nice names for plotting
results = pd.DataFrame(ST_adata.obsm['q05_cell_abundance_w_sf'])
results.columns = ST_adata.uns['mod']['factor_names']

# normalize rows from 0 to 1, convert to percentage
results = results.div(results.sum(axis=1), axis=0)
results.to_csv(os.path.join(output_folder, 'cell_abundance_q5.csv'), index=True)


# Find duplicate gene names in spatial data
from collections import Counter
mylist = ST_adata.var_names
print(f"repeat genes {[k for k,v in Counter(mylist).items() if v>1]}")
