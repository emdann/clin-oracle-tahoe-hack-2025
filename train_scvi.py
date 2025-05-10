import numpy as np
import scanpy as sc
import scvi
import torch

torch.set_float32_matmul_precision("high")

adata = sc.read_h5ad("/home/ubuntu/data/merged.h5ad")

scvi.model.LinearSCVI.setup_anndata(
    adata, layer=None, batch_key="batch"
)

model = scvi.model.LinearSCVI(
    adata, 
    n_latent=50,
    latent_distribution="ln",
    use_batch_norm=True,
)

model.train(
    max_epochs=200,
    early_stopping=True,
    early_stopping_patience=2,
    batch_size=16_834,
    load_sparse_tensor=True,
    plan_kwargs={
        "lr": 1e-3,
        "n_epochs_kl_warmup": 2,
    },
    datasplitter_kwargs={
        "pin_memory": False,
    },
    accelerator="gpu",
)

train_elbo = model.history["elbo_train"][1:]
test_elbo = model.history["elbo_validation"]

np.save('train_elbo.npy', np.array(train_elbo))
np.save('test_elbo.npy', np.array(test_elbo))

z_hat = model.get_latent_representation()
np.save('latent_representation.npy', z_hat.cpu().numpy() if isinstance(z_hat, torch.Tensor) else np.array(z_hat))

loadings = model.get_loadings()
loadings.to_csv("loadings.tsv", sep='\t')
