<div align="center">
<img src="https://github.com/C0nc/River/blob/main/figure/logo.png" width="200px">

**A Python package for identification Differential Spatial Expression Pattern (DESP) gene by interpretable deep learning from multi-slice spatial omics data.**

---

<p align="center">
  <a href="https://www.biorxiv.org/content/10.1101/2024.08.04.606512v1" target="_blank">Preprint</a>
</p>

</div>

River is able to identify Differential Spatial Expression Pattern (DSEP) across multi-slice dataset, and offers the downstream analysis based on obtained DSEP genes.
</p>
<p align="center">
  <img src="https://github.com/C0nc/River/blob/main/figure/pipeline.png" width="800px">
</p>

## Getting started


Please refer to the  
- Stereo-seq 3D dataset [Tutorial][link-tutorial_2] (Can be downloaded by `pysodb` package) 
- Stereo-seq development dataset [Tutorial][link-tutorial_3] (Can be downloaded by `pysodb` package)
- Slide-seq mouse diabetes disease dataset [Tutorial][link-tutorial_4]. (Can be downloaded by `pysodb` package)
- MIBI TNBC disease dataset [Tutorial][link-tutorial_5]. (Can be downloaded by `pysodb` package)
- CODEX lupus dataset [Tutorial][link-tutorial_6]. (Can be downloaded by `pysodb` package)

## Installation

1. Create a conda environment
```bash
conda create -n river python=3.8 -y && conda activate river

```
2. Install the River dependency
```bash
pip install scSLAT
python -c "import torch; print(torch.__version__)"
pip install pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv -f https://data.pyg.org/whl/torch-2.0.0+cu117.html  # replace torch and CUDA version to yours
pip install captum ipykernel 
```

Install the `pysodb` for efficient download processed Anndata in h5ad format (https://pysodb.readthedocs.io/en/latest/) 
Install the `CellCharter` for multi-slice co-clustering in Slide-seq analysis (https://github.com/CSOgroup/cellcharter)
## Contribution

If you found a bug or you want to propose a new feature, please use the [issue tracker][issue-tracker].

[issue-tracker]: https://github.com/C0nc/River/issues
[link-docs]: https://cellcharter.readthedocs.io
[link-api]: https://cellcharter.readthedocs.io/en/latest/api.html
[link-tutorial_1]: https://github.com/C0nc/River/blob/main/figure2.ipynb
[link-tutorial_2]: https://github.com/C0nc/River/blob/main/figure3.ipynb
[link-tutorial_3]: https://github.com/C0nc/River/blob/main/figure4.ipynb
[link-tutorial_4]: https://github.com/C0nc/River/blob/main/figure5.ipynb
[link-tutorial_5]: https://github.com/C0nc/River/blob/main/figure6.ipynb
[link-tutorial_6]: https://github.com/C0nc/River/blob/main/figure7.ipynb
