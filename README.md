<div align="center">
<img src="https://github.com/C0nc/River/blob/main/figure/logo.png" width="200px">

**A Python package for the Scalale and accurate identification condition-relevant niches from spatial omics data.**

---

<p align="center">
  <a href="https://doi.org/10.1101/2024.05.30.596656" target="_blank">Preprint</a>
</p>

</div>

Taichi is able to automatically identify condition-relevant niches, and offers the downstream analysis based on obtained niches.
</p>
<p align="center">
  <img src="https://github.com/C0nc/River/blob/main/figure/pipeline.png" width="800px">
</p>

## Getting started


Please refer to the  
- Stereo-seq 3D dataset [Tutorial][link-tutorial_2] (Can be downloaded by `pysodb` package) 
- Stereo-seq development dataset [Tutorial][link-tutorial_3] (Can be downloaded by `pysodb` package)
- MIBI TNBC disease dataset [Tutorial][link-tutorial_4]. (Can be downloaded by `pysodb` package)
- CODEX lupus dataset [Tutorial][link-tutorial_5]. (Can be downloaded by `pysodb` package)

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

Install the `pysodb` for efficient download processed Anndata in h5ad format (https://pysodb.readthedocs.io/en/latest/) if you want to run the DKD and CRC related analysis

## Contribution

If you found a bug or you want to propose a new feature, please use the [issue tracker][issue-tracker].

[issue-tracker]: https://github.com/C0nc/River/issues
[link-docs]: https://cellcharter.readthedocs.io
[link-api]: https://cellcharter.readthedocs.io/en/latest/api.html
[link-tutorial_1]: https://github.com/C0nc/River/blob/main/figure_2.ipynb
[link-tutorial_2]: https://github.com/C0nc/River/blob/main/figure_3.ipynb
[link-tutorial_3]: https://github.com/C0nc/River/blob/main/figure_4.ipynb
[link-tutorial_4]: https://github.com/C0nc/River/blob/main/figure_5.ipynb
[link-tutorial_5]: https://github.com/C0nc/River/blob/main/figure_6.ipynb
