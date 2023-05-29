# GeneLab utils

Some helper programs for [NASA GeneLab](https://genelab.nasa.gov/), such as `GL-download-GLDS-data` for downloading files from a specific OSD or GLDS ID, and `GL-get-workflow` for downloading workflows used by [GeneLab for processing datasets](https://github.com/nasa/GeneLab_Data_Processing).

---

## Conda install
The genelab-utils package should be installed with conda/mamba. If you are not familiar with conda, you can find an introduction [here](https://astrobiomike.github.io/unix/conda-intro) if wanted, and if you are not familiar with mamba, there is a super-short introduction on that same page [here](https://astrobiomike.github.io/unix/conda-intro#bonus-mamba-no-5) if wanted – it's definitely worth using mamba if you use conda at all :+1: 

```bash
conda install -c conda-forge -n base mamba

mamba create -n genelab-utils -c conda-forge -c bioconda -c defaults -c astrobiomike genelab-utils

conda activate genelab-utils
```

All programs are prefixed with `GL-` and have a help menu accessible with `-h`. Version info can be accessed with `GL-version`.

---

## Some example pages
- Programmatically downloading [GLDS data](https://genelab-data.ndc.nasa.gov/genelab/)
  - [`GL-download-GLDS-data`](https://hackmd.io/@astrobiomike/using-genelab-utils-to-download-GLDS-data)  
- Downloading GeneLab workflows
  - [`GL-get-workflow`](https://hackmd.io/@astrobiomike/using-genelab-utils-to-download-workflows)  

---
---
