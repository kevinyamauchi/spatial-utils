# spatial-utils
a collection of utilities to support analysis with [squidpy](https://github.com/theislab/squidpy).

# installation
## easy installation
You can use our environment.yaml file to create an anaconda environment with spatial-utils, squidpy, and all other dependencies. To do so, first download the environment.yaml file from this repo. Then, open a terminal, navigate to the directory you have downloaded the `environment.yaml` file to and run the following command

```bash
conda env create -f environment.yaml
```

This creates an environment called spatial-analysis. If you would like to use the environment, you can activate the environment with the following:

```bash
conda activate spatial-analysis
```

## development installation
Create an environment for working on your squidpy projects. If you are using anaconda, you can enter the following in your terminal:

```bash
conda create -n squidpy python=3.8
```

Activate your newly created squidpy environment. For anaconda you can enter:

```bash
conda activate squidpy
```

Navigate to the directory you would like to download the squidpy-utils to. For example, if you want to clone the repository to your Documents folder, you would enter:

```bash
cd ~/Documents
```

Clone the squidpy-utils repository

```bash
git clone https://github.com/kevinyamauchi/spatial-utils.git
```

Navigate to the spatial-utils directory

```bash
cd spatial-utils
```

Install `spatial-utiils`. If you would like to be able to add/edit functionality, install in editable mode:

```bash
pip install -e .
```

Otherwise, perform a standard installation:

```bash
pip install .
```

# Usage
### Loading a visium dataset

You can load a dataset that was processed with Kallisto into an AnnData object using the `load_visium_kallisto` function:

```python
from spatial_utils import load_visium_kallisto

adata = load_visium_kallisto(
    counts_table,
    gene_names,
    barcodes,
    tissue_positions_list,
    scale_factors,
    hires_im,
    lowres_im,
    library_id,
    chemistry_name
)
```
The input arguments are explained below. To see an example of the function being used, please see the notebook at `/examples/load_visium_from_kallisto.ipynb`

```
Parameters
----------
counts_table : str
        The path to the counts table output from Kallisto.
        This file usually has the extension ".mtx".
gene_names : str
    The path to the file containing the gene names for the
    Kallisto counts table. This file usually ends with "genes.txt".
barcodes : str
    The path to the file containing the spot barcodes for the
    Kallisto counts table. This file usually ends with "barcodes.txt".
tissue_positions_list : str
    The path to the file containing the coordinates of each barcode
    that is output from the 10X space ranger pipeline.
    This file is usually called: "tissue_positions_list.csv".
scale_factors : str
    The path to the file output by the 10X space ranger pipeline
    containing the scale factors that map the hires and
    lowres image to the original image.
    This file is usually called: "scalefactors_json.json"
hires_im : str
    The path to the hires image that is output from the
    10X space ranger pipeline.
    This file is usually called: tissue_hires_image.png
lowres_im : str
    The path to the lowres image that is output from the
    10X space ranger pipeline.
    This file is usually called: tissue_lowres_image.png
library_id  : str
    The unique identifier for the library that was sequenced.
chemistry_name : str
    The name of the chemistry used to create the library.
    The default value is: "Spatial 3' v1"
```