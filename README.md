[![Python package](https://github.com/Redx-Pharma/mmp_os/actions/workflows/testing.yaml/badge.svg)](https://github.com/Redx-Pharma/mmp_os/actions/workflows/testing.yaml)

[![Deploy static content to Pages](https://github.com/Redx-Pharma/mmp_os/actions/workflows/static.yaml/badge.svg)](https://github.com/Redx-Pharma/mmp_os/actions/workflows/static.yaml)

# MMP
Repository of methods for Matched Molecular Pairs (MMP) analysis using the open source MMPDB library. Please see [pages](https://reimagined-couscous-wgo85em.pages.github.io/) for more detail on the functions.

## Description
This repository contains codes to enable the MMPDB library to be used through a python interface. The code build the commands for the library and launches tasks as subprocesses to offload to the library.

The main function for the MMP analysis are contained within `src/MMP/runmmp.py`. We include a data prepartion module `prep.py`, a utilities module `utils.py` and a visualization module `viz.py`.

## Dependencies
* Python >3.10
* [MMPDB](https://github.com/rdkit/mmpdb) version 3.1
    * Git clone the above repo
    * Build a virtual env either:
        * `conda create -n mmp python=3.8`
        * `virtualenv venv -p 3.8`
    * Activate the virual env
        * `conda activate mmp`
        * `source /mmp/bin/activate` or on windows `source /mmp/scripts/activate`
    * Then install the MMPDB from the cloned directories top level directory using `pip install .`

## Install

This code is built as a python package which is locally installable.

Editable mode: Run the following from the top level directory
```bash
pip install -e .
```

Fixed mode: Run the following from the top level directory
```bash
pip install .
```

## Examples

#### Prepare data and define MMP rules
The first cell below provides a wrapper to prepare data for MMP. This will process the input csv file given above in the following way:
* Extract and validate the smiles strings then standardize the smiles format for use with RDKit
    * The standardized and validated smiles are stored in the file `$OUTPUT_FILENAME.smi`
    * These smiles and labels are added to the MMP rule database `$OUTPUT_FILENAME.mmpdb`
* Extract the properties from the dataset and store them
    * The properties are stored in the file `properties_$OUPUT_FILENAME.txt`
    * The properties are added to the MMP rule database `$OUTPUT_FILENAME.mmpdb`
    * Statistics are generated over the matched pairs for the change in each of the properties i.e. average and standard deviation over the matched pairs for the change in each property

```python
import numpy as np
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import PandasTools
import MMP
from MMP import runmmp, viz
import mmpdblib

input_data_csv = "/path/to/data.csv"
properties = ["solubility", "melting point", "logP", "logD"]
smiles_column = "SMILES"
labels_column = "Names"
output_file_name = "out_data"

runmmp.setup_mmp(
    input_data_csv,
    properties_columns=properties,
    smiles_column=smiles_column,
    labels_column=labels_column,
    output_filename=output_file_name,
    number_of_cuts=1
    )
```

#### Build New Molecules and Predict Props

```python
base = "c1ccccc1C"
base_mol = Chem.MolFromSmiles(base)

output = runmmp.generate_and_predict(
    base,
    "out_data.mmpdb",
    min_pairs=2,
    min_constant_size=10,
    max_variable_size=7,
    min_radius=1,
    properties_list="solubility,melting point,logP"
    )

vis_df = viz.show_mmp_diffs(base, output, save_filename="out_data_viz.xlsx")
```


## Help/FAQ

* The paper describing the MMPDB library can be found [online](https://pubs.acs.org/doi/10.1021/acs.jcim.8b00173_). This library is built as a wrapper around the MMPDB methods.

## Authors

* James L. McDonagh
