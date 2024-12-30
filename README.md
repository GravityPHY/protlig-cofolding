## Benchmarking co-folding algorithms
### Author: Hao Yu
### Email: imhaoyu@bu.edu
This repository contains my summer intern project at Novartis in 2024, which includes all the scripts and instructions I have written. Note that there is no involvement of commercialized content.


#### Infrastructure
##### Install Umol 
Go to `umol_conda` folder, 
```
#load modules
ml purge
ml Conda

#proxy
export HTTP_PROXY=http://nibr-proxy.global.nibr.novartis.net
export HTTPS_PROXY=https://nibr-proxy.global.nibr.novartis.net

conda env create -f environment.yml
```
I named the environment as `umol-test`, feel free to change the name in the `environment.yml` file. 
`conda activate umol-test` should activate the environment.  
Next, activate the enviroment, and installed openmmforcefields. (Installed from conda-forge will give you version 0.13.0, which missing some functions that will be used in the Umol framework.) We can download the specific version source code of openmmforcefields and install to the current conda environment. I have downloaded 0.12.0 version (from https://github.com/openmm/openmmforcefields/archive/refs/tags/0.12.0.tar.gz) and put in this repository, the file name is `openmmforcefields-0.12.0.tar.gz`.  

```console
conda activate umol-test

# this step is to get the path of python executable of the currecnt environment
which python3 # the <absolute path to python3> would show in the next line
<absolute path to python3> -m pip install openmmforcefields-0.12.0.tar.gz
```

Once the installation is complete, the umol conda environment is ready!  
If you want to use this environment in a Jupyter Notebook,
```
conda install ipykernel 
ipython kernel install --user --name=<name_for_kernel>
```
When you open a Jupyter Notebook, select the kernel name.  
Additional notes if version conflicts show up when building the conda environment.


- run ```mamba env create -f environment.yml```  to check what are the conflicts. Mostly the incompatiblity comes from opemmforcefields version, python version, numpy-base(already deleted it in the .yml in this repo), intel-openmp (already deleted it in the .yml in this repo).

##### Umol pipeline 
Prepare the input information in a `.json` file.  
Example
```
{"STRUCID":”1oda”,
"EXPRSEQUENCE":” GPDL…”,
"POCKETPOS":null,
"SMILES":"CC...OC=O",
"OUTPUTDIR":”./”,
"random":42}
```

Run prediction with wrapper function
```python3 .umol_conda/wrapper/predict_umol.py --json_path input.json```

#### Data
single chain single ligand benchmark set available at this [internal link](go/cofold-benchmark).  
Analysis scripts 
- Visualize tool (developming)

#### Results

#### Other notes
You might find this [note](https://confluence.prd.nibr.novartis.net/x/xBf2Ew) useful when trying to create a conda environment working with pytorch_geometric.  
You might find this [note](https://confluence.prd.nibr.novartis.net/display/~yuha1k/Lessons+for+Conda+environment+installation) helpful when creating an environment from .yml file.  
If you are trying to install packages through conda locally at the company computer, check this [note](https://confluence.prd.nibr.novartis.net/display/~yuha1k/Develop+Tools+Conda+and+Git)

