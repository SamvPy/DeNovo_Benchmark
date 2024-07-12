# Set up the de novo pipeline

NOTE: The pipeline uses NextFlow, meaning if you work on Windows, you will need to install WSL and run it there. (https://learn.microsoft.com/en-us/windows/wsl/install)

## Clone the github repository

Run in command line
```
git clone -b master https://github.com/SamvPy/DeNovo_Benchmark.git
```

## Download conda

Download conda from the conda documentation https://docs.anaconda.com/miniconda/ or run:

```bash
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh

~/miniconda3/bin/conda init bash
~/miniconda3/bin/conda init zsh
```

## Create the environments by running the bash script 'create_envs.sh'

1. Go to the source_code directory 

```bash
cd DeNovo_Benchmark/source_code
```

2. (optional if required) Download some required tools that would otherwise cause issues with environment installations:

```bash
sudo apt-get install gcc
```

3. (optional if required) Restart your shell and run the following

```bash
sudo apt-get update
sudo apt-get install build-essential
sudo apt-get install python3-dev
```

4. Create the de novo tool conda environments using the bash script

```bash
chmod +x create_envs.sh
./create_envs.sh
```


# Test the de novo pipeline

1. Running de novo tools

Go towards the nextflow pipelines directory

From the root directory (DeNovo_Benchmark) , run:
```bash
cd nf_pipelines/test
```

2. Update the nextflow.config file

To be able to run the nextflow pipeline, some paths need to be correctly configured.

Open the nextflow.config file in the test directory, and update the following variables
while leaving the others unchanged

- ``workDir``: the directory where nextflow will create folders to run each process in (can be used to check log files when the process files)
- ``config_directory``: All configuration files should be stored here (will most likely be $HOME/DeNovo_Benchmark/configs)
- ``data_directory``: The directory where all raw data and result data is stored. Will need subfolders ``mgf`` and ``results`` if 
not customely configured in parameters ``denovo_results_dir`` and ``mgf_files``
- ``checkpoint_directory``: Directory where all model checkpoints are stored. Model checkpoints are stored in subfolders with the name of the tool in lower case.
- ``params.instanovo_plus_run_script``: Copy the path to the script in 'DeNovo_Benchmark/package/instanovo_scripts/run_instanovoplus.py'

(optional)
- ``params.run``_tool: A boolean flag used to run a given tool
- ``params.serialize``: A boolean flag to parallelize (false) or run every tool in series (true). Can be set to true if mgf-files are large and RAM is limited.

3. Run a test nextflow pipeline to make sure everything works as expected

When you are at the path nf_pipelines/test, and nextflow is installed, run to make predictions for peptide sequences with supported de novo tools:
```bash
nextflow run run_denovo_tools_test.nf
```

Run to rescore/finetune predictions made with the nextflow workflow above:
```bash
nextflow run refine_denovo_results_test.nf
```

If either process fails, the following cases might be the reason for issues: 
- Check if the GPU is recognized by your system
- Check your cuda installation
- When out-of-memory issues, lower the batch_size parameters in the config files of the tools.
- Check whether the paths to the results and workdir directories exist. If not, create them.
- Check whether all conda environments are correctly installed with 
```bash
conda env list
```
You should see the following environments:
- casanovo_env
- instanovo_env
- contranovo_env
- spectralis_env
- pepnet_env
- novob_env
- denovo_analysis_env

If either is missing, create a new environment with this name by going to the correct folder in source_code and run the following:
```bash
conda env create -f "$env_file"
source activate "${env_name,,}_env"

# For casanovo
pip install casanovo

# For spectralis
pip install torch
pip install .

# For all other tools
pip install .
```

4. Run pipelines with your own files!

Make sure you convert your raw data towards mgf files!

Go to the folder nf_pipelines/tools to run de novo tools and get peptide sequence predictions or to nf_pipelines/refinement to post-process (rescore or finetune) de novo outputs.

5. API

If you are experienced in python, you can fully leverage the functionalities of the denovo_utils package.

Open a notebook and select the denovo_analysis_env environment which should have this package already installed.

```python
# To load denovo result files and work with them through the psm-utils package
from denovo_utils.parsers.converters import DenovoEngineConverter

parser = DenovoEngineConverter.select("<insert name of de novo tool here>")
psm_list = parser.parse(
    result_path="<path to de novo output>",
    mgf_path="<path to corresponding mgf-file>"
)
```

TODO: Create documentation for API