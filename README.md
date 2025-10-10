Preprint: https://www.biorxiv.org/content/10.1101/2025.08.19.671052v1


# Set up the de novo pipeline

NOTE: The pipeline uses NextFlow, meaning if you work on Windows, you will need to install WSL and run it there. (https://learn.microsoft.com/en-us/windows/wsl/install)

## Clone the github repository

Run in command line
```
git clone -b master https://github.com/SamvPy/DeNovo_Benchmark.git
```

## Download conda if not yet done so

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

### 1. (optional if required) Download some required tools that would otherwise cause issues with environment installations:

```bash
sudo apt-get install gcc
```

### 2. (optional if required) Restart your shell and run the following

```bash
sudo apt-get update
sudo apt-get install build-essential
sudo apt-get install python3-dev
```

### 3. Create the de novo tool conda environments using the bash script

```bash
chmod +x create_envs.sh
./create_envs.sh
```


# Test the de novo pipeline

### 1. Running de novo tools

Go towards the nextflow pipelines directory

From the root directory (DeNovo_Benchmark), run:
```bash
cd nf_pipelines/tools
```

### 2. Update the nextflow.config file

To be able to run the nextflow pipeline, some paths need to be correctly configured.

Open the nextflow.config file in the tools directory, and update the following variables
while leaving the others unchanged

- ``workDir``: the directory where nextflow will create folders to run each process in (can be used to check log files when the process fails)
- ``root_results``: Directory where the de novo sequencing results will be stored. Will create another subdirectory called ``denovo_output``
- ``params.mgf_files``: Directory where mgf-files are stored. Example: ``$HOME/DeNovo_Benchmark/test_data/mgf/*.mgf``
- ``params.serializer``: Path to an empty file used to run the de novo tools sequentially. Might be preferred when you run into 'out-of-memory' issues.
- ``params.serialize``: A boolean flag to parallelize (false) or run every tool in series (true). Can be set to true if mgf-files are large and RAM is limited.

(optional)
- ``params.run``_tool: A boolean flag used to run a given tool

Now, open the nextflow.config file in the tools subfolder and update the following variables while leaving the others unchanged

- ``root_checkpoints``: Directory where all model checkpoints are stored. Model checkpoints are stored in subfolders with the name of the tool in lower case.
- ``root_configs``: All configuration files should be stored here (will most likely be $HOME/DeNovo_Benchmark/configs)

(optional)
- ``params.maxforks_tool``: To set how many files should be processed in parallel for each tool. Lower when running into out-of-memory issues.


### 3. Run a test nextflow pipeline

Change the params.mgf_files in ./nf_pipelines/tools/nextflow.config to ``(absolute_prefix)/DeNovo_Benchmark/test_data/mgf/*.mgf

When you are at the path nf_pipelines/tools, and nextflow is installed, run to make predictions for peptide sequences with supported de novo tools:
```bash
nextflow run run_tools.nf
```

Eventhough this mgf only contains 5000 spectra, this can take quiet some time (+-15 minutes per tool depending on the tool on an NVIDIA GeForce RTX 4090). To speed up, you can turn serialization off yet requires much more memory.

To test whether the pipeline works as expected more quickly, change the ``params.mgf_files`` to ``$HOME/DeNovo_benchmark/test_data/test_nf/*.mgf``. Note that on these results, the rescoring pipeline cannot be tested, as there are too few PSMs.

### 4. Refine the de novo results

Make sure you make the appropriate changes to the following files:

- ``$HOME/DeNovo_benchmark/configs/spectralis/spectralis_config.yaml``: Here, you need to specify the path to the spectralis models (binreclass_model_path and scorer_path)

- ``$HOME/DeNovo_benchmark/nf_pipelines/refinement/nextflow.config``: Make sure all path variables point to the correct locations.

Now go inside the nf_pipelines/refinement folder.

Run to rescore/finetune predictions made with the nextflow workflow above:
```bash
nextflow run refine_denovo_results_test.nf
```

If either process fails, the following cases might be the reason for issues: 
- Check if the GPU is recognized by your system
- Check your cuda installation
- When out-of-memory issues, lower the batch_size parameters in the config files of the tools.
- Check whether the paths to the results and workdir directories exist. If not, create them (add the folder ``denovo_output`` to the folder specified in ``root_results``).
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
- pihelixnovo_env
- piprimenovo_env
- adanovo_env
- denovo_analysis_env

### 5. Rescore the results using MS2Rescore

Go to the nf_pipelines/rescoring folder and edit the nextflow.config file so that both the paths for the ms2rescore_script (in Denovo_Benchmark/scripts/rescoring_script.py) and the ms2rescore config have the correct paths. Also edit the ms2rescore config at Denovo_Benchmark/configs/rescoring/ms2rescore so that all the paths are correct and the appropriate de novo engines and post-processors are selected. Also provide the folder where the ground-truth results are stored and specify which file type this has with ground_truth_engine. All formats supported by ``psm-utils`` are supported as well as any de novo engine supported in the ``denovo_utils`` package.

### 6. Explore the results

In the notebooks folder, all notebook files can be explored that generated the figures for the manuscript. Start with the notebooks in the demo folder to create a Run object for you to explore.

### 7. Run pipelines with your own files!

Make sure you convert your raw data towards mgf files!

Go to the folder nf_pipelines/tools to run de novo tools and get peptide sequence predictions or to nf_pipelines/refinement to post-process (rescore or finetune) de novo outputs. After acquiring these results, you can run MS2Rescore-based rescoring and load them into Run objects for in-depth de novo benchmarking.

### 8. API

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
