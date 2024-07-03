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

Go to the source_code directory and run the bash script

```bash
cd source_code
chmod +x create_envs.sh
./create_envs.sh
```


# Run a de novo pipeline

TODO