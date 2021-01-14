#!/bin/bash

conda update conda -y

cd /application

env_file=$( find /application -name environment.yml )

env_name=$( cat ${env_file} | head -n 1 | cut -d ':' -f 2 | tr -d ' ' ) 

conda env create --file ${env_file}

export PATH=/opt/anaconda/envs/${env_name}/bin:/opt/anaconda/bin:/opt/anaconda/condabin:/usr/local/sbin:/usr/local/bin:/sbin:/bin:/usr/sbin:/usr/bin:$PATH

/opt/anaconda/envs/${env_name}/bin/python -m ipykernel install --name ${env_name}