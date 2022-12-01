

vagrant destroy && rm Vagrantfile && \
export VM=sylabs/singularity-3.0-ubuntu-bionic64 && \
vagrant init $VM && \
vagrant up && \
vagrant ssh


logout
vagrant ssh

conda create --name MedGen && conda activate MedGen
conda install -c bioconda snakemake samtools bedtools
