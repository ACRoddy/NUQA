# NUQA

NGS tool for Unsupervised analysis of fastQ files using Alignment-free

A software to support an alignment-free sequence comparison of large-scale fastq files using a k-mer counting methods paired with distance metrics, Jensen-Shannon divergence and Hellinger Distance.

## Installation

Dependencies:

  *[jellyfish](https://github.com/gmarcais/Jellyfish/releases) V2.2.6+ 
  *phylip V3.6+
  *gcc V6.5+
  *autoconf V2.69+
  *automake V1.15+
  *Eigen (where /Eigen file containing the headers is located at /usr/local/include for an example run check file 'Eigen_help.txt')


Once all dependencies have been installed you can proceed with the following steps:

```bash
git clone https://github.com/ARoddy/NUQA.git
cd NUQA
./configure
make
make install
```


## Usage

NUQA reads in a set of raw fastq files - ideally, these are DNA-seq files produced from multiple samples taken from the same patient. Based on these it will identify differences between the files and produce a distance matrix and a newick tree file

You can run NUQA using the following command

`bash NUQA.sh -k *int* -t *int* -d *distance* -p *prefix* /path/to/directory/*.fastq`

options:
  k - chosen k-mer length (Default:21)
  t - number of threads (Default:1)
  d - distance metric (j for Jensen-Shannon divergence, h for Hellinger distance, Default:j)
  p - prefix for output files (Default:NUQA)

##License

Released under GNU GPLv3
Copyright 2019 Aideen Roddy


