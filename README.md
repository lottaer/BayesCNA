# Sensitive detection of copy number alterations in low-pass liquid biopsy sequencing data

This repository contains the code used for the analysis in _Sensitive detection of copy number alterations in low-pass liquid biopsy sequencing data_. In liquid biopsy data, the sequenced DNA is a mixture of both healthy healthy DNA (without genomic alterations) and tumor DNA (containing genomic alterations), which reduces the signal and makes the detection of CNAs more difficult. Furthermore, the samples are contaminated by sequencing noise that further obscures the signal. We introduce BayesCNA, a method designed to improve signal extraction from low-pass liquid biopsy sequencing data, by utilizing a Bayesian changepoint detection algorithm, to enahnce signal extraction from such samples.

The workflow is presented in `analysis.R`

For further details, see the manuscipt [TBA]

## Details on data simulation

## Generating synthethic datasets

The code for simulating synthethic dataset is presented in `simulate_data.R`. The pipeline for generating and processing the data is presented in `analysis.R`.

### Generating raw sequencing reads

We run CNVsim (available here: https://github.com/NabaviLab/CNV-Sim) using Docker. Pre-defined alterations are defined in the `bed`-format. The script `simulate_bed.R` is used to produce the `bed`-files. 

The simulation pipeline is presented in `make_sample.sh`. To run the simulation file, the absolute path to the folder containing the required genome file (`chr1.fa`) should be provided

```
docker run --rm -v <absolute_local_path_to_input_directory>:/my_data nabavilab/cnv-sim \
         ./cnv-sim.py -o /my_data/$OUTPUT_DIR --read_length 100 --n_reads 0 --cnv_list /my_data/$BED_FILE genome /my_data/$1
```
The human genome reference is available for download here: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/

The script for estimating the purity of the tumor data is presented in `purity_script.R`. A generated sample is discared if the estimated purity is 50% (can be adjusted by changing `LOWER_LIMIT` in `make_sample.sh`). 

Finally, the script `make_sample_multiple.sh` is used a wrapper of `make_sample.sh` to produce multiple samples with purities [5, 10, ..., 30] at a certain coverage. The script takes the following arguments __in this particular order__
- Number of samples to produce (the number of outputted samples might be lower to this if samples are discarded due to too low purity)
- Path to reference file
- Desired coverage of output files
Example of simulating 10 samples of purity [5, 10, ..., 30] of coverage 0.15 and using reference file `chr1.fa` in folder `reference`
```
bash make_sample_multiple 10 reference/chr1.fa 0.15
```

### Downsampling bam-files

The pipeline used for producing the downsampled `bam`-files is presented in the script `downsample_bam.sh`
The script takes the following arguments, __in this particular order__
- Path to the `.bam`-file to subsample
- Million number of reads in the output file
- Name of the output file (without `-bam`-extension)
The scripts expects a folder `cell_line_subs`, where it outputs the downsampled file. Example of downsampling file to 10 million reads with output name `downsampled.bam`

```
sh downsample_bam.sh path/to/bam_file 10 downsampled
```

## Licence

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
