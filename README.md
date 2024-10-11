# MiniSim
A realisitic long reads RNA-seq data simulator <br>
MiniSim utilizes kernel density estimation (KDE) models trained based on real ONT data to estimate 5' and 3' truncation length. Reads simulated by miniSim better represent the 5' and 3' truncation length distribution real ONT data
![Compare_truncation_miniSim and nanosim](https://github.com/Tidesun/MiniSim/blob/a34026cbd907c93008c4c7b7dc693640b8f29b98/Supplementary_Fig_19.png)
## Requirements
### Dependency
```
Python>=3.6
Linux operating system
```
The software has been tested with following software version
```
Python==3.9.7
```
## Installation and set up pretrained models
### Installation
```
git clone https://github.com/Tidesun/MiniSim.git
cd MiniSim
python -m venv base
source base/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```
### Download pre-trained KDE models for the truncation profile
1. `cd MiniSim & source base/bin/activate`
2. Download `models.zip` from https://www.dropbox.com/sh/ju9go788dwq3znr/AACsVfaVXyPZF7YxeqEBq_8Ga?dl=0
3. `unzip models.zip`
## Usage
### Simulate long reads RNA-seq data
```
usage: python simulator.py [-h] --kde_model_path KDE_MODEL_PATH [--insertion_rate INSERTION_RATE] [--deletion_rate DELETION_RATE] [--substitution_rate SUBSTITUTION_RATE] [--coord_error_in_5_end COORD_ERROR_IN_5_END]
                    [--coord_error_in_3_end COORD_ERROR_IN_3_END] -expr EXPRESSION_PROFILE -cdna REFERENCE_TRANSCRIPTOME --num_reads NUM_READS [-t THREADS] -o OUTPUT

Simulate reads with kde model

optional arguments:
  -h, --help            show this help message and exit

required named arguments:
  --kde_model_path KDE_MODEL_PATH
                        Path of KDE model. For example models/kde_cDNA to use directRNA-ONT truncation model.
  --insertion_rate INSERTION_RATE
                        Insertion rate
  --deletion_rate DELETION_RATE
                        Deletion rate
  --substitution_rate SUBSTITUTION_RATE
                        Substitution rate
  --coord_error_in_5_end COORD_ERROR_IN_5_END
                        Coord randomness in 5' end
  --coord_error_in_3_end COORD_ERROR_IN_3_END
                        Coord randomness in 3' end
  -expr EXPRESSION_PROFILE, --expression_profile EXPRESSION_PROFILE
                        Expression profile
  -cdna REFERENCE_TRANSCRIPTOME, --reference_transcriptome REFERENCE_TRANSCRIPTOME
                        Reference transcriptome
  --num_reads NUM_READS
                        The number of simulated reads
  -t THREADS, --threads THREADS
                        Threads
  -o OUTPUT, --output OUTPUT
                        The path of output simulated reads
```
### Train KDE models from real long reads RNA-seq data
```
python build_kde_model/build_kde_model.py [-h] -gtf GTF_ANNOTATION_PATH -o OUTPUT_PATH -lrsam LONG_READ_SAM_PATH -t THREADS

Build kde model on long reads real data

optional arguments:
  -h, --help            show this help message and exit
  -gtf GTF_ANNOTATION_PATH, --gtf_annotation_path GTF_ANNOTATION_PATH
                        The path of annotation file
  -o OUTPUT_PATH, --output_path OUTPUT_PATH
                        The path of output directory
  -lrsam LONG_READ_SAM_PATH, --long_read_sam_path LONG_READ_SAM_PATH
                        The path of long read sam file mapping to reference genome
  -t THREADS, --threads THREADS
                        Threads
```


