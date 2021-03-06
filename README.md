SeqClone
==============================

seqclone identifies genes in single cell RNA-seq data which are clonally regulated.

Project Organization
------------

    │
    ├── data/    <--- https://drive.google.com/file/d/1pwDhRTRI-Oh_HisETlYkxtUThIepnVX5/view?usp=sharing
    │
    ├── figures/            <- Figures saved by scripts or notebooks.
    │
    ├── notebooks/          <- Jupyter notebooks. Naming convention is a short `-` delimited 
    │                         description, a number (for ordering), and the creator's initials,
    │                        e.g. `initial-data-exploration-01-hg`.
    │
    ├── output/             <- Manipulated data, logs, etc.
    │
    ├── tests/              <- Unit tests.
    │
    ├── seqclone/     <- Python module with source code of this project.
    │
    ├── environment.yml     <- conda virtual environment definition file.
    │
    ├── LICENSE
    │
    ├── Makefile            <- Makefile with commands like `make environment`
    │
    ├── README.md           <- The top-level README for developers using this project.
    │
    └── tox.ini             <- tox file with settings for running tox; see tox.testrun.org


--------


Set up
------------

Install the virtual environment with conda and activate it:

```bash
$ conda env create -f environment.yml
$ conda activate seqclone 
```

Install `seqclone` in the conda virtual environment:

```bash
$ pip install --editable .
```

Decompress data

```bash 
$ gunzip data/CombinedDivision.h5ad.gz
```
Or simply download it from the Google Drive link


Usage
------------
Path to the data is specified in the .ini file. Data can be downloaded from Google Drive.
```bash
$ conda activate seqclone
$ cd seqclone/ 
$ python CloneStats.py test.ini
``` 
