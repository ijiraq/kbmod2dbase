# kbmod2dbase
Reformat files coming from the kbmod outputs of the pipeline into the two-line format needed by dbase, the TNO astrometry database system developed for CFEPS and operated at UBC.


## Install

I recommend using a virtual environment, something like `venv`, to keep various Python packages cleanly separated.   Once you have activated your `venv` installation is via pip:

`pip install git+https://github.com/CLASSY-survey/kbmod2dbase.git`

## Usage

```
kbmod2tnodb --help

usage: kbmod2tnodb [-h] [--log-level {DEBUG,INFO,ERROR}] filename field

positional arguments:
  filename
  field

options:
  -h, --help            show this help message and exit
  --log-level {DEBUG,INFO,ERROR}
```

*filename* is the name of the input `kbmod` output file, and field is the name of the survey field whose `kbmod` file is being converted to *tnodb* format.  

The resulting *tnodb* lines are sent to standard out.

For each line of the `kbmod` file, the position, time, and rate information are used to produce an astrometric line in MPC 80 column format at three epochs.  The epoch of the first line is the time of the kbmod line and the two other times that are 1.5 and 3 hours later.  The measurement information (e.g. frame number of the exposure, likelihood value from `kbmod`, pixel locations, etc.) is stored in a comment line.  The format of the comment line is understood by the dbase ingester. 
