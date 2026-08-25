"""
Convert an input in either KBMOD discovery or tracking format to the TNOdb format.
The output is written to stdout
"""
from .reader import kbmod_file_iterator
from astropy import units
import argparse
import logging
import sys


def write_observations_as_tnodb_records(field, filename, nobs=3):
    """
    Write a triplet of ephemeris lines in tnodb format for each kbmod formatted observation in filename

    Each ephemeris line is offset by 0, 1.5, and 3 hours from the base observation time.

    :param field: The name of the survey field, used to determine the provisional designation
    :param filename: The name of the kbmod formatted input file
    :param nobs: The number of observations to generate per source
    """
    for record in kbmod_file_iterator(field, filename):
        dt = 1.5*units.hour
        with open(record.provisional_name+".tnodb", 'a') as f:
            [f.write(record.offset(i*dt).observation.to_tnodb() + "\n") for i in range(nobs)]


def main():
    main_parser = argparse.ArgumentParser()
    main_parser.add_argument('filename', type=str)
    main_parser.add_argument('field', type=str)
    main_parser.add_argument('--log-level',
                             choices=['DEBUG', 'INFO', 'ERROR'],
                             default='INFO')
    main_parser.add_argument('--num-of-obs', help="Number of observations to generate per source",
                             type=int, default=3)
    args = main_parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level))
    write_observations_as_tnodb_records(args.field, args.filename, nobs=args.num_of_obs)


if __name__ == '__main__':
    main()
