"""
Convert an input in either KBMOD discovery or tracking format to the TNOdb format.
The output is written to stdout
"""
import kbmod2dbase
from astropy import units
import argparse
import logging
import sys


def write_observations_as_tnodb_records(field, filename):
    """
    Write a triplet of ephemeris lines in tnodb format for each kbmod formatted observation in filename

    Each ephemeris line is offset by 0, 1.5, and 3 hours from the base observation time.

    :param field: The name of the survey field, used to determine the provisional designation
    :param filename: The name of the kbmod formatted input file
    """
    for record in kbmod2dbase.kbmod_file_iterator(field, filename):
        [sys.stdout.write(record.offset(dt).observation.to_tnodb() + "\n") for dt in
         [0*units.hour, 1.5*units.hour, 3*units.hour]]


def main():
    main_parser = argparse.ArgumentParser()
    main_parser.add_argument('filename', type=str)
    main_parser.add_argument('field', type=str)
    main_parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'ERROR'], default='INFO')
    args = main_parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level))
    write_observations_as_tnodb_records(args.field, args.filename)


if __name__ == '__main__':
    main()
