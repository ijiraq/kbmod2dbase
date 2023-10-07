from .reader import DiscoveryFile, TrackingFile
from astropy import units
import argparse
import logging
import sys


def write_triplet_of_observation(record, file_out):
    for dt in [0*units.hour, 1.5*units.hour, 3*units.hour]:
        file_out.write(str(record.offset(dt).observation.to_tnodb()) + "\n")


def write_records(file_in, file_out):
    for record in file_in:
        write_triplet_of_observation(record, file_out)


def run():
    main_parser = argparse.ArgumentParser()
    main_parser.add_argument('filename', type=str)
    main_parser.add_argument('field', type=str)
    main_parser.add_argument('output', type=str, nargs='?', default=None)
    main_parser.add_argument('--tracking', action="store_true", default=False)
    main_parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'ERROR'], default='INFO')
    args = main_parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level))

    cls = args.tracking and TrackingFile or DiscoveryFile
    file_out = args.output is None and sys.stdout or open(args.output, 'w')
    with cls(args.field, args.filename) as file_in:
        write_records(file_in, file_out)


if __name__ == '__main__':
    run()
