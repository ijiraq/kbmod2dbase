import numpy
from mp_ephem import BKOrbit
import math
from .reader import DiscoveryFile, TrackingFile
from astropy import units
import argparse
import logging
import sys


def observation_triplet_from_kbmod_recrod(record):
    for dt in [0*units.hour, 1.5*units.hour, 3*units.hour]:
        yield record.offset(dt).observation


def load_orbits(file_in):
    orbits = {}
    for record in file_in:
        observations = []
        for observation in observation_triplet_from_kbmod_recrod(record):
            observations.append(observation)
        orbit = BKOrbit(observations)
        orbits[orbit.name] = orbit
    return orbits


def run():
    main_parser = argparse.ArgumentParser()
    main_parser.add_argument('discovery_file', type=str)
    main_parser.add_argument('tracking_file', type=str)
    main_parser.add_argument('field', type=str)
    main_parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'ERROR'], default='INFO')
    args = main_parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level))

    discovery_orbits = load_orbits(DiscoveryFile(args.field, args.discovery_file))
    tracking_orbits = load_orbits(TrackingFile(args.field, args.tracking_file))
    for discovery in discovery_orbits:
        for tracking in tracking_orbits:
            obs = discovery_orbits[discovery].observations.copy()
            obs.extend(tracking_orbits[tracking].observations)
            orbit = BKOrbit(obs)
            if not 20*units.au < orbit.distance < 300*units.au or orbit.a > 20*units.au:
                continue
            residuals = []
            for ob in obs:
                orbit.predict(ob.date)
                residuals.append(ob.coordinate.separation(orbit.coordinate).to('arcsec').value)
            if numpy.max(residuals) < 1:
                print(discovery, tracking)
                print(orbit.summarize())


if __name__ == '__main__':
    run()
