"""
Various classes to read the output of the KBMOD astrometry code and create Observation objects

Create an Observation from the 'rough' kbmod astrometry file

The current format of the discovery files is series of entries like the following:

detection #, distance estimate (au), magnitude, visit # (0/1/2), chip, # of detections
        x1 y1 dx1 dy1 mjd1 ra1 dec1 cos(dec)*dra1 ddec1
        X2 y2 dx2 dy2 mjd2 ra2 dec2 cos(dec)*dra2 ddec2
        x3 y3 dx3 dy3 mjd3 ra3 dec3 cos(dec)*dra3 ddec3

The first line is a header line, and the following lines are the measurements of the object.

The current format of the tracking files is series of entries like the following:

chip index mjd x y dx/dt dy/dt mag likelihood ra dec cos(dec)*ra_arc_rate/dt dec_arc_rate/dt

#
"""
import string
import re
import inspect
from collections import OrderedDict

import orderdict
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.units import Quantity
from mp_ephem.ephem import Observation
from dataclasses import dataclass, field
from abc import ABC, abstractmethod

KEYS = string.ascii_uppercase + string.ascii_lowercase + string.digits


class MissingColumnsError(Exception):
    pass


def year_to_letter(year) -> chr:
    """Convert a year value to a letter"""
    return chr(year - 2000 + ord('A'))


def index_to_key(index: int):
    """Return a two character string for the index"""
    return KEYS[index // 62] + KEYS[index % 62]


def key_to_index(key):
    """Return the index for a two character string"""
    return KEYS.index(key[0]) * 62 + KEYS.index(key[1])


def chip_to_key(chip: int):
    """Return a single character string for the chip number"""
    return KEYS[chip]


def key_to_chip(key):
    """Return the chip number for a single character string"""
    return KEYS.index(key)


@dataclass(frozen=False, kw_only=True, slots=True)
class KBModRecord:
    """
    A class to hold info about an object measured by
    the classy pipeline using KBMod
    """
    survey_field: str  # "Name of the CLASSY field"
    chip: int  # CCD number
    index: int  # "Index of the object in the field"
    mjd: Quantity  # "MJD of the observation"
    x: Quantity  # "X position of the object on the CCD"
    y: Quantity  # "Y position of the object on the CCD"
    dx: Quantity  # "X velocity of the object on the CCD"
    dy: Quantity  # "Y velocity of the object on the CCD"
    ra: Quantity  # "Right ascension of the object"
    dec: Quantity  # "Declination of the object"
    ra_arc_rate: Quantity  # "rate of RA motion"
    dec_arc_rate: Quantity  # "rate of DEC motion"
    mag: Quantity = Quantity(99.99, 'mag')  # "Magnitude of the object"
    merr: Quantity = Quantity(0.99, 'mag')  # "Magnitude error of the object"
    likelihood: int = -1  # "Likelihood of the object"
    coord: SkyCoord = field(init=False)
    date: Time = field(init=False)
    frame: str = field(init=False)
    observation: Observation = field(init=False)

    def __post_init__(self) -> None:
        self.date = Time(self.mjd, format='mjd', precision=6)
        self.coord = SkyCoord(ra=self.ra, dec=self.dec, unit='deg', frame='icrs', obstime=self.date)
        self.frame = f"{self.survey_field}{self.date.strftime('%y%m%d')}{int(self.chip):02d}"
        self.observation = Observation(provisional_name=self.provisional_name,
                                       frame=self.frame,
                                       survey_code='C',
                                       mag=self.mag.to('mag').value,
                                       mag_err=self.merr.to('mag').value,
                                       xpos=self.x.to('pixel').value,
                                       ypos=self.y.to('pixel').value,
                                       ra=self.ra,
                                       dec=self.dec,
                                       date=self.date.mpc,
                                       band='r',
                                       observatory_code=568,
                                       comment="",
                                       likelihood=self.likelihood)

    @property
    def provisional_name(self) -> str:
        year_key = year_to_letter(self.date.datetime.year)
        day = self.date.datetime.timetuple().tm_yday
        ccd_key = chip_to_key(int(self.chip))
        detection_index_key = index_to_key(int(self.index))
        return f"{self.survey_field}{year_key}{day:03d}{ccd_key}{detection_index_key}"

    def offset(self, dt: Quantity) -> 'KBModRecord':
        """
        Return a new KBModRecord with the position offset by the velocity times dt
        """
        new_coord = self.coord.spherical_offsets_by(self.ra_arc_rate * dt,
                                                    self.dec_arc_rate * dt)
        return KBModRecord(survey_field=self.survey_field,
                           chip=self.chip,
                           index=self.index,
                           mjd=self.mjd + dt,
                           x=self.x + self.dx * dt,
                           y=self.y + self.dy * dt,
                           dx=self.dx,
                           dy=self.dy,
                           mag=self.mag,
                           merr=self.merr,
                           likelihood=self.likelihood,
                           ra=new_coord.ra,
                           dec=new_coord.dec,
                           ra_arc_rate=self.ra_arc_rate,
                           dec_arc_rate=self.dec_arc_rate)

    def __str__(self):
        return f"{self.provisional_name} {self.date} " \
               f"{self.ra} {self.dec} {self.ra_arc_rate} {self.dec_arc_rate} " \
               f"{self.mag} {self.merr} {self.likelihood}"


class KBModFileIterator(ABC):
    """
    Open an iterator over a KBMOD file, there are two types of KBMOD files, Discovery and Tracking
    """

    def __init__(self, survey_field, filename):
        self._line = None
        self._object_info = {}
        self.filename = filename
        self.survey_field = survey_field
        self._file_object = None

    @property
    @abstractmethod
    def observation_row_columns(self) -> OrderedDict:
        pass

    @property
    @abstractmethod
    def object_row_columns(self) -> OrderedDict:
        pass

    @property
    def line(self):
        if self._line is None:
            self.next_line()
        return self._line

    def is_observation_line(self) -> bool:
        return len(self.line.strip().split()) == len(self.observation_row_columns)

    def is_object_line(self) -> bool:
        return len(self.line.strip().split()) == len(self.object_row_columns)

    def __enter__(self):
        self._file_object = open(self.filename, 'r')
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._file_object.close()
        self._file_object = None

    @property
    def file_object(self):
        if self._file_object is None:
            self._file_object = open(self.filename, 'r')
        return self._file_object

    def parse_line(self, column_definitions) -> dict:
        """
        Parse a line from a KBMOD file
        """
        raw_values = self.line.strip().split()
        line_contains_expected_number_of_columns = len(column_definitions) == len(raw_values)
        if not line_contains_expected_number_of_columns:
            # Wrong number of records for a discovery observation
            raise MissingColumnsError(f"Expected {len(column_definitions)} columns, got {len(raw_values)}\n"
                                      f"{column_definitions}\n"
                                      f"{self.line}")
        quantities = [self.convert_into_quantity_if_possible(raw_value, unit) for (raw_value, unit)
                      in zip(self.caste_based_on_decimal_point(raw_values),
                             column_definitions.values())]
        return {column_name: quantity for (column_name, quantity) in zip(column_definitions.keys(), quantities)}

    def caste_based_on_decimal_point(self, raw_values: list) -> list:
        return [func(raw_value) for (raw_value, func) in
                zip(raw_values, ['.' in raw_value and float or int for raw_value in raw_values])]

    def convert_into_quantity_if_possible(self, value, unit) -> (Quantity, int, float):
        if unit is None:
            return value
        return Quantity(value, unit)

    def get_measure(self) -> dict:
        """Return the next measurement of an object from the kbmod input file"""
        self.next_line()
        if self.is_object_line():
            self._object_info = self.parse_line(self.object_row_columns)
            self.next_line()
        measure = self._object_info.copy()
        measure.update(self.parse_line(self.observation_row_columns))
        return measure

    def next_line(self) -> None:
        self._line = self.file_object.readline()
        if self._line == "":
            raise StopIteration
        if self._line.strip().startswith("#"):
            self.next_line()

    def __iter__(self):
        return self

    def __next__(self):
        """Return the next 'object' from the kbmod input file, or raise StopIteration"""
        try:
            measure = self.get_measure()
        except MissingColumnsError:
            raise StopIteration
        if measure is None:
            raise StopIteration
        # gather all the observations of this source into a KBModRecord
        init_arguments = inspect.signature(KBModRecord).parameters.keys()
        kwargs = dict((argument, measure[argument]) for argument in init_arguments if argument in measure)
        return KBModRecord(survey_field=self.survey_field, **kwargs)


class TrackingFile(KBModFileIterator):
    """
    Class to loop over tacking observation file from classy.
    """

    @property
    def observation_row_columns(self) -> OrderedDict:
        return OrderedDict((('chip', None),
                            ('index', None),
                            ('mjd', 'day'),
                            ('x', 'pixel'),
                            ('y', 'pixel'),
                            ('dx', 'pixel/day'),
                            ('dy', 'pixel/day'),
                            ('mag', 'mag'),
                            ('likelihood', None),
                            ('ra', 'degree'),
                            ('dec', 'degree'),
                            ('ra_arc_rate', 'arcsec/hour'),
                            ('dec_arc_rate', 'arcsec/hour')))

    @property
    def object_row_columns(self) -> OrderedDict:
        return OrderedDict()


class DiscoveryFile(KBModFileIterator):
    """
    Read in a Detection file.  Creates an iterator that returns sets of KBModRecords for the Discovery file
    """
    OBJ_START_PATTERN = (r'\s*(?P<id>\d+)'
                         r'\s+(?P<dist>\d+(\.\d*)?)'
                         r'\s+(?P<mag>\d+(\.\d*)?)'
                         r'\s+(?P<visit>\d+)'
                         r'\s+(?P<chip>\d+)'
                         r'\s+(?P<ndet>\d+)\s*')
    OBJ_START_LINE = re.compile(OBJ_START_PATTERN)
    OBJ_START_COLUMNS = "index dist mag visit chip ndet".split()

    @property
    def object_row_columns(self) -> OrderedDict:
        return OrderedDict((('index', None),
                            ('dist', 'au'),
                            ('mag', 'mag'),
                            ('visit', None),
                            ('chip', None),
                            ('ndet', None)))

    @property
    def observation_row_columns(self) -> OrderedDict:
        return OrderedDict((('x', 'pixel'),
                            ('y', 'pixel'),
                            ('dx', 'pixel/hour'),
                            ('dy', 'pixel/hour'),
                            ('mag', 'mag'),
                            ('mjd', 'day'),
                            ('ra', 'degree'),
                            ('dec', 'degree'),
                            ('ra_arc_rate', 'arcsec/hour'),
                            ('dec_arc_rate', 'arcsec/hour')))
