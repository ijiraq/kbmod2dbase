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

And the "Grouping" format is currently:
Index, Detection key, optional comment, [G/B] (Good/Bad)
The first line is a header line, and the following lines are the measurements of the object.
        ID X Y RA DEC MJD RA_arc_rate DEC_arc_rate Number likelihood ObjID
#
"""
import string
import inspect
from abc import ABC, abstractmethod
from collections import OrderedDict
from copy import deepcopy

from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.units import Quantity
from mp_ephem.ephem import Observation
from dataclasses import dataclass, field


# KEYS = string.digits + string.ascii_uppercase + string.ascii_lowercase
KEYS = string.digits + string.ascii_uppercase

class MissingColumnsError(Exception):
    pass


def year_to_letter(year) -> chr:
    """Convert a year value to a letter"""
    return chr(year - 2000 + ord('A'))


def index_to_key(index: int):
    """Return a two character string for the index"""
    if not index:
        return "0"
    b = len(KEYS)
    return index_to_key(index//b).lstrip("0") + KEYS[index % b]


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
    provisional_name: str # name of the object
    chip: int = -1 # CCD number
    index: int  # "Index of the object in the field"
    mjd: Quantity  # "MJD of the observation"
    x: Quantity  # "X position of the object on the CCD"
    y: Quantity  # "Y position of the object on the CCD"
    dx: Quantity = None # "X velocity of the object on the CCD"
    dy: Quantity = None # "Y velocity of the object on the CCD"
    ra: Quantity  # "Right ascension of the object"
    dec: Quantity  # "Declination of the object"
    ra_arc_rate: Quantity  # "rate of RA motion"
    dec_arc_rate: Quantity  # "rate of DEC motion"
    mag: Quantity = Quantity(99.9, 'mag')  # "Magnitude of the object"
    merr: Quantity = Quantity(9.99, 'mag')  # "Magnitude error of the object"
    likelihood: int = -1  # "Likelihood of the object"
    flag: str = ""  # "Flag(G/B) on this object"
    detkey: str = ""  # "Detection key for the object"
    objid: int = 0  # "Index of the object in the field"
    comment: str = field(default=None)
    coord: SkyCoord = field(init=False)
    date: Time = field(init=False)
    frame: str = field(init=False)
    band: str = field(init=False, default='r')
    observatory_code: str = field(init=False, default='568')

    def __post_init__(self) -> None:
        self.date = Time(self.mjd, format='mjd', precision=5)
        self.coord = SkyCoord(ra=self.ra, dec=self.dec, unit='deg', frame='icrs', obstime=self.date)
        self.frame = f"{self.survey_field}{self.date.strftime('%y%m%d')}{int(self.chip):02d}"
        scale = Quantity(0.185, 'arcsec/pixel')
        self.dx = self.dx if self.dx is not None else self.ra_arc_rate/scale
        self.dy = self.dy if self.dy is not None else self.dec_arc_rate/scale
        comment = self.comment if self.comment else ""
        self.comment = f"{self.flag} {self.detkey} {self.objid} {comment.replace("_"," ")}"

    def offset(self, dt: Quantity) -> 'KBModRecord':
        """
        Return a new KBModRecord with the position offset by the velocity times dt
        """
        new_coord = self.coord.spherical_offsets_by(self.ra_arc_rate * dt,
                                                    self.dec_arc_rate * dt)
        record = deepcopy(self)
        record.ra = new_coord.ra
        record.dec = new_coord.dec
        record.x = record.x + self.dx*dt
        record.y = record.y + self.dy*dt
        record.date = record.date + dt
        return record

    def __str__(self):
        return f"{self.provisional_name} {self.date} " \
               f"{self.ra} {self.dec} {self.ra_arc_rate} {self.dec_arc_rate} " \
               f"{self.mag} {self.merr} {self.likelihood}"

    @property
    def observation(self) -> Observation:
        return Observation(provisional_name=self.provisional_name,
                           null_observation=self.flag == "B",
                           frame=self.frame,
                           survey_code='C',
                           mag=self.mag.to('mag').value,
                           mag_err=self.merr.to('mag').value,
                           xpos=self.x.to('pixel').value,
                           ypos=self.y.to('pixel').value,
                           ra=self.ra,
                           dec=self.dec,
                           date=self.date.mpc,
                           band=self.band,
                           observatory_code=self.observatory_code,
                           comment=self.comment,
                           likelihood=self.likelihood)


class KBModFileIterator(ABC):
    """
    Open an iterator over a KBMOD file, there are two types of KBMOD files, Discovery and Tracking
    """
    FIRST_ROW = OrderedDict()
    OBSERVATION_COLUMNS = OrderedDict()
    OBJECT_COLUMNS = OrderedDict()

    def __init__(self, survey_field, filename):
        self._line = None
        self._object_info = {}
        self.filename = filename
        self.survey_field = survey_field
        self._file_object = None

    @abstractmethod
    def get_provisional_name(self, measure) -> str:
        raise NotImplementedError

    @property
    def observation_row_columns(self) -> OrderedDict:
        return self.OBSERVATION_COLUMNS

    @property
    def object_row_columns(self) -> OrderedDict:
        return self.OBJECT_COLUMNS

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
            raise MissingColumnsError(f"Expected {len(column_definitions)} columns, got {len(raw_values)}\n"
                                      f"{column_definitions}\n"
                                      f"{self.line}")
        quantities = [self.convert_into_quantity_if_possible(raw_value, unit) for (raw_value, unit)
                      in zip(self.caste_based_on_decimal_point(raw_values),
                             column_definitions.values())]
        return {column_name: quantity for (column_name, quantity) in zip(column_definitions.keys(), quantities)}

    def caste_based_on_decimal_point(self, raw_values: list) -> list:
        values = []
        for row in raw_values:
            func = float if '.' in row else int
            try:
                values.append(func(row))
            except ValueError:
                values.append(row)
        return values

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
        measure['provisional_name'] = self.get_provisional_name(measure)
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
       0  15  59815.34616  1120.00  4131.00  -367.97   132.79  26.74     8.93  335.070908  -12.000795   -2.807   -1.012
    """

    OBSERVATION_COLUMNS = OrderedDict((('chip', None),
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

    OBJECT_COLUMNS = OrderedDict()
    FIRST_ROW = OBSERVATION_COLUMNS

    def get_provisional_name(self, measure) -> str:
        date = Time(measure['mjd'], format='mjd', precision=5)
        year_key = year_to_letter(date.datetime.year)
        day = date.datetime.timetuple().tm_yday
        ccd_key = chip_to_key(int(measure['chip']))
        detection_index_key = index_to_key(int(measure['index'])).zfill(5)
        return f"{self.survey_field}{year_key}{day:03d}{ccd_key}{detection_index_key}"


class DiscoveryFile(KBModFileIterator):
    """
    Read in a Detection file.  Creates an iterator that returns sets of KBModRecords for the Discovery file
    """
    OBJECT_COLUMNS = OrderedDict((('index', None),
                               ('dist', 'au'),
                               ('mag', 'mag'),
                               ('visit', None),
                               ('chip', None),
                               ('ndet', None)))

    OBSERVATION_COLUMNS = OrderedDict((('x', 'pixel'),
                                       ('y', 'pixel'),
                                       ('dx', 'pixel/day'),
                                       ('dy', 'pixel/day'),
                                       ('mag', 'mag'),
                                       ('mjd', 'day'),
                                       ('ra', 'degree'),
                                       ('dec', 'degree'),
                                       ('ra_arc_rate', 'arcsec/hour'),
                                       ('dec_arc_rate', 'arcsec/hour'),
                                       ('likelihood', None)))

    FIRST_ROW = OBJECT_COLUMNS

    def get_provisional_name(self, measure) -> str:
        detection_index_key = index_to_key(measure['index']).zfill(3)
        return f"{self.survey_field}{measure['index']:04d}"


class GroupDiscoveryFile(DiscoveryFile):
    """
    Read in a Detection file.  Creates an iterator that returns sets of KBModRecords for the Discovery file
    """
    OBJECT_COLUMNS = OrderedDict((('index', None),
                                  ('detkey', None),
                                  ('flag', None),
                                  ('comment', None)))

    # Observation Line (objid is optional; rejected sources often omit it):
    #         0 384.0 4007.0 335.067704 -11.980512 59814.330278 -0.62 -0.37 19 23.64 21785381
    #         0 274.0 1736.0 335.079198 -11.857041 59813.360306 -3.05 -0.94 10 5.90
    OBSERVATION_COLUMNS = OrderedDict((
        ('counter', None),
        ('x', 'pixel'),
        ('y', 'pixel'),
        ('ra', 'degree'),
        ('dec', 'degree'),
        ('mjd', 'day'),
        ('ra_arc_rate', 'arcsec/hour'),
        ('dec_arc_rate', 'arcsec/hour'),
        ('number', None),
        ('likelihood', None),
        ('objid', None)))

    FIRST_ROW = OBJECT_COLUMNS

    def _is_group_object_header(self, parts: list) -> bool:
        """True for headers like '58 N13', '101 N123 night 1 & 3 ... G'."""
        return (len(parts) >= 2
                and parts[0].lstrip('-').isdigit()
                and parts[1].startswith('N')
                and parts[1][1:].isdigit())

    def is_observation_line(self) -> bool:
        parts = self.line.strip().split()
        if self._is_group_object_header(parts):
            return False
        n_cols = len(parts)
        # With or without trailing objid
        return n_cols in (len(self.observation_row_columns) - 1, len(self.observation_row_columns))

    def is_object_line(self) -> bool:
        """
        Object headers are index + night key (N13/N123/...), optionally a comment and G/B.
        Observation rows are purely numeric with 10 or 11 columns.
        """
        return self._is_group_object_header(self.line.strip().split())

    def parse_line(self, column_definitions) -> dict:
        """
        Override the parse_line method to handle the group discovery file format
        to first check if line is an object line or an observation line
        and if an object line then trim it down to the object columns

        Possible Object Line Forms:
        58 N13
        98 N23 WARNING: MULTI GROUP!
        2646 N123 night 1 & 2 sources dont match WARNING: MULTI GROUP!
        2633 N13 WARNING: MULTI GROUP! G
        2624 N13 G
        2162 N123 G

        First value is the index
        Second is a key indicating which nights the source was detected on.
        Optional middle tokens are a comment of variable length
        Optional last token is G or B indicating Good/Bad after inspection
        """
        if self.is_object_line():
            raw_values = self.line.strip().split()
            if raw_values[-1] in ('G', 'B'):
                flag = raw_values[-1]
                comment_parts = raw_values[2:-1]
            else:
                flag = '_'
                comment_parts = raw_values[2:]
            comment = '_'.join(comment_parts) if comment_parts else '_'
            self._line = f"{raw_values[0]} {raw_values[1]} {flag} {comment}"
        elif self.is_observation_line():
            raw_values = self.line.strip().split()
            if len(raw_values) == len(self.observation_row_columns) - 1:
                # Pad missing objid
                self._line = f"{self.line.strip()} 0"
        return super().parse_line(column_definitions)

def _get_first_line(filename) -> list:
    with open(filename, 'r') as fobj:
        while True:
            first_line = fobj.readline().strip()
            if first_line.startswith("#"):
                continue
            break
    return first_line.split()


def _looks_like_group_file(filename: str) -> bool:
    """Group headers look like: '<index> N13 ...' or '<index> N123 ... G'."""
    first = _get_first_line(filename)
    return (len(first) >= 2
            and first[0].lstrip('-').isdigit()
            and first[1].startswith('N')
            and first[1][1:].isdigit())

def kbmod_file_iterator(survey_field: str, filename: str) -> KBModFileIterator:
    """
    Parse a KBMOD file and return a list of KBModRecord objects
    """
    if _looks_like_group_file(filename):
        return GroupDiscoveryFile(survey_field, filename)
    for cls in KBModFileIterator.__subclasses__():
        if cls is GroupDiscoveryFile:
            continue
        if len(cls.FIRST_ROW) == len(_get_first_line(filename)):
            try:
                return cls(survey_field, filename)
            except TypeError:
                pass
    raise ValueError(f"Unknown file format for {filename}")

