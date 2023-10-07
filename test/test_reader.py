import math
from tempfile import NamedTemporaryFile
from unittest import TestCase

from astropy import units
from astropy.coordinates import SkyCoord
from astropy.units import Quantity

from kbmod2dbase.reader import TrackingFile, KBModRecord, DiscoveryFile


class TestTrackingFile(TestCase):

    def setUp(self) -> None:
        self.tracking_line = ("35   1  59796.40457   640.00  3194.00   573.18  -207.56  23.51 311.24"
                              "  334.531556  -12.492674   -4.377   -1.571")
        self.tracking_file_handle = NamedTemporaryFile('w')
        self.tracking_file_handle.write(self.tracking_line + "\n")
        self.tracking_file_handle.flush()

    def test_get_measure(self):
        dec_rate = float(self.tracking_line.split()[12])
        ra_rate = float(self.tracking_line.split()[11])

        with TrackingFile('test', self.tracking_file_handle.name) as d:
            obs1 = next(d)
            for dt in [1.5*units.hour, 3*units.hour]:
                obs2 = obs1.offset(dt)
                dRA = (obs2.observation.coordinate.ra - obs1.observation.coordinate.ra).to('arcsec').value
                dDEC = (obs2.observation.coordinate.dec - obs1.observation.coordinate.dec).to('arcsec').value
                # print(math.cos(obs1.observation.coordinate.dec.to('rad').value)*dRA/dt.to('hour').value)
                self.assertAlmostEqual(math.cos(obs1.observation.coordinate.dec.to('rad').value)*dRA,
                                       ra_rate * dt.to('hour').value, 3)
                self.assertAlmostEqual(dDEC, dec_rate * dt.to('hour').value, 3)
            for record in d:
                self.assertAlmostEqual(record.x, 640.0 * units.pixel)
                self.assertAlmostEqual(record.ra, 334.531556 * units.degree)
                self.assertAlmostEqual(record.dec, -12.492674 * units.degree)
                self.assertAlmostEqual(record.mag, 23.51 * units.mag)
                self.assertAlmostEqual(record.observation.coordinate.dec, -12.492674 * units.degree)
                self.assertAlmostEqual(record.observation.coordinate.ra, 334.531556 * units.degree)
                self.assertAlmostEqual(
                    (record.ra_arc_rate - (-4.377 * units.arcsec / units.hour)).to('arcsec/hour').value,
                    ra_rate, 3)
                break

    def tearDown(self) -> None:
        self.tracking_file_handle.close()


class TestKBModRecord(TestCase):

    def setUp(self) -> None:
        self.kbmod_record = KBModRecord(survey_field='test', chip=Quantity(35, None),
                                        index=Quantity(1, None),
                                        mjd=Quantity(59796.40457, 'day'),
                                        x=Quantity(640.00, 'pixel'),
                                        y=Quantity(3194.00, 'pixel'),
                                        dx=Quantity(573.18, 'pixel/hour'),
                                        dy=-Quantity(207.56, 'pixel/hour'),
                                        mag=Quantity(23.51, 'mag'),
                                        merr=Quantity(0.99, 'mag'),
                                        likelihood=Quantity(311.24, None),
                                        ra=334.531556 * units.degree,
                                        dec=-12.492674 * units.degree,
                                        ra_arc_rate=-4.377 * units.arcsec / units.hour,
                                        dec_arc_rate=-1.571 * units.arcsec / units.hour)

    def test_offset(self):
        offset_record = self.kbmod_record.offset(1 * units.hour)
        dra = (offset_record.ra - self.kbmod_record.ra) * (math.cos(self.kbmod_record.dec.to('rad').value)) / units.hour
        ddec = (offset_record.dec - self.kbmod_record.dec) / units.hour
        self.assertAlmostEqual(self.kbmod_record.ra_arc_rate.to('arcsec/hour').value,
                               dra.to('arcsec/hour').value, 3)
        self.assertAlmostEqual(self.kbmod_record.dec_arc_rate.to('arcsec/hour').value,
                               ddec.to('arcsec/hour').value, 3)


class TestDiscoveryFile(TestCase):

    def setUp(self) -> None:
        self.discovery_lines = """   59 42.28  26.78 2 33 2
122.0 2929.0 361.48 -149.56 26.37 59817.336190 334.425922 -12.645317 -2.791127 -1.150713
83.0 2963.0 354.07 -166.35 27.18 59814.326828 334.481756 -12.622597 -2.732894 -1.279644"""

        self.discovery_file_obj = NamedTemporaryFile('w')
        self.discovery_file_obj.writelines(self.discovery_lines)
        self.discovery_file_obj.flush()

    def test_get_measure(self):
        with DiscoveryFile('test', self.discovery_file_obj.name) as d:
            for record in d:
                print(record)
                print(type(record))
                self.assertAlmostEqual(record.x, 122.0 * units.pixel)
                self.assertAlmostEqual(record.ra, 334.425922 * units.degree)
                self.assertAlmostEqual(record.dec, -12.645317 * units.degree)
                self.assertAlmostEqual(record.mag, 26.37 * units.mag)
                self.assertAlmostEqual(
                    (record.ra_arc_rate - (-2.791127 * units.arcsec / units.hour)).to('arcsec/hour').value,
                    0.0)
                break

    def tearDown(self) -> None:
        self.discovery_file_obj.close()

