import math
from tempfile import NamedTemporaryFile
from unittest import TestCase

from astropy import units
from astropy.coordinates import SkyCoord
from astropy.units import Quantity

from kbmod2dbase.reader import TrackingFile, KBModRecord, DiscoveryFile, GroupDiscoveryFile


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
        self.discovery_lines = """   4 236.80  25.89 0 2 2
          788.0 4319.0 -75.67 -33.80 26.39 59963.288713 120.360236 21.571612 -0.588165 0.259929 5.2
          1390.0 4394.0 -61.65 -28.30 25.39 59965.314781 120.352643 21.575596 -0.479325 0.217934 5.1
        """
        self.discovery_file_obj = NamedTemporaryFile('w')
        self.discovery_file_obj.writelines(self.discovery_lines)
        self.discovery_file_obj.flush()

    def test_get_measure(self):
        with DiscoveryFile('test', self.discovery_file_obj.name) as d:
            for record in d:
                print(record)
                print(type(record))
                self.assertAlmostEqual(record.x, 788.0 * units.pixel)
                self.assertAlmostEqual(record.ra, 120.360236 * units.degree)
                self.assertAlmostEqual(record.dec, 21.571612 * units.degree)
                self.assertAlmostEqual(record.mag, 26.39 * units.mag)
                self.assertAlmostEqual(
                    (record.ra_arc_rate - (-0.588165 * units.arcsec / units.hour)).to('arcsec/hour').value,
                    0.0)
                break

    def tearDown(self) -> None:
        self.discovery_file_obj.close()

class TestGroupFile(TestCase):

    def setUp(self) -> None:
        self.discovery_lines = """3 N123 WARNING: MULTI GROUP! G 
        0 963.0 2257.0 335.115552 -11.883817 59813.360306 -1.18 -0.21 57 35.66 -1
        0 1174.0 2170.0 335.107898 -11.88629 59814.330278 -1.1 -0.43 65 123.56 72979860
        0 1762.0 1925.0 335.084753 -11.894749 59817.33964 -1.21 -0.47 85 144.38 72979860
        """
        self.discovery_file_obj = NamedTemporaryFile('w')
        self.discovery_file_obj.writelines(self.discovery_lines)
        self.discovery_file_obj.flush()

    def test_get_measure(self):
        with GroupDiscoveryFile('test', self.discovery_file_obj.name) as d:
            for record in d:
                self.assertAlmostEqual(record.x, 963.0 * units.pixel)
                self.assertAlmostEqual(record.ra, 335.115552 * units.degree)
                self.assertAlmostEqual(record.dec, -11.883817 * units.degree)
                self.assertAlmostEqual(record.mjd, 59813.360306 * units.day)
                self.assertAlmostEqual(record.likelihood, 35.66)
                self.assertEqual(record.flag, 'G')
                self.assertEqual(record.objid, -1)
                self.assertAlmostEqual(
                    (record.ra_arc_rate - (-1.18 * units.arcsec / units.hour)).to('arcsec/hour').value,
                    0.0)
                break

    def tearDown(self) -> None:
        self.discovery_file_obj.close()


class TestGroupRejectedFile(TestCase):
    """Rejected group files omit G/B flags and the trailing objid column."""

    def setUp(self) -> None:
        self.discovery_lines = """58 N13
	0 274.0 1736.0 335.079198 -11.857041 59813.360306 -3.05 -0.94 10 5.90
	0 150.0 1763.0 335.000477 -11.886194 59817.33964 -2.89 -1.13 7 9.00
98 N23 WARNING: MULTI GROUP!
	8 519.0 2074.0 334.181652 -11.884561 59814.330278 -2.65 -0.9 34 6.99
	8 454.0 1992.0 334.124488 -11.901284 59817.33964 -2.73 -0.89 31 7.95
"""
        self.discovery_file_obj = NamedTemporaryFile('w')
        self.discovery_file_obj.writelines(self.discovery_lines)
        self.discovery_file_obj.flush()

    def test_get_measure(self):
        with GroupDiscoveryFile('AS1', self.discovery_file_obj.name) as d:
            records = list(d)
        self.assertEqual(len(records), 4)
        first = records[0]
        self.assertAlmostEqual(first.x, 274.0 * units.pixel)
        self.assertAlmostEqual(first.ra, 335.079198 * units.degree)
        self.assertAlmostEqual(first.dec, -11.857041 * units.degree)
        self.assertAlmostEqual(first.mjd, 59813.360306 * units.day)
        self.assertAlmostEqual(first.likelihood, 5.90)
        self.assertEqual(first.flag, '_')
        self.assertEqual(first.objid, 0)
        self.assertEqual(first.detkey, 'N13')
        last = records[-1]
        self.assertEqual(last.index, 98)
        self.assertEqual(last.detkey, 'N23')
        self.assertIn('MULTI', last.comment)

    def tearDown(self) -> None:
        self.discovery_file_obj.close()

