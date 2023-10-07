import math
from tempfile import NamedTemporaryFile
from unittest import TestCase

from astropy import units
from astropy.coordinates import SkyCoord
from astropy.units import Quantity

from kbmod2dbase.reader import TrackingFile, KBModRecord, DiscoveryFile


class TestTrackingFile(TestCase):

    def setUp(self) -> None:
        tracking_lines = (
            """35   1  59796.40457   640.00  3194.00   573.18  -207.56  23.51   311.24  334.531556  -12.492674   -4.377   -1.571
   8   4  59796.40457    92.00  2757.00  -565.84   226.80  23.31   288.57  334.460907  -11.802132   -4.346   -1.788
  13   3  59796.40457  2047.00  1465.00  -290.26   103.01  24.17   118.90  335.011083  -11.991661   -2.258   -0.813
   0  12  59796.40457   965.00  2787.00  -244.77   129.25  24.79    59.14  335.399034  -11.799053   -1.862   -0.991""")
        self.tracking_fobj = NamedTemporaryFile('w')
        self.tracking_fobj.writelines(tracking_lines)
        self.tracking_fobj.flush()

    def test_get_measure(self):
        with TrackingFile('test', self.tracking_fobj.name) as d:
            for record in d:
                print(record)
                print(type(record))
                self.assertAlmostEqual(record.x, 640.0 * units.pixel)
                self.assertAlmostEqual(record.ra, 334.531556 * units.degree)
                self.assertAlmostEqual(record.dec, -12.492674 * units.degree)
                self.assertAlmostEqual(record.mag, 23.51 * units.mag)
                self.assertAlmostEqual(
                    (record.ra_arc_rate - (-4.377 * units.arcsec / units.hour)).to('arcsec/hour').value,
                    0.0)
                break

    def tearDown(self) -> None:
        self.tracking_fobj.close()


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

