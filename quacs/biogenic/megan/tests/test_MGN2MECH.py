import unittest
import warnings
import numpy as np
from quacs.biogenic.megan.src.MGN2MECH import MGN2MECH
class TestMapping(unittest.TestCase):
    def setUp(self):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            self.m=MGN2MECH()
    def test_dimensions(self):
        self.assertEqual(self.m.to_species(np.ones(19)).shape,(201,))
        self.assertEqual(self.m.convert(np.ones((19,2,3))).shape,(len(self.m.mechanism_names),2,3))
    def test_isoprene(self):
        x=np.zeros(19);x[0]=7
        s=self.m.to_species(x)
        self.assertAlmostEqual(s[0],7)
        y=self.m.convert(x)
        self.assertAlmostEqual(y[self.m.mechanism_names.index('ISOP')],7)
    def test_acetone(self):
        x=np.zeros(19);x[self.m.class_names.index('ACTO')]=3
        y=self.m.convert(x)
        self.assertAlmostEqual(y[self.m.mechanism_names.index('KET')],3)
        self.assertAlmostEqual(y[self.m.mechanism_names.index('PAR')],6)
    def test_wrong_shape(self):
        with self.assertRaises(ValueError):self.m.convert(np.ones(18))
    def test_diagnostics(self):
        self.assertIn(48,self.m.duplicate_species_ids)
        self.assertTrue(self.m.duplicate_mechanism_rows)
if __name__=='__main__':unittest.main()
