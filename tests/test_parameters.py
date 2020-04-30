import unittest

from ..parameters import Element, Configuration, Cowan  # noqa


class TestElement(unittest.TestCase):
    def setUp(self):
        self.element = Element("Fe", "3+")

    def test_valence_occupancy(self):
        self.assertEqual(self.element.valence_occupancy, 5)

    def test_valence_subshell(self):
        self.assertEqual(self.element.valence_subshell, "3d")


class TestConfiguration(unittest.TestCase):
    def setUp(self):
        self.conf = Configuration("1s1,3d6")
        self.element = Element("Fe")
        self.cowan = Cowan(self.element, self.conf)
        self.conf.energy, *self.conf.atomic_parameters = self.cowan.get_parameters()

    def test_subshells(self):
        self.assertCountEqual(self.conf.subshells, ["1s", "3d"])

    def test_energy(self):
        self.assertAlmostEqual(self.conf.energy, -27438.6132941, places=6)

    def test_parameters(self):
        f2 = self.conf.atomic_parameters["F2(3d,3d)"]
        f4 = self.conf.atomic_parameters["F4(3d,3d)"]
        g2 = self.conf.atomic_parameters["G2(1s,3d)"]
        zeta = self.conf.atomic_parameters["ζ(3d)"]

        self.assertEqual(f2, 12.7364)
        self.assertEqual(f4, 7.9634)
        self.assertEqual(g2, 0.0655)
        self.assertEqual(zeta, 0.0753)


def suite():
    loadTests = unittest.defaultTestLoader.loadTestsFromTestCase
    test_suite = unittest.TestSuite()
    test_suite.addTest(loadTests(TestElement))
    test_suite.addTest(loadTests(TestConfiguration))
    return test_suite


if __name__ == "__main__":
    unittest.main(defaultTest="suite")
