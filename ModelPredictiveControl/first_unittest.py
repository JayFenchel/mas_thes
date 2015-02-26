__author__ = 'jayf'

from QuadraticProgram import QuadraticProgram
import unittest



class MyTestCase(unittest.TestCase):

    def setUp(self):
        # sys = AirCraft()  # later nicht mehr nötig
        self.test_qp = QuadraticProgram()

    # Quadratic Program
    def test_set_constraints(self):

        self.test_qp.set

        self.assertEqual(True, False)


if __name__ == '__main__':
    unittest.main()
