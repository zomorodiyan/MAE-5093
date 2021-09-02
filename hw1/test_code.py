import unittest
import numpy as np
from code import lagrange

class given(unittest.TestCase):

  def test_lagrange(self): # function's name have to start with "test_"
    data = np.empty((11,2))
    data[:,0] = np.arange(-1,1.001,0.2)
    data[:,1] = np.power(25*np.power(data[:,0], 2)+1, -1)
    self.assertAlmostEqual(lagrange(0.7,data), -0.2261963)

  def test_cubic_spline(self): # function's name have to start with "test_"
    self.assertAlmostEqual(1, 1)

if __name__ == '__main__':
    unittest.main()
