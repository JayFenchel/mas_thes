from MyMath import forward_substitution
from MyMath import backward_substitution
import numpy as np

# A=np.array([[1., 0., 0.],
#            [2., 1., 0.],
#            [1., 2., 1.]])
#
# b = np.array([[2., 3., 1.],
#      [2., 3., 1.],
#      [2., 4., 1.]])
#
# b = np.array([[2., 3.],
#      [2., 3.],
#      [2., 4.]])
#
# b = np.array([[2.],
#      [2.],
#      [2.]])

A = np.array([[  0.23996015,   0., 0.,   0., 0.],
              [ -0.37221757,   1., 0.,   0., 0.],
              [ -0.99008755,   0., 0.13885973,   0., 0.],
              [-48.93540655, 64.1, 2.39923411,   1., 0.],
              [  0.,           0., 0.,           0., 1.]])

b = np.array([[-1.2346445 ,1],
              [-1.43828223,1],
              [-4.48282454,1],
              [-1.79989043,1],
              [1.,1]])

x_f = forward_substitution(A, b)
if (abs(np.dot(A, x_f)-b)).sum() < 1e-10:
    print("forward_substitution did it right, error =", (abs(np.dot(A, x_f)-b)).sum())
else:
    print("forward_substitution failed")

x_b = backward_substitution(A.T, b)
if (abs(np.dot(A.T, x_b)-b)).sum() < 1e-10:
    print("backward_substitution did it right, error =", (abs(np.dot(A.T, x_b)-b)).sum())
else:
    print("backward_substitution failed")