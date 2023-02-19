import datetime

from numpy.typing import NDArray
import numpy as np



## Precession Polynomials (Arc-Seconds)
# Equation 5.16
precession_x = np.polynomial.polynomial.Polynomial(
    [-0.016617, 2004.191898, -0.4297829, -0.19861834]
)
precession_y = np.polynomial.polynomial.Polynomial(
    [-0.006951,  -0.025896, -22.4072747, 0.0019005900]
)

# Nutation Polynomials (Arc-Seconds)
# Equation 5.43

# Mean_anomaly of the moon
moon_anomoly = np.polynomial.polynomial.Polynomial(
    []
)