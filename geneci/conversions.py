from datetime import datetime

from numpy.typing import NDArray
import numpy as np

ARC_SECONDS_TO_RADIANS = np.pi / 180 / 3600
DERIVATIVE_CONSTANT = 1.00273781191135448 / 3600 / 24  


def rotation_matrix_ecef_to_eci(
    julian_century: float
) -> NDArray[np.float64]:
    """Return ecef to eci rotation matrix for a given time.

    Calculation is

    P_eci = R(t) @ P_ecef


    Args:
        julian_date (float): Time in Julian centuries.


    Returns:
        rotation_matrx (NDArray[float]): 3x3 Matrix to rotate ECEF to ECI
    """
    # Angular motion of the earth
    earth_rotation_angle = np.pi * (
        0.7790572732640 + 1.00273781191135448 * 36525.0 * julian_century
    )
    earth_matrix = np.eye(3)
    earth_matrix[0, 0] = earth_matrix[1, 1] = np.cos(earth_rotation_angle)
    earth_matrix[1, 0] = np.sin(earth_matrix)
    earth_matrix[0, 1] = -1 * earth_matrix[1, 0]

    # Precession / Nutation rotation matrix (Eq. 5.10)
    gcrs_x, gcrs_y = compute_celestial_positions(julian_century)
    a = 0.5 + 0.125 * (gcrs_x*gcrs_x + gcrs_y*gcrs_y)
    
    pn_matrix = np.array([
        [1-a*gcrs_x*gcrs_x,  -a*gcrs_x*gcrs_y, gcrs_x],
        [ -a*gcrs_x*gcrs_y, 1-gcrs_y*gcrs_y  , gcrs_y],
        [ -gcrs_x         ,  -gcrs_y         , 1-a*(gcrs_x*gcrs_x + gcrs_y*gcrs_y)]
    ])

    # Return the rotation
    return pn_matrix @ earth_matrix


## The monstrosity that is the nutation / precession
def compute_celestial_positions(
    julian_century: float
)  -> tuple[float, float]:
    """Compute the x-y components of the celestial pole in earth reference frame.

    Args:
        julian_date (float): Time in Julian centuries.

    Returns:
        celestial_x (float): x-component of the pole vector
        celestial_y (float): y-component of the pole vector

    Notes:
        See Equation 5.16 with Table 5.2a / 5.2b. Supplemental material has
        all 2000 parameters or so. Download zip file from website
    """
    celestial_x = precession_x(julian_century)
    celestial_y = precession_y(julian_century)

    # Update for nutation (Coeffients are micro-arc-seconds)
    omega = moon_ascension(julian_century) * ARC_SECONDS_TO_RADIANS
    D =  moon_elongation(julian_century) * ARC_SECONDS_TO_RADIANS
    F = moon_longitude(julian_century) * ARC_SECONDS_TO_RADIANS
    l_prime = sun_anomoly(julian_century) * ARC_SECONDS_TO_RADIANS

    # Precompute reoccuring arguments
    f_omega_d = 2 * (F + omega - D)

    celestial_x += 1e-6 * ((
        -6844318.44 * np.sin(omega) - 523908.04 * np.sin(f_omega_d) 
        - 90552.22 * np.sin(2*(F+omega)) + 82168.76 * np.sin(2*omega)
        + 58707.02 * np.sin(l_prime)
    ) + julian_century * (
        205833.11 * np.cos(omega) + 12814.01  * np.cos(f_omega_d)
    ))

    celestial_y += 1e-6 * ((
        9205236.26 * np.cos(omega) + 573033.42 * np.cos(f_omega_d) 
        + 97846.69 * np.cos(2*(F+omega)) - 89618.24 * np.cos(2*omega)
        + 22438.42 * np.cos(l_prime-f_omega_d)
    ) + julian_century * (
        153041.79 * np.sin(omega) + 11714.49  * np.sin(f_omega_d)
    ))

    return celestial_x * ARC_SECONDS_TO_RADIANS, celestial_y * ARC_SECONDS_TO_RADIANS


def utc_time_to_julian_date(
    utc_time: datetime
) -> float:
    """Convert UTC time to Julian date.

    This calculation is only valid for days after March 1900.

    Args:
        utc_time (datetime): The observation time as a datetime object

    Returns:
        julian_date (float): The observation time as a julian date.
    """
    year, month, day = utc_time.year, utc_time.month, utc_time.day
    julian_date = (
        367 * year - 7 * (year + (month + 9) // 12) // 4 + 275 * month // 9 + day + 1721013.5
    )

    # update with the frational day
    julian_date +=  (
        utc_time.hour + utc_time.minute / 60 
        + (utc_time.second  + 1e-6 * utc_time.microsecond) / 3600 
    ) / 24

    return julian_date


## Precession Polynomials (Arc-Seconds)
# Equation 5.16
precession_x = np.polynomial.polynomial.Polynomial(
    [-0.016617, 2004.191898, -0.4297829, -0.19861834]
)
precession_y = np.polynomial.polynomial.Polynomial(
    [-0.006951,  -0.025896, -22.4072747, 0.00190059]
)

# Nutation Polynomials (Arc-Seconds)
# Equation 5.43

# Mean_anomaly of the moon (l)
moon_anomoly = np.polynomial.polynomial.Polynomial(
    [485868.249036, 1717915923.217800, 31.879200, 0.05163500]
)

# Mean_anomaly of the sun (l-prime)
sun_anomoly = np.polynomial.polynomial.Polynomial(
    [1287104.793048, 129596581.048100, -0.55320]
)

# Moon thing 1 (F)
moon_longitude = np.polynomial.polynomial.Polynomial(
    [335779.526232, 1739527262.8478, -12.7512]
)

# Moon elongation from sun (D)
moon_elongation = np.polynomial.polynomial.Polynomial(
    [1072260.703692, 1602961601.209000, -6.3706]
)

# Moon ascension node (Omega)
moon_ascension = np.polynomial.polynomial.Polynomial(
    [450160.398036, -6962890.5431, 7.4722]
)