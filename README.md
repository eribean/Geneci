# Geneci

TLDR: ECI / ECEF conversion tool for nearish earth transformations.

Geneci (the 'g' is pronounced the same as in gif) is an ECI <----> ECEF conversion tool. The motivation  behind this package is a lightweight python implementation that provides *Good ENough ECI* transformations that account for a moderate amount of precession and nutation **without** the bulk of other conversion tools and a very simple interface. Many of the milli arc-second terms are dropped and polar motion is not included. In addition, it assumes UT1 == UTC in the earth rotational rate calculations. These approximations are often sufficient when looking from space (~1000km HAE) down toward earth. It is not intended for astronomical use since it lacks, on purpose, the necessary accuracy.

Velocity conversions are supported with the precession / nutation rotational time derivative contributions assumed to be zero. In other words, only the earth rotation is included in the time-derivative velocity calculation.

## Installation

A PyPi package is coming soon, but for now it can be built from source.

```sh
pip install .
```

## Usage

There are two primary functions in Geneci.

### Conversion from ECI to ECEF

The function ```eci_to_ecef``` takes an ECI position vector, an observation time, and an optional ECI velocity and converts into an ECEF position / velocity.

```python
from datetime import datetime, timedelta
import numpy as np
import geneci

# Choose a stationary position
HAE = 1000000
radius = 6378000 + HAE #  Earth radius plus 1000km

# Fake ECI positions
stationary = np.tile([radius * np.cos(np.pi / 4),  0, radius * np.sin(np.pi/4)], (100, 1))
velocity = np.zeros_like(stationary)

observation = np.array([datetime(year=2014, month=2, day=8) + timedelta(hours=offset)
                        for offset in np.linspace(0, 24, 100)])

ecef = np.asarray([geneci.eci_to_ecef(stationary[ndx], observation[ndx], velocity[ndx]) 
                   for ndx in range(100)])

# Velocity in ECEF is non-zero since the earth is spinning
ecef[2, 1]

# >>  array([-217.18217551, 311.70176875, 0.])
```

### Conversion from ECEF to ECI

The function ```ecef_to_eci``` takes an ECEF position vector, an observation time, and an optional ECEF velocity and converts into an ECI position / velocity. It is essentially the inverse of the above transformation with the same function signature.

### Rotation Matrix

The ```ecef_to_eci``` rotation matrix can be obtained by supplying a time in julian centries. This function is meant to be a helper function but is helpful if doing bulk transformations at the same observation time. The ```eci_to_ecef``` rotation matrix is the simply the transpose of the ```ecef_to_eci``` matrix.

```python
from datetime import datetime, timedelta

import geneci
from geneci.conversions import DERIVATIVE_MATRIX

# Convert UTC time to Julian Century
observation_time =  datetime(year=2021, month=3, day=19, hour=7, second=41)
julian_date = geneci.utc_time_to_julian_date(observation_time)
julian_century = (julian_date - 2451545.0) / 36525.0

# Compute the rotation matrix
ecef_to_eci_rotation = geneci.rotation_matrix_ecef_to_eci(julian_century)

# Time derivative of the rotation matrix
time_derivative = ecef_to_eci_rotation @ DERIVATIVE_MATRIX
```

## Reference

The calculations were taken from IERS 2010 convention. The document is found [here](https://iers-conventions.obspm.fr/content/tn36.pdf).

If you require higher fidelity calculations, try [pyerfa](https://github.com/liberfa/pyerfa).
