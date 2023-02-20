# Geneci

TLDR: ECI / ECEF conversion tool for nearish earth calculations.

Geneci (the 'g' is pronounced the same as in gif) is an ECI <----> ECEF conversion tool.
The motivation behind this package is a lightweight python implementation that provides
*Good ENough ECI* transformations that account for a moderate amount of precession and
nutation **without** the bulk of the other conversion tools. Many of the milli arc-second terms
are dropped and polar motion is not included. In addition, it assumes UT1 == UTC in the earth
rotational rate calculations. These approximations are often sufficient when looking from
space (~1000km HAE) down toward earth. It is not intended for astronomical use since it lacks,
on purpose, the necessary accuracy.

Velocity conversions are supported with the precession / nutation rotational time derivative
contributions assumed to be zero. In other words, only the earth rotation is included in
the time-derivative velocity calculation.

## Installation

A PyPi package is coming soon, but for now it can be built from source.

```sh
pip install .
```

## Usage

There are two primary functions in Geneci.

1. ecef_to_eci
2. eci_to_ecef
