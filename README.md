# Cold-atoms-calculator
A calculator for use in a cold atom experiment. Can compute: Zeeman+hyperfine structure, AC polarizability, magnetic and optical traps, optical lattices, elastic and inelastic collision rates in a trapped gas, and more.


Written from about 2022 to 2024, while I was working on the KCs experiment in Innsbruck. It has been a long time since I checked all the calculations, so there may be some bugs or incorrect outputs. I'd recommend checking outputs with your own calculations before relying on it too much. 


The scattering data tab can load files containing scattering length a or inelastic scattering coefficient k2 versus magnetic field B. It reads files which start with the isotope: e.g. '39K_*.txt'. These files should have a structure as follows:

Pinpointing Feshbach resonances and testing Efimov universalities in 39 K, PRR5, 013174&linear
B (G)&a (a0)
1.000100000000000017e-01 -4.554707963893424250e+01
1.100099999999999967e-01 -4.555665393661984552e+01
.....

The first line should contain the source for the data or Feshbach resonance parameters used to produce the data, followed by the specification of the y-scale. Use 'linear' for scattering lengths and 'log' for inelastic collision coefficients.
The second line contains the axis labels. The rest of the file contains the data.


