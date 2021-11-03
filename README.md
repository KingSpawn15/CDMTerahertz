
## Data types

1.	‘pwd\Polarization\111_STEM_x1.5k_StaticSpot_10nJ_wTimeScan20210617061959\20210617061959_1_PAngle=0000_td=188495.dm4’

111 is the crystal orientation. 111 gives OR, 1-10 gives only Photo-Dember.
PAngle is the waveplate angle before the laser enters the microscope.
(Theta_pol = 2*theta_WavePlate)
td is the time delay between laser and electron beam. 188495 mean 1884.95 ps
the absolute value is not important.

2.	‘pwd\111 Power scans\30nJ 58deg\20210620010207_001.dm4’
Scan over laser pulse energy (30nJ). Waveplate angle is const at 58 degrees.
Each folder contains 3-5 iterations of the SAME experiment (001-005)
You also have the individual files for specific time delays (not really needed).

3.	‘pwd\111 Pulse duration\SCMP370_50fs\20210618113507_001.dm4’
Changing the laser pulse dispersion (50-350fs) (not the spectral band-width) and acquiring time scans like in item no. 2.
This measurement doesn’t change the frequency components inside the pulse, so at this time using the results doesn’t make much sense. Need to double-check that.