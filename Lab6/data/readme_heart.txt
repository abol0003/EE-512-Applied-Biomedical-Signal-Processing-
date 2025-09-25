The files heart_#.dat contain three columns: 

1st column
----------
Time (in milliseconds) between successive heartbeats (also called RR-intervals).
It is the inverse of the cardiac rhythm.
The raw time series, by nature, is irregularly sampled. Here, through an interpolation,
(standard for these time series) it is regularly sampled. 

2nd column
----------
Average blood pressure at each heartbeat (in mm Hg). This time series is also regularly
re-sampled.

3rd column
----------
Instantaneous lung volume (ILV, arbitrary units), which measures the respiration.

All three signals have a sampling frequency of 4 Hz.


Cardiovascular regulation is mediated by the autonomous nerve system (ANS). RR-intervals and
blood pressure are of course correlated. Respiration influences both of them directly 
(at respiration frequency, around 0.25 Hz) and through ANS control loop. This second
effect gives rise to a spectral peak around 0.1 Hz.

A standard classification for cardiac signals is:

0 - 0.04 Hz   : Very Low Frequencies (VLF). Mainly influenced by hormonal regulation.
0.04 - 0.15 Hz: Low Frequencies (LF). Influenced by the sympathetic and parasympathetic
branches of the ANS. Sympathetic branch increases heart rate, parasympathetic branch 
decreases it.
0.15 - 0.4 Hz : High Frequencies (HF). Mainly influenced by parasympathetic branch.

Description of the signals
--------------------------

Files Heart_1.dat and Heart_2.dat are extracted from a protocol on the influence of
alcohol on ANS cardiac regulation.

Heart_1.dat: patient condition at rest. Respiration rate controled.
Heart_2.dat: same patient, but 40 minutes after alcohol ingestion. 
