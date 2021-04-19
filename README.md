# sth-AMIS
Runs Adaptive Multiple Importance Sampling algorithm for STH (Ascaris) and resamples to create parameter and simulation output files for each IU.

# required input files
1) maps data

2) scen and group reference file

3) python code of model- multiple files

4) prior files- data and scripts

5) also need to set up appropriate folders for output

# Possible improvements to code
1) To improve the efficiency of the code, could parallelize weighting step.  This step takes a long time.
 
2) To improve quality of AMIS algorithm, could add cumulative distribution difference as stopping criteria and for adapting proposal distribution.
