# Multi_Pulse_Complex_Ginzburg_Landau_Equation_2D
Codes for the paper "Dynamic interactions and equilibrium configurations of pulses in the two-dimensional complex quintic Ginzburg-Landau equation" by M.R. Turner and D.J.B. Lloyd 2026

To run the codes you will need to first install the chebfun libraries for MATLAB at 
https://www.chebfun.org/download/
The two main codes to run are
pulsecode_initial_pulse.m
main_N_pulse_run.m

pulsecode_initial_pulse.m  - Generates the data file for the form of the initial pulse. The key parameter to vary in this code is br=Real(beta). 
** run this code first to generate the pulsedata2.mat file which contains the initial pulse information to use in main_N_pulse-run.m.

main_N_pulse_run.m - Uses the form of the initial pulse adnt time integrates the projected system. For this N=2 case the main parameter to vary r0, which is the stating separation distance of the pulses.
