This folder contains the Matlab and OpenBUGS scripts and functions, as well the necessary data files, to implement the parameter
estimation methods in Colenutt et al. "Quantifying the transmission of foot-and-mouth disease virus in cattle via a contaminated
environment" (see TextS1)

MATLAB REQUIREMENTS AND CODE (see Text S1.1-1.4)
The scripts/functions were run using Matlab version 2019b and require the Statistics and Machine Learning and
Parallel Computing toolboxes. However, they can be easily adapted to run without the Parallel Computing toolbox
by changing the "parfor" loop in the ParEst function to a "for" loop.

ParEst.m - loads the data, implements the adaptive Metropolis scheme for each model/parameterisation
           and computes the DIC
Lhood.m - computes the log likelihood and prior for the input parameters

EnvironmentalTransmissionData.mat - Matlab data file containing the data analysed by ParEst.m; note the same data are also
                                    provided in the Excel file DataS1.xlsx


OpenBUGS CODE (see Text S1.6)

EstimateDecayRates_qPCR.odc - OpenBUGS (version 3.2.3) script to estimated decay rates for FMDV RNA in environmental samples