This folder contains the Matlab and OpenBUGS scripts and functions, as well the necessary data files, to implement the
parameter estimation methods in:

Colenutt et al. (2020) Quantifying the transmission of foot-and-mouth disease virus in cattle via a contaminated
environment. mBio 11, e00381-20 (https://doi.org/10.1128/mBio.00381-20)

MATLAB REQUIREMENTS AND CODE (see Text S1.1-1.4)
The scripts/functions were run using Matlab version 2019b and require the Statistics and Machine Learning and
Parallel Computing toolboxes. However, they can be easily adapted to run without the Parallel Computing toolbox
by changing the "parfor" loop in the ParEst function to a "for" loop.

ParEst.m - loads the data, implements the adaptive Metropolis scheme for each model/parameterisation
           and computes the DIC
Lhood.m - computes the log likelihood and prior for the input parameters

EnvironmentalTransmissionData.mat - Matlab data file containing the data analysed by ParEst.m; note the same
                                    data are also provided in the Excel file Data Set S1


OpenBUGS CODE (see Text S1.6)
EstimateDecayRates_qPCR.odc - OpenBUGS (version 3.2.3) script to estimated decay rates for FMDV RNA in environmental
                              samples
