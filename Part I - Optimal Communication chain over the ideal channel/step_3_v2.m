%%   =================== Step 3 v2 - Simulation for BER Curve Generation Using Functions ===================
%
%       Purpose: This script generates a BER curve for the communication chain over an AWGN 
%       channel across a range of Eb/N0 values, as part of the. It refines the simulation 
%       by modularizing the code into functions for improved organization and reusability.
%
%       Context: simulating the same DVB-C-like chain as previous 
%       iterations (symbol mapping, RRC filtering, AWGN). The use of functions (e.g., 
%       generateBERData, plotBERCurve) enhances maintainability while producing identical BER 
%       results to Step 3 v1. It maintains the 5 Msymb/s rate and 0.2 roll-off factor specified 
%       in the project.
%
%       Relation to Other Iterations: Evolves from Step 3 v1 by restructuring the simulation 
%       into modular components, retaining the same core logic. It builds on Step 2 v1’s 
%       baseline and Step 2 v2’s noise handling. The modularity 
%       distinguishes it from v1’s monolithic design.
%
%       Outputs: Produces a BER curve plot comparing simulated and theoretical performance, 
%       consistent with earlier iterations but with a cleaner implementation.
%
%   ====================================================================================================

clear; close all; clc;
addpath('functions');
Nbps    = 4;
params  = initParameters(Nbps);
displayParameters(params);
ber_datas = generateBERData(params);
plotBERCurve(ber_datas, params);