%% ======= Step 1_2 - Assessing the Impact of Synchronization Errors - BER Comparison =========
% 
% Applying a Modular approach to Synchronization errors Addition
% and Plotting the BER curve for different CFO and phase offset values using function
%
%% ============================================================================================

clear; close all; clc;
addpath('../Part I - Optimal Communication chain over the ideal channel/functions');
addpath('functions');
Nbps    = 4;
params  = initParameters(Nbps);
displayParameters(params);
ber_datas = generateBERDataOffset(params, 10, 0);
plotBERCurve(ber_datas, params);