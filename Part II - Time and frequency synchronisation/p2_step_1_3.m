%% ======= Step 1_3 - Assessing the Impact of Synchronization Errors - BER Comparison =========
% 
% Applying a Modular approach to Synchronization errors Addition
% and Plotting the BER curve for different CFO and phase offset values using functions
%
%% ============================================================================================

clear; close all; clc;
addpath('../Part I - Optimal Communication chain over the ideal channel/functions');
addpath('functions');

Nbps    = 4;
params  = initParameters_v2(Nbps);
displayParameters(params);
ber_datas = generateBERData_v2(params, 20, 0);
plotBERCurve(ber_datas, params);