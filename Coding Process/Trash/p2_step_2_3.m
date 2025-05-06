%% ======= Step 2_v3 - Assessing the Impact of Synchronization Errors - Plotting BER =========
% 
% Applying a Modular approach to Synchronization errors Addition
% and Plotting the BER curve for different CFO and phase offset values using functions
%
%% ============================================================================================

clear; close all; clc;
addpath('../Part I - Optimal Communication chain over the ideal channel/p1_functions');
addpath('p2_functions');

% -- Load Simulation Parameters ---
Nbps    = 2;
params  = initParameters_v2(Nbps);
displayParameters(params);
    
% ---- CFO and SCO Parameters ----
Fc = 600e6;                                 % Carrier frequency in Hz
delta_cfo_ppm  = 0.09 * 1e-6 * Fc;            % Frequency offset in Hz (1 ppm)
phi_0 = 0;                              % Phase offset in rad
sample_time_offset = 0;                      % Sample offset (0 samples)

% --- Plots ---
ber_datas = generateBERDataWithSyncErrors(params, delta_cfo_ppm, phi_0, sample_time_offset);
plotBERCurve(ber_datas, params);