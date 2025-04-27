%% ======= Step 1_3 - Assessing the Impact of Synchronization Errors - BER Comparison =========
% 
% Applying a Modular approach to Synchronization errors Addition
% and Plotting the BER curve for different CFO and phase offset values using functions
%
%% ============================================================================================

clear; close all; clc;
addpath('../Part I - Optimal Communication chain over the ideal channel/functions');
addpath('functions');

% -- Load Simulation Parameters ---
Nbps    = 4;
params  = initParameters_v2(Nbps);
displayParameters(params);

% ---- CFO and SCO Parameters ----
Fc = 600e6;                                 % Carrier frequency in Hz
delta_cfo_hz    = 1 * 1e-6 * Fc;            % Frequency offset in Hz (1 ppm)
phi_0           = 0;                        % Phase offset in rad
sco_samples_shift = 0;                      % Sample offset (0 samples)

% --- Plots ---
ber_datas = generateBERData_v2(params, delta_cfo_hz, phi_0, sco_samples_shift);
plotBERCurve(ber_datas, params);