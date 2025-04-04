%% =================== Simulate and Plot BER for Multiple Nbps ===================
%
%   Purpose: This script orchestrates the simulation of a digital communication
%   chain (like DVB-C) for several different modulation orders (defined by Nbps)
%   over an AWGN channel. It generates Bit Error Rate (BER) data for each Nbps
%   across a range of Eb/N0 values and plots all the resulting BER curves
%   on a single graph for comparison against theoretical predictions.
%
%   Workflow:
%       1. Define the Nbps values to simulate (e.g., 2 for QPSK, 4 for 16-QAM).
%       2. Loop through each Nbps value:
%           a. Initialize simulation parameters using `initParameters`.
%           b. Run the BER simulation using `generateBERData`.
%           c. Store the parameters and BER results.
%       3. Call `plotMultipleBERCurves` to plot all stored results together.
%
%   Outputs:
%       - Console output showing simulation progress for each Nbps.
%       - A figure containing BER vs. Eb/N0 plots for all simulated Nbps values,
%         comparing simulated points against theoretical curves.
%
%
%   Chain: Bits -> Map -> Upsample -> Tx RRC Filter -> (Ideal Channel + Noise) -> Rx RRC Filter -> Downsample -> Demap -> Bits
% ==============================================================================

clear; close all; clc;
addpath('functions'); % Ensure your functions folder is in the path

%% ========================== Simulation Setup ============================
% Define the Nbps values to simulate 
nbps_to_simulate = [1 2 3 4 5 6 7 8 9 10 11 12]; % Nbps values to simulate
fprintf('=============================================================================================================================================\n');
fprintf('Starting Multi-Nbps BER Simulation\n');
fprintf('Nbps values to simulate: %s\n', mat2str(nbps_to_simulate));
fprintf('=============================================================================================================================================\n\n');


%% ========================== Data Storage ============================
num_simulations = length(nbps_to_simulate);
all_ber_results = cell(1, num_simulations);      % To store BER data vectors
all_params      = cell(1, num_simulations);      % To store the params struct for each run


%% ========================== Simulation Loop ============================
for i = 1:num_simulations
    current_nbps = nbps_to_simulate(i);
    fprintf('==============================================================================================\n');
    fprintf('==== Starting Simulation for Nbps = %d ====\n', current_nbps);
    fprintf('==============================================================================================\n');
    params = initParameters(current_nbps);      % Initialize parameters for the current Nbps
    ber_data = generateBERData(params);         % Generate BER data for the current Nbps
    
    % --- Store Results ---
    all_ber_results{i}  = ber_data;       % We use cell arrays to store results because the length of the BER data may vary
    all_params{i}       = params;         % We use cell arrays to store the params struct for each run because it contains various fields
    
end

fprintf('===============================================\n');
fprintf('All simulations complete. Generating combined plot...\n');
fprintf('===============================================\n\n');

%% ========================== Plotting Combined Results ============================
% Plot all BER results on a single graph for comparison
hFig = plotMultipleBERCurves(nbps_to_simulate, all_ber_results, all_params);

fprintf('===============================================\n');
fprintf('Script Finished.\n');
fprintf('===============================================\n');