% =================== Master Script for Multi-BER Curve Simulation ===================
%
%   Orchestrates BER simulations for a range of modulation orders (defined by Nbps)
%   and plots all results on a single comparison figure.
%
% ========================================================================================

clear; close all; clc;

%% ========================================== Simulation Setup ==========================================

% --- Define the Nbps values (bits per symbol) to simulate ---
Nbps_values_to_simulate = [1, 2, 4, 6]; % Example: BPSK, 4-QAM, 16-QAM, 64-QAM
% Nbps_values_to_simulate = [1, 3, 5]; % Example: BPSK, 8-PAM, 32-PAM

num_simulations = length(Nbps_values_to_simulate);
all_results = cell(1, num_simulations); % Cell array to store results structs

fprintf('============================================\n');
fprintf('Starting Multi-BER Simulation Run\n');
fprintf('============================================\n');
fprintf('Nbps values to simulate: %s\n', mat2str(Nbps_values_to_simulate));
fprintf('--------------------------------------------\n');


%% ========================================== Simulation Loop ==========================================

for i = 1:num_simulations
    current_Nbps = Nbps_values_to_simulate(i);

    % --- 1. Get Parameters for this Nbps ---
    % It's good practice to get a fresh parameter set for each run
    params = coreParameters(current_Nbps);

    % --- Optional: Override specific common parameters here if needed ---
    % Example: Ensure all simulations use the exact same EbN0 range
    % params.simulation.EbN0_min_dB = -2;
    % params.simulation.EbN0_max_dB = 18;
    % params.simulation.EbN0_step_dB = 1;
    % params.simulation.iterations_per_EbN0 = 10; % Maybe increase iterations

    fprintf('\n[%d/%d] Running Simulation for Nbps = %d (%d-%s)\n', ...
            i, num_simulations, current_Nbps, params.modulation.ModulationOrder, upper(params.modulation.ModulationType));

    % --- 2. Run the Simulation for this Parameter Set ---
    [EbN0_dB, ber_sim, bits_sim, errors_acc] = simulateSingleBERCurve(params);

    % --- 3. Store Results ---
    results_struct = struct();
    results_struct.EbN0_dB = EbN0_dB;
    results_struct.BER_sim = ber_sim;
    results_struct.params = params; % Store the params used for this run
    results_struct.total_bits = bits_sim; % Optional: store for analysis
    results_struct.total_errors = errors_acc; % Optional: store for analysis

    all_results{i} = results_struct;

    fprintf('Finished simulation for Nbps = %d.\n', current_Nbps);
    fprintf('--------------------------------------------\n');

end % End loop over Nbps values

fprintf('\n============================================\n');
fprintf('All Simulations Complete.\n');
fprintf('============================================\n');


%% ========================================== Plotting Results ==========================================

fprintf('\nPlotting combined BER curves...\n');

if ~isempty(all_results)
    plotMultiBERCurves(all_results);
    fprintf('Plotting complete.\n');
else
    fprintf('No results were generated to plot.\n');
end

fprintf('============================================\n');