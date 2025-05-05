%%   =================== Step 2 - Simulation of a Communication Chain over a Noisy Channel ===================
%
%       This script simulates a complete communication chain, including the
%       transmitter, channel, and receiver. The simulation uses a random bitstream
%       as input, which is modulated, transmitted through a noisy channel, and
%       received.
%       The script also includes plotting functions to visualize the results.
%
%   ===========================================================================================================

clear; close all; clc;

%% ========================== Simulation Parameters ============================
fprintf('Initializing Simulation Parameters...\n');
params = initializeParameters(); % Initialize parameters

% --- Extract parameters ---
Nbps        = params.modulation.Nbps;
NumBits     = params.timing.NumBits;
ModType     = params.modulation.ModulationType;
ModOrder    = params.modulation.ModulationOrder;
OSF         = params.sampling.OversamplingFactor;
SymRate     = params.timing.SymbolRate;
Fs          = params.sampling.SamplingFrequency;
Beta        = params.filter.RolloffFactor;
NumTaps     = params.filter.NumFilterTaps;
BitRate     = params.timing.BitRate;

% --- Define Fixed Eb/N0 and Averaging Iterations ---
iterations_for_accuracy = params.simulation.iterations_per_EbN0; 
EbN0dB                  = 10 ; % Set the single, fixed Eb/N0 value (dB)

fprintf('------------------------------------\n');
fprintf('Accurate BER Test Setup:\n');
fprintf('  Modulation: %d-%s, Eb/N0 = %.1f dB\n', ModOrder, upper(ModType), EbN0dB);
fprintf('  %d Bits per iteration, %d iterations\n', NumBits, iterations_for_accuracy);
fprintf('  Total bits simulated = %d\n', NumBits * iterations_for_accuracy);
fprintf('------------------------------------\n');

%% ======================== Transmitter Chain (Executed Once) =====================
fprintf('\nExecuting Transmitter Chain (Once)...\n');

% -------- 1. Generate Source Data Bits --------
bit_tx = randi([0, 1], 1, NumBits);

% -------- 2. Symbol Mapping --------
symb_tx = mapping(bit_tx, Nbps, ModType);

% -------- 3. Upsampling --------
symb_tx_up = upSampler(symb_tx, OSF).'; % Transpose for filter

% -------- 4. Pulse Shaping (Transmit Filtering) --------
fprintf('  Generating RRC filter (plots characteristics)...\n');
h_rrc = rrcFilter(Beta, SymRate, OSF, NumTaps); % Generates filter & its plots
signal_tx = applyFilter(symb_tx_up, h_rrc, NumTaps);

% -------- 5. Calculate Tx Signal Power and Eb (Once) --------
signalPower_tx = mean(abs(signal_tx).^2);
Eb_tx          = signalPower_tx / BitRate;
fprintf('  Tx Eb = %.2e\n', Eb_tx);
fprintf('Transmitter Processing Complete.\n');

%% ================= Simulation Loop for Averaging at Fixed Eb/N0 ==================
fprintf('\nStarting %d iterations for averaging at Eb/N0 = %.1f dB...\n', iterations_for_accuracy, EbN0dB);

total_bit_errors = 0; % Initialize error accumulator

for iter = 1:iterations_for_accuracy
    % --- 1. AWGN Channel ---
    % Add a new random noise instance each iteration
    signal_tx_noisy = addAWGN(signal_tx, Eb_tx, EbN0dB, Fs);

    % --- 2. Receiver Chain ---
    signal_rx_filtered = applyFilter(signal_tx_noisy, h_rrc, NumTaps);
    symb_rx = downSampler(signal_rx_filtered, OSF).';
    bit_rx = demapping(symb_rx, Nbps, ModType)';

    % --- 3. Calculate and Accumulate Errors ---
    num_errors_iter = sum(bit_tx ~= bit_rx);
    total_bit_errors = total_bit_errors + num_errors_iter;

end % End of averaging loop

fprintf('Averaging Loop Complete.\n');

%% ==================== Performance Evaluation (Averaged Result) =====================
fprintf('\n------------------------------------\n');
fprintf('Final Performance (Averaged over %d runs at Eb/N0 = %.1f dB):\n', iterations_for_accuracy, EbN0dB);
fprintf('------------------------------------\n');

% -------- Calculate Average BER --------
total_bits_simulated = NumBits * iterations_for_accuracy;

if total_bits_simulated > 0
    ber_avg = total_bit_errors / total_bits_simulated;
else
    ber_avg = NaN;
    total_bit_errors = 0;
end

fprintf('   Total bits simulated: %d\n', total_bits_simulated);
fprintf('   Total bit errors:   %d\n', total_bit_errors);
fprintf('   Average BER:        %.3e\n', ber_avg);

%% ====================== Generate Plots (Based on Last Iteration) =======================
% Plots show a snapshot from the *last* run. The BER value above is the average.
fprintf('\nGenerating Plots (showing results from the LAST iteration)...\n');

% --- Plotting Setup ---
plotSetup = configPlots(params);

% --- Plot Constellations (Tx vs Rx from last iteration) ---
plotConstellation_Tx_Rx(ModOrder, ModType, symb_tx, symb_rx, plotSetup);

% --- Plot Bitstreams (Tx vs Rx from last iteration) ---
plotBitstream_Tx_Rx(bit_tx, bit_rx, plotSetup.bits_to_plot, plotSetup);

% --- Plot Power Spectral Density (PSD from last iteration) ---
plotPSD_Tx_Rx(signal_tx, signal_rx_filtered, Fs, plotSetup);

fprintf('\nPlotting complete. Simulation finished.\n');
fprintf('------------------------------------\n');
