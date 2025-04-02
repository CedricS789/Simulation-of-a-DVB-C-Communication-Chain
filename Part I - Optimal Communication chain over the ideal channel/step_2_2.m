%%   =================== Step 2.2 - Simulation with Fixed Noise Level ===================
%
%       Simulates the communication chain over an AWGN channel at a single,
%       fixed Eb/N0. Averages BER over multiple iterations for accuracy.
%
%   =====================================================================================

clear; close all; clc;

%% ========================== Simulation Parameters ============================
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
iterations_for_accuracy = params.simulation.iterations_per_EbN0; % Number of runs to average BER
EbN0dB                  = 100 ; % The single Eb/N0 value (dB) for this test

fprintf('\n\n------------------------------------');
fprintf('\nAccurate BER Test Setup:');
fprintf('\n------------------------------------');
fprintf('\n  Modulation: %d-%s, Eb/N0 = %.1f dB', ModOrder, upper(ModType), EbN0dB);
fprintf('\n  %d Bits per iteration, %d iterations', NumBits, iterations_for_accuracy);
fprintf('\n  Total bits simulated = %d', NumBits * iterations_for_accuracy);


%% ======================== Transmitter Chain (Executed Once) =====================
fprintf('\n  Executing Transmitter Chain ...');

% -------- 1. Generate Source Data Bits --------
bit_tx = randi([0, 1], 1, NumBits);

% -------- 2. Symbol Mapping --------
symb_tx = mapping(bit_tx, Nbps, ModType);

% -------- 3. Upsampling --------
symb_tx_up = upSampler(symb_tx, OSF).'; % Transpose for filter

% -------- 4. Pulse Shaping (Transmit Filtering) --------
fprintf('\n  Generating RRC filter...');
h_rrc = rrcFilter(Beta, SymRate, OSF, NumTaps); % Generates filter & its plots
signal_tx = applyFilter(symb_tx_up, h_rrc, NumTaps);

% -------- 5. Calculate Tx Signal Power and Eb (Once) --------
signalPower_tx = mean(abs(signal_tx).^2);
Eb_tx          = signalPower_tx / BitRate; % Energy per bit (Ptx / Rb)
fprintf('\n  Tx Eb = %.2e', Eb_tx);
fprintf('\n  Transmitter Processing Complete.');


%% ================= Simulation Loop for Averaging at Fixed Eb/N0 ==================
fprintf('\n  Starting %d iterations for averaging at Eb/N0 = %.1f dB...', iterations_for_accuracy, EbN0dB);

total_bit_errors = 0; % Accumulator for bit errors across iterations

for iter = 1:iterations_for_accuracy
    % --- 1. AWGN Channel ---
    % Add a new random noise instance each iteration
    signal_tx_noisy = addAWGN(signal_tx, Eb_tx, EbN0dB, Fs);

    % --- 2. Receiver Chain ---
    signal_rx_filtered = applyFilter(signal_tx_noisy, h_rrc, NumTaps); % Matched Filter
    symb_rx = downSampler(signal_rx_filtered, OSF).'; % Downsample
    bit_rx = demapping(symb_rx, Nbps, ModType); % Demap
    bit_rx = bit_rx(:).'; % Reshape to a vector

    % --- 3. Calculate and Accumulate Errors ---
    num_errors_iter = sum(bit_tx ~= bit_rx);
    total_bit_errors = total_bit_errors + num_errors_iter; % Accumulate errors

end % End of averaging loop

fprintf('\n  Averaging Loop Complete.');


%% ==================== Performance Evaluation (Averaged Result) =====================
fprintf('\n\n------------------------------------');
fprintf('\nFinal Performance (Averaged over %d runs at Eb/N0 = %.1f dB):', iterations_for_accuracy, EbN0dB);
fprintf('\n------------------------------------');

% -------- Calculate Average BER --------
total_bits_simulated = NumBits * iterations_for_accuracy;

if total_bits_simulated > 0
    ber_avg = total_bit_errors / total_bits_simulated;
else
    ber_avg = NaN;
    total_bit_errors = 0; % Ensure consistency if no bits simulated
end

fprintf('\n   Total bits simulated: %d', total_bits_simulated);
fprintf('\n   Total bit errors:     %d', total_bit_errors);
fprintf('\n   Average BER:          %.3e', ber_avg);


%% ====================== Generate Plots (Based on Last Iteration) =======================
% Plots show a snapshot from the *last* run. The BER value above is the average.
fprintf('\n\n------------------------------------');
fprintf('\nGenerating Plots (showing results from the LAST iteration)...');
fprintf('\n------------------------------------');

bits_to_plot = min(params.timing.NumBits, 100 * params.modulation.Nbps); % Example value

% --- Plot Constellations (Tx vs Rx from last iteration) ---
plotConstellation_Tx_Rx(ModOrder, ModType, symb_tx, symb_rx);

% --- Plot Bitstreams (Tx vs Rx from last iteration) ---
plotBitstream_Tx_Rx(bit_tx, bit_rx, bits_to_plot);

% --- Plot Power Spectral Density (PSD from last iteration) ---
plotPSD_Tx_Rx(signal_tx, signal_rx_filtered, Fs);

fprintf('\n\n------------------------------------');
fprintf('\nPlotting complete.');
fprintf('\n------------------------------------');