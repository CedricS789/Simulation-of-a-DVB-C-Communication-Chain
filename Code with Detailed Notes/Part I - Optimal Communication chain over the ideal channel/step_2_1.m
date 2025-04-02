%%   =================== Step 2 - Simulation of a Communication Chain over an Ideal Channel ===================
%
%       This script simulates a complete communication chain, including the
%       transmitter, channel, and receiver. The simulation uses a random bitstream
%       as input, which is modulated, transmitted through an ideal channel, and
%       received.
%       The script also includes plotting functions to visualize the results.
%
%   ===========================================================================================================

clear; close all; clc;


%% ========================================== Simulation Parameters  ==========================================
params = initializeParameters();
Nbps = params.modulation.Nbps;
NumBits = params.timing.NumBits;
ModType = params.modulation.ModulationType;
ModOrder = params.modulation.ModulationOrder;
OSF = params.sampling.OversamplingFactor;
SymRate = params.timing.SymbolRate;
Fs = params.sampling.SamplingFrequency;
Ts = params.sampling.SamplePeriod;
Beta = params.filter.RolloffFactor;
NumTaps = params.filter.NumFilterTaps;



%% ========================================== Transmitter Chain ==========================================
fprintf('\n\n------------------------------------\n');
fprintf('Transmitter Processing...\n');
fprintf('------------------------------------\n');

% -------- 1. Generate Source Data Bits --------
bit_tx = randi([0, 1], 1, NumBits);
fprintf('  Generated %d random bits.\n', NumBits);

% -------- 2. Symbol Mapping --------
symb_tx = mapping(bit_tx, Nbps, ModType);
fprintf('  Mapped bits to %d %d-%s symbols.\n', length(symb_tx), ModOrder, upper(ModType));

% -------- 3. Upsampling --------
symb_tx_up = upSampler(symb_tx, OSF).'; % Transpose for filter conv
fprintf('  Upsampled symbols by a factor of %d.\n', OSF);

% -------- 4. Pulse Shaping (Transmit Filtering) --------
% Generate Time-domain impulse response of the RRC Filter
h_rrc = rrcFilter(Beta, SymRate, OSF, NumTaps);
fprintf('  Generated RRC filter (beta=%.2f, %d taps).\n', Beta, NumTaps);
% Apply Filter using the dedicated function (handles conv + delay)
signal_tx = applyFilter(symb_tx_up, h_rrc, NumTaps);
fprintf('  Applied RRC pulse shaping filter.\n');



%% ========================================== Ideal Channel ==========================================
fprintf('\n\n------------------------------------\n');
fprintf('Simulating Ideal Channel (No Noise or Distortion).\n');
fprintf('------------------------------------\n');
signal_channel_output = signal_tx;



%% ========================================== Receiver Chain ==========================================
fprintf('\n\n------------------------------------\n');
fprintf('Receiver Processing...\n');
fprintf('------------------------------------\n');


% -------- 1. Matched Filtering --------
signal_rx = applyFilter(signal_channel_output, h_rrc, NumTaps);
fprintf('  Applied RRC matched filter.\n');

% -------- 2. Downsampling (Symbol Timing Recovery) --------
symb_rx = downSampler(signal_rx, OSF).'; % Transpose to row vector
fprintf('  Downsampled filtered signal by %d to recover symbols.\n', OSF);

% -------- 3. Symbol Demapping --------
bit_rx = demapping(symb_rx, Nbps, ModType)'; % Transpose to row vector
fprintf('  Demapped received symbols back to bits.\n');



%% ========================================== Performance Evaluation ==========================================
% Calculate Bit Error Rate

fprintf('\n\n------------------------------------\n');
fprintf('Simulation Complete.\n');
fprintf('------------------------------------\n');



%% ========================================== Plots ==========================================
% --- Plotting Setup ---
bits_to_plot = min(params.timing.NumBits, 100 * params.modulation.Nbps);
plotConstellation_Tx_Rx(ModOrder, ModType, symb_tx, symb_rx);
plotBitstream_Tx_Rx(bit_tx, bit_rx, bits_to_plot);
plotFilterCharacteristics(h_rrc, Beta, Fs, OSF);
plotPSD_Tx_Rx(signal_tx, signal_rx, Fs);

fprintf('\n\n------------------------------------\n');
fprintf('Plotting complete.\n');
fprintf('------------------------------------\n');
