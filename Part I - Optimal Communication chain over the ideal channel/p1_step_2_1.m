%%   =================== Step 2_1 - Simulation over an Ideal Channel (No Noise) ===================
%
%       Purpose: This script simulates a DVB-C communication chain (transmitter to receiver) 
%       over an ideal channel without noise. It implements symbol mapping, 
%       Root Raised Cosine (RRC) pulse shaping, and matched filtering to verify the chain’s 
%       functionality and filter properties, such as inter-symbol interference (ISI) cancellation.
%
%       Outputs: Generates plots for constellation diagrams, bit streams, filter characteristics, 
%       and power spectral density (PSD) to visualize the signal at key stages.
%
%   ======================================================================================

clear; close all; clc;
addpath('functions'); 

%% ========================================== Simulation Parameters  ==========================================
Nbps    = 1;
params  = initParameters(Nbps);
NumBits = params.timing.NumBits;
ModType = params.modulation.ModulationType;
ModOrder = params.modulation.ModulationOrder;
OSF     = params.sampling.OversamplingFactor;
SymRate = params.timing.SymbolRate;
Fs      = params.sampling.SamplingFrequency;
Ts      = params.sampling.SamplePeriod;
Beta    = params.filter.RolloffFactor;
NumTaps = params.filter.NumFilterTaps;
displayParameters(params);


%% ========================================== Transmitter Chain ==========================================
fprintf('\n\n------------------------------------');
fprintf('\nTransmitter Processing...');
fprintf('\n------------------------------------');

% -------- 1. Generate Source Data Bits --------
bit_tx = randi([0, 1], 1, NumBits).';
fprintf('\n  Generated %d random bits.', NumBits);

% -------- 2. Symbol Mapping --------
symb_tx = mapping(bit_tx, Nbps, ModType);
fprintf('\n  Mapped bits to %d %d-%s symbols.', length(symb_tx), ModOrder, upper(ModType));

% -------- 3. Upsampling --------
symb_tx_up = upSampler(symb_tx, OSF).'; % Transpose for filter
fprintf('\n  Upsampled symbols by a factor of %d.', OSF);

% -------- 4. Pulse Shaping (Transmit Filtering) --------
h_rrc = rrcFilter(Beta, SymRate, OSF, NumTaps);
fprintf('\n  Generated RRC filter (beta=%.2f, %d taps).', Beta, NumTaps);
signal_tx = applyFilter(symb_tx_up, h_rrc, NumTaps);
fprintf('\n  Applied RRC pulse shaping filter.');


%% ========================================== Ideal Channel ==========================================
fprintf('\n\n------------------------------------');
fprintf('\nSimulating Ideal Channel (No Noise or Distortion).');
fprintf('\n------------------------------------');
signal_channel_output = signal_tx; % Output = Input


%% ========================================== Receiver Chain ==========================================
fprintf('\n\n------------------------------------');
fprintf('\nReceiver Processing...');
fprintf('\n------------------------------------');

% -------- 1. Matched Filtering --------
signal_rx = applyFilter(signal_channel_output, h_rrc, NumTaps); % Use same RRC filter
fprintf('\n  Applied RRC matched filter.');

% -------- 2. Downsampling (Symbol Timing Recovery) --------
symb_rx = downSampler(signal_rx, OSF).'; % Transpose back to row vector
fprintf('\n  Downsampled filtered signal by %d to recover symbols.', OSF);

% -------- 3. Symbol Demapping --------
bit_rx = demapping_v2(symb_rx, Nbps, ModType); % Transpose back to row vector
bit_rx = bit_rx(:).'; % Reshape to a vector
fprintf('\n  Demapped received symbols back to bits.');


%% ========================================== Plots ==========================================
bits_to_plot = min(params.timing.NumBits, 100 * params.modulation.Nbps);    % If NumBits is too large, plot only 100*Nbps bits

plotConstellation_Tx_Rx(ModOrder, ModType, symb_tx, symb_rx);
plotBitstream_Tx_Rx(bit_tx, bit_rx, bits_to_plot);
plotFilterCharacteristics(h_rrc, Beta, Fs, OSF); % Visualize RRC and combined RC pulse
plotPSD_Tx_Rx(signal_tx, signal_rx, Fs);

fprintf('\n\n------------------------------------');
fprintf('\nPlotting complete.');
fprintf('\n------------------------------------');