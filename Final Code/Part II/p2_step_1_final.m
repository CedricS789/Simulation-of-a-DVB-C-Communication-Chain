%% =================== Step 1_1_final - Assessing the Impact of Synchronization Errors - CFO and Phase Offset ===================
%   Introduce Errors: The synchronization errors affect the signal as it's received,    
%   before any receiver processing. So, the point to introduce these errors mathematically  
%   is after the transmitter's pulse shaping (signal_tx)     
%   and before the receiver's matched filter (signal_rx)
%
%   Chain: Bits -> Map -> Upsample -> h_rrc filter-> signal_tx -> Apply Sync Errors -> Add AWGN -> signal_rx_filtered -> ...
%   
%   This code introduce a simple single CFO error to the transmitted signal, and plot of the constellation diagram
% =================================================================================================================================

clear; close all; clc;
addpath('../Part I - Optimal Communication chain over the ideal channel/p1_functions');
addpath('p2_functions')

%% ========================================== Load Simulation Parameters  ==========================================
Nbps    = 4;
params  = initParameters(Nbps);
NumBits = params.timing.NumBits;
ModType = params.modulation.ModulationType;
ModOrder= params.modulation.ModulationOrder;
OSF     = params.sampling.OversamplingFactor;
SymRate = params.timing.SymbolRate;
BitRate = params.timing.BitRate;
Fs      = params.sampling.SamplingFrequency;
Ts      = params.sampling.SamplePeriod;
Beta    = params.filter.RolloffFactor;
NumTaps = params.filter.NumFilterTaps;
iterations = params.simulation.iterations_per_EbN0;
displayParameters(params);

% ---- CFO Parameters ----
Fc = 600e6;                                 % Carrier frequency in Hz
delta_cfo_ppm   = 1 * 1e-6 * Fc;            % Frequency offset in Hz (1 ppm)
delta_omega     = 2 * pi * delta_cfo_ppm;    % Frequency offset in rad/s
phi_0           = 0;                        % Phase offset in rad


%% ========================================== Communication Chain ==========================================
% --- Transmitter  ---
bit_tx      = randi([0, 1], 1, NumBits).';
symb_tx     = mapping(bit_tx, Nbps, ModType);
symb_tx_up  = upSampler(symb_tx, OSF).';
h_rrc       = rrcFilter(Beta, SymRate, OSF, NumTaps);
signal_tx   = applyFilter(symb_tx_up, h_rrc, NumTaps);
signalPower_tx  = mean(abs(signal_tx).^2);
Eb              = signalPower_tx / BitRate;

% -- Introduce Noise --
EbN0dB     = 1000;                                                % Very Low Noise
signal_tx_noisy = addAWGN(signal_tx, Eb, EbN0dB, OSF, SymRate);

% --- Introduce CFO and phase offset ---
num_samples_tx  = length(signal_tx_noisy);                        % Number of samples in the transmitted signal
time_vector     = (0 : num_samples_tx - 1).' * Ts;                % The TA insisted on this
offset_signal   = exp(1j * (delta_omega * time_vector + phi_0));  % Create the offset signal
signal_tx_offset   = signal_tx_noisy .* offset_signal;            % Apply CFO to the transmitted signal

% --- Receiver Chain ---
signal_rx  = applyFilter(signal_tx_offset, h_rrc, NumTaps);
symb_rx    = downSampler(signal_rx, OSF).';
bit_rx     = demapping(symb_rx, Nbps, ModType); 
bit_rx     = bit_rx(:).';  


%% ====================== Generate Plots  =======================
bits_to_plot = min(params.timing.NumBits, 100 * Nbps); 
plotConstellation_Tx_Rx(ModOrder, ModType, symb_tx, symb_rx);
% plotBitstream_Tx_Rx(bit_tx, bit_rx, bits_to_plot);
% plotFilterCharacteristics(h_rrc, Beta, Fs, OSF);
% plotPSD_Tx_Rx(signal_tx, signal_rx, Fs);
% plotBasebandFrequencyResponse(signal_tx, signal_rx, Fs);
