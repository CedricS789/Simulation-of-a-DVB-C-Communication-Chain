%% =================== Step 1_1 - Assessing the Impact of Synchronization Errors - Assessing the effects of CFO ===================
%   Introduce Errors: The synchronization errors affect the signal as it's received,    
%   before any receiver processing. So, the point to introduce these errors mathematically  
%   is after the transmitter's pulse shaping (signal_tx)     
%   and before the receiver's matched filter (signal_rx)
%
%   Chain: Bits -> Map -> Upsample -> h_rrc filter-> signal_tx -> Apply Sync Errors -> Add AWGN -> signal_rx_filtered -> ...
%   Metrics: Constellation diagram, and BER
%   
%   This code introduce a simple single CFO error to the transmitted signal, and plot of the constellation diagram
% =================================================================================================================================

clear; close all; clc;
addpath('../Part I - Optimal Communication chain over the ideal channel/functions');
addpath('functions');


%% ========================================== Simulation Parameters  ==========================================
Nbps    = 4;
params  = initParameters(Nbps);
NumBits = params.timing.NumBits;
ModType = params.modulation.ModulationType;
ModOrder = params.modulation.ModulationOrder;
OSF     = params.sampling.OversamplingFactor;
SymRate = params.timing.SymbolRate;
BitRate     = params.timing.BitRate;
Fs      = params.sampling.SamplingFrequency;
Ts      = params.sampling.SamplePeriod;
Beta    = params.filter.RolloffFactor;
NumTaps = params.filter.NumFilterTaps;
iterations = params.simulation.iterations_per_EbN0;
displayParameters(params);


%% ========================================== Chain ==========================================
% --- Transmitter Chain ---
bit_tx      = randi([0, 1], 1, NumBits).';
symb_tx     = mapping(bit_tx, Nbps, ModType);
symb_tx_up  = upSampler(symb_tx, OSF).';
h_rrc       = rrcFilter(Beta, SymRate, OSF, NumTaps);
signal_tx   = applyFilter(symb_tx_up, h_rrc, NumTaps);
signalPower_tx  = mean(abs(signal_tx).^2);
Eb              = signalPower_tx / BitRate;

% ---- CFO Parameters and Applying CFO ----
delta_cfo_hz    = 1000;                                     % Frequency offset in Hz
delta_omega     = 2 * pi * delta_cfo_hz;                    % Frequency offset in rad/s
num_samples_tx  = length(signal_tx);                        % Number of samples in the transmitted signal
time_vector     = (0 : num_samples_tx - 1) * Ts;            % The TA insisted on this
cfo_shift       = exp(1j * delta_omega * time_vector);      % CFO signal
signal_tx_cfo   = signal_tx .* cfo_shift;                   % Apply CFO to the transmitted signal

% -- Introduce Noise --
EbN0dB     = 50 ;
signal_tx_noisy = addAWGN(signal_tx, Eb, EbN0dB, OSF, SymRate);

% --- Receiver Chain ---
signal_rx  = applyFilter(signal_tx_noisy, h_rrc, NumTaps);
symb_rx    = downSampler(signal_rx, OSF).';
bit_rx     = demapping_v2(symb_rx, Nbps, ModType); 
bit_rx     = bit_rx(:).';  

%% ====================== Generate Plots (Based on Last Iteration) =======================
bits_to_plot = min(params.timing.NumBits, 100 * Nbps); 
plotConstellation_Tx_Rx(ModOrder, ModType, symb_tx, symb_rx);
plotBitstream_Tx_Rx(bit_tx, bit_rx, bits_to_plot);