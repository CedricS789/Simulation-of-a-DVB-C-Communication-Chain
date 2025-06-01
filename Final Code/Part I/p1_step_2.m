%% =================== p1_step_2 - Simulation over an Ideal Channel ===================
clear; close all; clc;
addpath('p1_functions');


%% =================== Load Simulation Parameters  ===================
Nbps = 6;
params = initParameters(Nbps);
displayParameters(params);
NumBits = params.timing.NumBits;
ModType = params.modulation.ModulationType;
ModOrder= params.modulation.ModulationOrder;
OSF = params.sampling.OversamplingFactor;
SymRate = params.timing.SymbolRate;
BitRate = params.timing.BitRate;
Fs = params.sampling.SamplingFrequency;
Ts = params.sampling.SamplePeriod;
Beta = params.filter.RolloffFactor;
NumTaps = params.filter.NumFilterTaps;
iterations = params.simulation.iterations_per_EbN0;



%% =================== Communication Chain ===================
% --- Transmitter  ---
bit_tx = randi([0, 1], 1, NumBits).';
symb_tx = mapping(bit_tx, Nbps, ModType);
symb_tx = upSampler(symb_tx, OSF).';
h_rrc = rrcFilter(Beta, SymRate, OSF, NumTaps);
signal_tx = applyFilter(symb_tx, h_rrc, NumTaps);
signalPower_tx = mean(abs(signal_tx).^2);
Eb = signalPower_tx / BitRate;

% -- Introduce Noise --
EbN0dB = 10;
signal_tx_noisy = addAWGN(signal_tx, Eb, EbN0dB, OSF, SymRate);

% --- Receiver Chain ---
signal_rx = applyFilter(signal_tx_noisy, h_rrc, NumTaps);
symb_rx = downSampler(signal_rx, OSF);
bit_rx = demapping(symb_rx, Nbps, ModType); 
bit_rx = bit_rx(:);  


%% =================== Generate Plots  ===================
bits_to_plot = min(params.timing.NumBits, 100 * Nbps); 
plotConstellation_Tx_Rx(ModOrder, ModType, symb_tx, symb_rx);
plotBitstream_Tx_Rx(bit_tx, bit_rx, bits_to_plot);
plotFilterCharacteristics(h_rrc, Beta, Fs, OSF);
plotPSD_Tx_Rx(signal_tx, signal_rx, Fs);
plotBasebandFrequencyResponse(signal_tx, signal_rx, Fs);
