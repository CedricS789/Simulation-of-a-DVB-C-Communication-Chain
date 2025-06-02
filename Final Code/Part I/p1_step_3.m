%%   =================== p1_step_2 - Simulation for BER Curve Generation ===================
clear; close all; clc;
addpath('p1_functions'); 


%% =================== Load Simulation Parameters  ===================
Nbps = 4;
params = initParameters(Nbps);
displayParameters(params);
NumBits = params.timing.NumBits;
ModType = params.modulation.ModulationType;
ModOrder= params.modulation.ModulationOrder;
OSF = params.sampling.OversamplingFactor;
SymRate = params.timing.SymbolRate;
Tsymb = params.timing.SymbolPeriod;
BitRate = params.timing.BitRate;
Fs = params.sampling.SamplingFrequency;
Ts = params.sampling.SamplePeriod;
Beta = params.filter.RolloffFactor;
NumTaps = params.filter.NumFilterTaps;

% --- BER Curve Parameters ---
EbN0_min_dB = params.simulation.EbN0_min_dB;            % Minimum Eb/N0 value in dB
EbN0_max_dB = params.simulation.EbN0_max_dB;            % Maximum Eb/N0 value in dB
EbN0_step_dB = params.simulation.EbN0_step_dB;          % Step size for Eb/N0 sweep in dB
iterations = params.simulation.iterations_per_EbN0;     % Iterations for averaging BER at each Eb/N0 point
EbN0_domain_dB = params.simulation.EbN0_domain_dB;      % Range of Eb/N0 values to simulate (dB)
num_EbN0_points = length(EbN0_domain_dB);               % Number of points on the BER curve


%% =================== Communication Chain ===================
% --- Transmitter  ---
bit_tx = randi([0, 1], 1, NumBits).';
symb_tx = mapping(bit_tx, Nbps, ModType);
symb_tx_up = upSampler(symb_tx, OSF).';
g_rrc = rrcFilter(Beta, SymRate, OSF, NumTaps);
signal_tx = applyFilter(symb_tx_up, g_rrc, NumTaps);
signalPower_tx = mean(abs(signal_tx).^2);
Eb = signalPower_tx / BitRate;

% -- Introduce Noise --
EbN0dB = 10;
signal_tx_noisy = addAWGN(signal_tx, Eb, EbN0dB, OSF, SymRate);

% --- Receiver Chain ---
signal_rx = applyFilter(signal_tx_noisy, g_rrc, NumTaps);
symb_rx_down = downSampler(signal_rx, OSF).';
bit_rx = demapping(symb_rx_down, Nbps, ModType); 
bit_rx = bit_rx(:);  


%% =================== Generate Plots  ===================
bits_to_plot = min(params.timing.NumBits, 100 * Nbps); 
plotConstellation_Tx_Rx(ModOrder, ModType, symb_tx, signal_tx_noisy);
plotConstellation_Tx_Rx(ModOrder, ModType, symb_tx, signal_rx);
plotConstellation_Tx_Rx(ModOrder, ModType, symb_tx, symb_rx_down);
plotBitstream_Tx_Rx(bit_tx, bit_rx, bits_to_plot);
plotPSD_Tx_Rx(signal_tx, signal_rx, Fs);
plotBasebandFrequencyResponse(signal_tx, signal_rx, Fs);



%% ========================= BER (Multiple Modulations) ======================
nbps_values = [1 2, 3, 4, 5 6, 7, 8];
g_rrc = rrcFilter(Beta, SymRate, OSF, NumTaps);
all_ber_data_modulation = zeros(num_EbN0_points, length(nbps_values));

for idx_nbps = 1:length(nbps_values)
    current_nbps = nbps_values(idx_nbps);
    params = initParameters(current_nbps);
    NumBits = params.timing.NumBits;
    ModType = params.modulation.ModulationType;
    ModOrder= params.modulation.ModulationOrder;
    OSF = params.sampling.OversamplingFactor;
    SymRate = params.timing.SymbolRate;
    Tsymb = params.timing.SymbolPeriod;
    BitRate = params.timing.BitRate;
    Fs = params.sampling.SamplingFrequency;
    Ts = params.sampling.SamplePeriod;
    Beta = params.filter.RolloffFactor;
    NumTaps = params.filter.NumFilterTaps;

    % --- BER Curve Parameters ---
    EbN0_min_dB = params.simulation.EbN0_min_dB;            % Minimum Eb/N0 value in dB
    EbN0_max_dB = params.simulation.EbN0_max_dB;            % Maximum Eb/N0 value in dB
    EbN0_step_dB = params.simulation.EbN0_step_dB;          % Step size for Eb/N0 sweep in dB
    iterations = params.simulation.iterations_per_EbN0;     % Iterations for averaging BER at each Eb/N0 point
    EbN0_domain_dB = params.simulation.EbN0_domain_dB;      % Range of Eb/N0 values to simulate (dB)
    num_EbN0_points = length(EbN0_domain_dB);               % Number of points on the BER curve

    current_mod_order = 2^current_nbps;

    ber_data_one_modulation = zeros(num_EbN0_points, 1);
    for idx_EbN0 = 1:num_EbN0_points
        EbN0dB = params.simulation.EbN0_domain_dB(idx_EbN0);
        bit_tx = randi([0, 1], 1, NumBits).'; 
        symb_tx = mapping(bit_tx, current_nbps, ModType);
        symb_tx = upSampler(symb_tx, OSF).'; 
        signal_tx = applyFilter(symb_tx, g_rrc, NumTaps); 
        signalPower = mean(abs(signal_tx).^2);

        current_bit_rate = SymRate * current_nbps; 
        Eb = signalPower / current_bit_rate;

        total_bit_errors_point = 0;
        total_bits_sim_point = 0;

        time_vector = (0:length(signal_tx)-1).' * Ts; 
        time_vector_symb = (0:length(symb_tx)-1).' * Tsymb; 

        for iter = 1:iterations 
            signal_tx_noisy = addAWGN(signal_tx, Eb, EbN0dB, OSF, SymRate); 
            signal_rx = applyFilter(signal_tx_noisy, g_rrc, NumTaps);

            symb_rx_down = downSampler(signal_rx, OSF); 
            bit_rx = demapping_v2(symb_rx_down, current_nbps, ModType); 

            num_errors_iter = sum(bit_tx ~= bit_rx);
            bits_iter = length(bit_tx);
            total_bit_errors_point = total_bit_errors_point + num_errors_iter;
            total_bits_sim_point = total_bits_sim_point + bits_iter;
        end
        ber_data_one_modulation(idx_EbN0) = total_bit_errors_point / total_bits_sim_point;
        fprintf('\n  Nbps = %d, Eb/N0 = %5.1f dB : Bits = %8d, Errors = %6d, BER = %.3e', ...
                  current_nbps, EbN0dB, total_bits_sim_point, total_bit_errors_point, ber_data_one_modulation(idx_EbN0));
    end
    all_ber_data_modulation(:, idx_nbps) = ber_data_one_modulation;
end
% Plot
plotBERCurvesModulation(all_ber_data_modulation, params, nbps_values);