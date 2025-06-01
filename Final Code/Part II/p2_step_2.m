%% ========================= p2_step_2 - CFO, Phase Offset, Time offset  =========================Add commentMore actions
close all; clear; clc;
addpath('../Part I/p1_functions');
addpath('p2_functions')


%% ========================= Load Simulation Parameters  ==========================================
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
EbN0_min_dB         = params.simulation.EbN0_min_dB;                % Minimum Eb/N0 value in dB
EbN0_max_dB         = params.simulation.EbN0_max_dB;                % Maximum Eb/N0 value in dB
EbN0_step_dB        = params.simulation.EbN0_step_dB;               % Step size for Eb/N0 sweep in dB
iterations          = params.simulation.iterations_per_EbN0;        % Iterations for averaging BER at each Eb/N0 point
EbN0_domain_dB      = params.simulation.EbN0_domain_dB;             % Range of Eb/N0 values to simulate (dB)
num_EbN0_points     = length(EbN0_domain_dB);                       % Number of points on the BER curve



%% ========================= Communication Chain =====================
% ---- CFO and Sample Time Offset Parameters ----
Fc = 600e6;                                   % Carrier frequency in Hz
ppm = 1;
delta_cfo = ppm * 1e-6 * Fc;                  % Frequency offset in Hz (0.08 ppm)
phi_0 = 0*randi([0 12]) * (pi/12);            % Phase offset in rad
timing_offset_norm = 0*0.04;                  % Normalized timing offset (% of symbol period)
initial_offset_samples = round(timing_offset_norm * OSF); % Initial offset in samples (1 sample)

% --- Transmitter  ---
bit_tx = randi([0, 1], 1, NumBits).';
symb_tx = mapping(bit_tx, Nbps, ModType);
symb_tx_up = upSampler(symb_tx, OSF).';
h_rrc = rrcFilter(Beta, SymRate, OSF, NumTaps);
signal_tx = applyFilter(symb_tx_up, h_rrc, NumTaps);
signalPower_tx = mean(abs(signal_tx).^2);
Eb = signalPower_tx / BitRate;

time_vector = (0 : length(signal_tx) - 1).' * Ts;
time_vector_symb = (0 : length(symb_tx) - 1).' * Tsymb;

% -- Introduce Noise --
EbN0dB = 1e10;
signal_tx_noisy = addAWGN(signal_tx, Eb, EbN0dB, OSF, SymRate);

% --- Introduce CFO, Phase Offset, and Sample Time Offset ---
signal_tx_distorted = signal_tx_noisy .* exp(1j * ((2 * pi * delta_cfo) * time_vector + phi_0));
signal_tx_distorted = circshift(signal_tx_distorted, initial_offset_samples);

% --- Receiver Chain ---
signal_rx  = applyFilter(signal_tx_distorted, h_rrc, NumTaps);
symb_rx_down = downSampler(signal_rx, OSF);
symb_rx = symb_rx_down .* exp(-1j * (2 * pi * delta_cfo * time_vector_symb));
bit_rx = demapping(symb_rx_down, Nbps, ModType); 


%% =========================  Plots  ================================
bits_to_plot = min(params.timing.NumBits, 100 * Nbps); 
plotConstellation_Tx_Rx(ModOrder, ModType, symb_tx, symb_rx_down);
plotBitstream_Tx_Rx(bit_tx, bit_rx, bits_to_plot);
plotFilterCharacteristics(h_rrc, Beta, Fs, OSF);
plotPSD_Tx_Rx(signal_tx, signal_rx, Fs);
plotBasebandFrequencyResponse(signal_tx, signal_rx, Fs);



% %% ========================= BER (Multiple CFO)  ======================
% ppm_values = [0 2 10 100 200 300 400 500];
% Fc = 600e6;
% phi_0 = 0;
% timing_offset_norm = 0;
% initial_offset_samples = 0;
% 
% % --- Pre-allocate results array ---
% all_ber_data = zeros(num_EbN0_points, length(ppm_values));
% 
% fprintf('\n\n========================================');
% fprintf('\n BER (Multiple CFO) Simulation Setup      ');
% fprintf('\n========================================');
% fprintf('\nModulation: %d-%s, Eb/N0: %.1f-%.1f dB (step = %.1f dB).', ModOrder, upper(ModType), EbN0_min_dB, EbN0_max_dB, EbN0_step_dB);
% fprintf('\nBits per Tx Block: %d, Iterations per Eb/N0: %d', NumBits, iterations);
% fprintf('\nPPM Values for CFO: %s', mat2str(ppm_values));
% fprintf('\n========================================');
% 
% h_rrc = rrcFilter(Beta, SymRate, OSF, NumTaps);
% 
% for idx_ppm = 1:length(ppm_values)
%     current_ppm = ppm_values(idx_ppm);
%     delta_cfo = current_ppm * 1e-6 * Fc;
%     fprintf('\n\nSimulating for CFO = %g ppm, phi_0 = %.f, t_0 = %.fT).', current_ppm, phi_0, timing_offset_norm);
% 
%     ber_data_one_cfo = zeros(num_EbN0_points, 1);
% 
%     for idx_EbN0 = 1:num_EbN0_points
%         EbN0dB = EbN0_domain_dB(idx_EbN0);
% 
%         bit_tx = randi([0, 1], 1, NumBits).';
%         symb_tx = mapping(bit_tx, Nbps, ModType);
%         symb_tx_up = upSampler(symb_tx, OSF).';
%         signal_tx = applyFilter(symb_tx_up, h_rrc, NumTaps);
%         signalPower = mean(abs(signal_tx).^2);
%         Eb = signalPower / BitRate;
% 
%         total_bit_errors_point = 0;
%         total_bits_sim_point = 0;
% 
%         time_vector = (0:length(signal_tx)-1).' * Ts;
%         time_vector_symb = (0:length(symb_tx)-1).' * Tsymb;
% 
%         for iter = 1:iterations
%             signal_tx_noisy = addAWGN(signal_tx, Eb, EbN0dB, OSF, SymRate);
%             signal_tx_distorted = signal_tx_noisy .* exp(1j * (2 * pi * delta_cfo * time_vector + phi_0));
%             signal_tx_distorted = circshift(signal_tx_distorted, initial_offset_samples);
% 
%             signal_rx = applyFilter(signal_tx_distorted, h_rrc, NumTaps);
%             symb_rx_down = downSampler(signal_rx, OSF);
%             symb_rx_down = symb_rx_down .* exp(-1j * (2 * pi * delta_cfo * time_vector_symb));
%             bit_rx = demapping(symb_rx_down, Nbps, ModType);
% 
%             num_errors_iter = sum(bit_tx ~= bit_rx);
%             bits_iter = length(bit_tx);
%             total_bit_errors_point = total_bit_errors_point + num_errors_iter;
%             total_bits_sim_point = total_bits_sim_point + bits_iter;
%         end
%         ber_data_one_cfo(idx_EbN0) = total_bit_errors_point / total_bits_sim_point;
%         fprintf('\n  PPM = %5.1f, Eb/N0 = %5.1f dB : Bits = %8d, Errors = %6d, BER = %.3e', ...
%                   current_ppm, EbN0dB, total_bits_sim_point, total_bit_errors_point, ber_data_one_cfo(idx_EbN0));
%     end
%     all_ber_data(:, idx_ppm) = ber_data_one_cfo;
% end
% fprintf('\n================================================================================\n');
% 
% % Plot
% plotBERCurvesCFO(all_ber_data, params, ppm_values);



%% ========================= BER (Multiple time offset)  ======================
time_offset_norm_values = [0, 0.05]; 
delta_cfo = 0; 
phi_0 = 0;     

all_ber_data_time_offset = zeros(num_EbN0_points, length(time_offset_norm_values));

fprintf('\n\n====================================================');
fprintf('\n BER (Multiple Sample Time Offsets) Simulation Setup  ');
fprintf('\n====================================================');
fprintf('\nModulation: %d-%s, Eb/N0: %.1f-%.1f dB (step = %.1f dB).', ModOrder, upper(ModType), EbN0_min_dB, EbN0_max_dB, EbN0_step_dB);
fprintf('\nBits per Tx Block: %d, Iterations per Eb/N0: %d', NumBits, iterations);
fprintf('\nNormalized Time Offset Values (fraction of Tsymb): %s', mat2str(time_offset_norm_values));
fprintf('\n(Note: Plotted as %% of Tsymb).');
fprintf('\nCFO fixed at: %g Hz, Phase Offset fixed at: %.2f rad', delta_cfo, phi_0);
fprintf('\n====================================================');

for idx_offset = 1:length(time_offset_norm_values)
    current_time_offset_norm = time_offset_norm_values(idx_offset);
    current_initial_offset_samples = round(current_time_offset_norm * OSF); 
    
    fprintf('\n\nSimulating for Time Offset = %.3f Tsymb (samples: %d).', current_time_offset_norm, current_initial_offset_samples);
    
    ber_data_one_offset = zeros(num_EbN0_points, 1);
    for idx_EbN0 = 1:num_EbN0_points
        EbN0dB = EbN0_domain_dB(idx_EbN0);
        
        bit_tx = randi([0, 1], 1, NumBits).'; 
        symb_tx = mapping(bit_tx, Nbps, ModType); 
        symb_tx_up = upSampler(symb_tx, OSF).'; 
        signal_tx = applyFilter(symb_tx_up, h_rrc, NumTaps); 
        
        signalPower = mean(abs(signal_tx).^2);
        Eb = signalPower / BitRate; 
        
        total_bit_errors_point = 0;
        total_bits_sim_point = 0;
        
        time_vector = (0:length(signal_tx)-1).' * Ts; 
        time_vector_symb = (0:length(symb_tx)-1).' * Tsymb; 

        for iter = 1:iterations
            signal_tx_noisy = addAWGN(signal_tx, Eb, EbN0dB, OSF, SymRate);
            
            signal_tx_distorted = signal_tx_noisy .* exp(1j * (2 * pi * delta_cfo * time_vector + phi_0));
            signal_tx_distorted = circshift(signal_tx_distorted, current_initial_offset_samples);
            
            signal_rx = applyFilter(signal_tx_distorted, h_rrc, NumTaps);
            symb_rx = downSampler(signal_rx, OSF);
            
            symb_rx = symb_rx .* exp(-1j * (2 * pi * delta_cfo * time_vector_symb));
            
            bit_rx = demapping(symb_rx, Nbps, ModType);
            
            num_errors_iter = sum(bit_tx ~= bit_rx);
            bits_iter = length(bit_tx); 

            total_bit_errors_point = total_bit_errors_point + num_errors_iter;
            total_bits_sim_point = total_bits_sim_point + bits_iter;
        end
        if total_bits_sim_point > 0
            ber_data_one_offset(idx_EbN0) = total_bit_errors_point / total_bits_sim_point;
        else
            ber_data_one_offset(idx_EbN0) = 1; 
        end
        fprintf('\n t_0 = %.1fTsymb, Eb/N0 = %5.1f dB : Bits = %8d, Errors = %6d, BER = %.3e', ...
                  current_time_offset_norm, EbN0dB, total_bits_sim_point, total_bit_errors_point, ber_data_one_offset(idx_EbN0));
    end
    all_ber_data_time_offset(:, idx_offset) = ber_data_one_offset;
end
fprintf('\n================================================================================\n');

plotBERCurvesTimeOffset(all_ber_data_time_offset, params, time_offset_norm_values);