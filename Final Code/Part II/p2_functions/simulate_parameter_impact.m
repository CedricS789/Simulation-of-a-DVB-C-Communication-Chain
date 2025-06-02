function [std_cfo_errors_ppm, std_toa_errors_samples] = simulateParameterImpact(params, analysis_mode, values_to_iterate, fixed_param_value)

    std_cfo_errors_ppm = zeros(length(params.EbN0_domain_dB_plot), length(values_to_iterate));
    std_toa_errors_samples = zeros(length(params.EbN0_domain_dB_plot), length(values_to_iterate));

    for idx_val = 1:length(values_to_iterate)
        val_iterated = values_to_iterate(idx_val);
        current_N_sim = 0; 
        current_K_sim = 0; 

        if strcmpi(analysis_mode, 'N_variant')
            current_N_sim = val_iterated;
            current_K_sim = fixed_param_value;
            fprintf('\n Simulating for N = %d, K = %d\n', current_N_sim, current_K_sim);
            NumBits_sim = (params.pilot_position + current_N_sim + current_K_sim + 150) * params.Nbps;
        elseif strcmpi(analysis_mode, 'K_variant')
            current_N_sim = fixed_param_value;
            current_K_sim = val_iterated;
            fprintf('\n Simulating for N = %d, K = %d\n', current_N_sim, current_K_sim);
            NumBits_sim = (params.pilot_position + current_N_sim + current_K_sim + 150) * params.Nbps;
        else
            error('Invalid analysis_mode. Must be "N_variant" or "K_variant".');
        end

        for idx_EbN0 = 1:length(params.EbN0_domain_dB_plot)
            current_EbN0dB = params.EbN0_domain_dB_plot(idx_EbN0);
            fprintf('   Eb/N0 = %.1f dB: ', current_EbN0dB);
            temp_cfo_errors_hz = zeros(params.num_iterations_for_std_dev, 1);
            temp_toa_errors = zeros(params.num_iterations_for_std_dev, 1);

            for iter_std = 1:params.num_iterations_for_std_dev
                if mod(iter_std, round(params.num_iterations_for_std_dev/10)) == 0 && iter_std ~= params.num_iterations_for_std_dev && params.num_iterations_for_std_dev >=10, fprintf('.'); end
                
                bit_tx = randi([0, 1], 1, NumBits_sim).';
                symb_tx = mapping(bit_tx, params.Nbps, params.ModType);
                
                if params.pilot_position + current_N_sim - 1 > length(symb_tx)
                     error('Pilot construction error: symb_tx too short for pilot_seq. NumBits_sim might be too small.');
                end
                pilot_seq = symb_tx(params.pilot_position : params.pilot_position + current_N_sim - 1);
                
                signal_tx = upSampler(symb_tx, params.OSF).';
                signal_tx_filtered = applyFilter(signal_tx, params.g_rrc, params.NumTaps);
                
                signalPower_tx = mean(abs(signal_tx_filtered).^2);
                Eb = signalPower_tx / params.BitRate;
                
                time_vector = (0 : length(signal_tx_filtered) - 1).' * params.Ts;
                
                signal_tx_noisy = addAWGN(signal_tx_filtered, Eb, current_EbN0dB, params.OSF, params.SymRate);
                
                signal_tx_distorted = signal_tx_noisy .* exp(1j * (2 * pi * params.delta_cfo * time_vector + params.phi_0));
                
                signal_rx = applyFilter(signal_tx_distorted, params.g_rrc, params.NumTaps);
                
                symb_rx_down = downSampler(signal_rx, params.OSF);
                
                [toa, delta_cfo_hat] = frameFreqAcquisition(pilot_seq, symb_rx_down, current_K_sim, params.Tsymb);
                
                temp_cfo_errors_hz(iter_std) = params.delta_cfo - delta_cfo_hat;
                temp_toa_errors(iter_std) = params.pilot_position - toa;
            end
            fprintf(' done.\n');
            std_cfo_errors_ppm(idx_EbN0, idx_val) = std(temp_cfo_errors_hz) / (1e-6 * params.Fc);
            std_toa_errors_samples(idx_EbN0, idx_val) = std(temp_toa_errors);
        end
    end
end