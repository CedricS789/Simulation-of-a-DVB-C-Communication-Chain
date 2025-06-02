function cfoRobustness(params, pilot_pos, EbN0_domain_dB_plot, num_iterations_for_std_dev, ppm_values, N_fixed, K_fixed, timing_offset_percent)

    Nbps = params.modulation.Nbps;
    ModType = params.modulation.ModulationType;
    ModOrder = params.modulation.ModulationOrder;
    OSF = params.sampling.OversamplingFactor;
    NumTaps = params.filter.NumFilterTaps;
    BitRate = params.timing.BitRate;
    Ts = params.sampling.SamplePeriod;
    SymRate = params.timing.SymbolRate;
    Tsymb = params.timing.SymbolPeriod;
    Beta = params.filter.RolloffFactor;
    NumBits = params.timing.NumBits;

    pilot_position = pilot_pos;
    Fc = 600e6;
    phi_0 = 0;

    g_rrc = rrcFilter(Beta, SymRate, OSF, NumTaps);
    initial_offset_samples = round(timing_offset_percent * OSF);

    std_cfo_errors_ppm = zeros(length(EbN0_domain_dB_plot), length(ppm_values));
    std_toa_errors_samples = zeros(length(EbN0_domain_dB_plot), length(ppm_values));

    for idx_ppm_loop = 1:length(ppm_values)
        current_ppm = ppm_values(idx_ppm_loop);
        delta_cfo = current_ppm * 1e-6 * Fc;
        fprintf('\n%d-%s: N=%d, K=%d, CFO=%dppm, t0=%.2fTsymb', ModOrder, upper(ModType), N_fixed, K_fixed, current_ppm, timing_offset_percent);

        for idx_EbN0 = 1:length(EbN0_domain_dB_plot)
            current_EbN0dB = EbN0_domain_dB_plot(idx_EbN0);
            fprintf('\n Eb/N0 = %.1f dB: ', current_EbN0dB);
            temp_cfo_errors_hz = zeros(num_iterations_for_std_dev, 1);
            temp_toa_errors = zeros(num_iterations_for_std_dev, 1);

            for iter_std = 1:num_iterations_for_std_dev
                bit_tx = randi([0, 1], 1, NumBits).';
                symb_tx = mapping(bit_tx, Nbps, ModType);
                pilot = symb_tx(pilot_position : pilot_position + N_fixed - 1);
                symb_tx_up = upSampler(symb_tx, OSF).';
                signal_tx_filtered = applyFilter(symb_tx_up, g_rrc, NumTaps);
                signalPower_tx = mean(abs(signal_tx_filtered).^2);
                Eb = signalPower_tx / BitRate;
                time_vector = (0 : length(signal_tx_filtered) - 1).' * Ts;
                signal_tx_noisy = addAWGN(signal_tx_filtered, Eb, current_EbN0dB, OSF, SymRate);
                signal_tx_cfo_distorted = signal_tx_noisy .* exp(1j * (2 * pi * delta_cfo * time_vector + phi_0));
                signal_tx_time_distorted = circshift(signal_tx_cfo_distorted, initial_offset_samples);
                signal_rx_matched_filtered = applyFilter(signal_tx_time_distorted, g_rrc, NumTaps);
                symb_rx_down = downSampler(signal_rx_matched_filtered, OSF);
                [toa_hat, delta_cfo_hat] = frameFreqAcquisition(pilot, symb_rx_down, K_fixed, Tsymb);
                temp_cfo_errors_hz(iter_std) = delta_cfo - delta_cfo_hat;
                temp_toa_errors(iter_std) = pilot_position - toa_hat;
            end
            std_cfo_errors_ppm(idx_EbN0, idx_ppm_loop) = std(temp_cfo_errors_hz) / (1e-6 * Fc);
            std_toa_errors_samples(idx_EbN0, idx_ppm_loop) = std(temp_toa_errors);
        end
         fprintf('\n');
    end

    % Plotting
    title_str_toa = sprintf('Robustness to CFO on ToA estimate (%d-%s, N=%d, K=%d)', ModOrder, upper(ModType), N_fixed, K_fixed);
    figure;
    set(gcf, 'Name', title_str_toa);
    clf;
    for idx_ppm_plot = 1:length(ppm_values)
        % Using 't_0' for timing offset in legend, consistent with images like Fig. 17, 18
        displayName = sprintf('CFO = %d ppm, $t_0 = %.2fT_{symb}$', ppm_values(idx_ppm_plot), timing_offset_percent);
        plot(EbN0_domain_dB_plot, std_toa_errors_samples(:,idx_ppm_plot), '-o', 'DisplayName', displayName, 'LineWidth', 1, 'MarkerSize', 6);
        hold on;
    end
    grid on;
    xlabel('Ratio $E_b/N_0$ [dB]', 'Interpreter', 'latex', 'FontSize', 30);
    ylabel('Std. Dev. of ToA Error (samples)', 'Interpreter', 'latex', 'FontSize', 30);
    title(title_str_toa, 'Interpreter', 'latex', 'FontSize', 30);
    legend('show', 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 30);
    ax = gca;
    ax.GridAlpha = 0.5;
    ax.XTick = EbN0_domain_dB_plot(1):2:EbN0_domain_dB_plot(end);
    set(gca, 'FontSize', 30);
    hold off;
    box on;

    title_str_cfo = sprintf('Robustness to CFO on CFO estimate (%d-%s, N=%d, K=%d)', ModOrder, upper(ModType), N_fixed, K_fixed);
    figure;
    set(gcf, 'Name', title_str_cfo);
    clf;
    for idx_ppm_plot = 1:length(ppm_values)
        displayName = sprintf('CFO = %d ppm, $t_0 = %.2fT_{symb}$', ppm_values(idx_ppm_plot), timing_offset_percent);
        plot(EbN0_domain_dB_plot, std_cfo_errors_ppm(:,idx_ppm_plot), '-o', 'DisplayName', displayName, 'LineWidth', 1, 'MarkerSize', 6);
        hold on;
    end
    grid on;
    xlabel('Ratio $E_b/N_0$ [dB]', 'Interpreter', 'latex', 'FontSize', 30);
    ylabel('Std. Dev. of Frequency Error (ppm)', 'Interpreter', 'latex', 'FontSize', 30);
    title(title_str_cfo, 'Interpreter', 'latex', 'FontSize', 30);
    legend('show', 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 30);
    ax = gca;
    ax.GridAlpha = 0.5;
    ax.XTick = EbN0_domain_dB_plot(1):2:EbN0_domain_dB_plot(end);
    set(gca, 'FontSize', 30);
    hold off;
    box on;
end