function analyzeAvgWindow(params, pilot_pos, EbN0_domain_dB_plot, num_iterations_for_std_dev, ppm, K_values, N_fixed)

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
    delta_cfo = ppm * 1e-6 * Fc;
    phi_0 = 0;
    g_rrc = rrcFilter(Beta, SymRate, OSF, NumTaps);

    std_cfo_errors_ppm = zeros(length(EbN0_domain_dB_plot), length(K_values));
    std_toa_errors_samples = zeros(length(EbN0_domain_dB_plot), length(K_values));

    for idx_K_loop = 1:length(K_values)
        current_K = K_values(idx_K_loop);
        fprintf('\n N = %d, K = %d', N_fixed, current_K);
        
        for idx_EbN0 = 1:length(EbN0_domain_dB_plot)
            current_EbN0dB = EbN0_domain_dB_plot(idx_EbN0);
            fprintf('\n Eb/N0 = %.1f dB: ', current_EbN0dB);
            temp_cfo_errors_hz = zeros(num_iterations_for_std_dev, 1);
            temp_toa_errors = zeros(num_iterations_for_std_dev, 1);

            for iter_std = 1:num_iterations_for_std_dev
                bit_tx = randi([0, 1], 1, NumBits).';
                symb_tx = mapping(bit_tx, Nbps, ModType);
                pilot = symb_tx(pilot_position : pilot_position + N_fixed - 1);
                signal_tx = upSampler(symb_tx, OSF).';
                signal_tx_filtered = applyFilter(signal_tx, g_rrc, NumTaps);
                signalPower_tx = mean(abs(signal_tx_filtered).^2);
                Eb = signalPower_tx / BitRate;
                time_vector = (0 : length(signal_tx_filtered) - 1).' * Ts;
                signal_tx_noisy = addAWGN(signal_tx_filtered, Eb, current_EbN0dB, OSF, SymRate);
                signal_tx_distorted = signal_tx_noisy .* exp(1j * (2 * pi * delta_cfo * time_vector + phi_0));
                signal_rx = applyFilter(signal_tx_distorted, g_rrc, NumTaps);
                symb_rx_down = downSampler(signal_rx, OSF);
                [toa, delta_cfo_hat] = frameFreqAcquisition(pilot, symb_rx_down, current_K, Tsymb);
                temp_cfo_errors_hz(iter_std) = delta_cfo - delta_cfo_hat;
                temp_toa_errors(iter_std) = pilot_position - toa;
            end
            std_cfo_errors_ppm(idx_EbN0, idx_K_loop) = std(temp_cfo_errors_hz) / (1e-6 * Fc);
            std_toa_errors_samples(idx_EbN0, idx_K_loop) = std(temp_toa_errors);
        end
    end

    title_str_cfo = sprintf('Impact of K parameter on CFO estimation (N=%d, %d-%s)', N_fixed, ModOrder, upper(ModType));
    title_str_toa = sprintf('Impact of K parameter on ToA estimation (N=%d, %d-%s)', N_fixed, ModOrder, upper(ModType));
    figure;
    set(gcf, 'Name', title_str_cfo);
    clf;
    for idx_K_plot = 1:length(K_values)
        displayName = sprintf('N = %d, K = %d', N_fixed, K_values(idx_K_plot));
        plot(EbN0_domain_dB_plot, std_cfo_errors_ppm(:,idx_K_plot), '-o', 'DisplayName', displayName, 'LineWidth', 1, 'MarkerSize', 6);
        hold on;
    end
    grid on;
    xlabel('Ratio $E_b/N_0$ [dB]', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('Std. Dev. of Frequency Error (ppm)', 'Interpreter', 'latex', 'FontSize', 20);
    title(title_str_cfo, 'Interpreter', 'latex', 'FontSize', 30);
    legend('show', 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 30);
    ax = gca;
    ax.XTick = EbN0_domain_dB_plot(1):2:EbN0_domain_dB_plot(end);
    set(gca, 'FontSize', 20);
    hold off;
    figure;
    set(gcf, 'Name', title_str_toa);
    clf;
    for idx_K_plot = 1:length(K_values)
        displayName = sprintf('N = %d, K = %d', N_fixed, K_values(idx_K_plot));
        plot(EbN0_domain_dB_plot, std_toa_errors_samples(:,idx_K_plot), '-o', 'DisplayName', displayName, 'LineWidth', 1, 'MarkerSize', 6);
        hold on;
    end
    grid on;
    xlabel('Ratio $E_b/N_0$ [dB]', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('Std. Dev. of ToA Error (in number of samples)', 'Interpreter', 'latex', 'FontSize', 20);
    title(title_str_toa, 'Interpreter', 'latex', 'FontSize', 30);
    legend('show', 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 30);
    ax = gca;
    ax.XTick = EbN0_domain_dB_plot(1):2:EbN0_domain_dB_plot(end);
    set(gca, 'FontSize', 20);
    hold off;
end