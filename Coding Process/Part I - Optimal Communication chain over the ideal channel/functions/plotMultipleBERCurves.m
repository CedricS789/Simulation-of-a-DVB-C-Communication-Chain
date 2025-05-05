function hFig = plotMultipleBERCurves(nbps_values, ber_results_cell, params_cell)
    hFig = figure();                % Create a new figure window
    axBER = axes('Parent', hFig);   % Create axes for BER plot
    hold(axBER, 'on');              % Hold on to the current axes
    set(axBER, 'YScale', 'log');    % Set Y-axis to logarithmic scale

    % --- Dynamic Title ---
    title_str = sprintf('BER Comparison ($\\beta=%.2f$, OSF=%d, Filter Taps=%d)', ...
                        params_cell{1}.filter.RolloffFactor, ...
                        params_cell{1}.sampling.OversamplingFactor, ...
                        params_cell{1}.filter.NumFilterTaps);

    % --- Plotting Loop ---
    for i = 1:length(nbps_values)
        % --- Get Data ---
        current_params  = params_cell{i};                               % Current params
        EbN0_domain_dB  = current_params.simulation.EbN0_domain_dB;     % Eb/N0 domain
        ModType         = current_params.modulation.ModulationType;     % Modulation type
        ModOrder        = current_params.modulation.ModulationOrder;    % Modulation order

        % ---- Calculate Theoretical BER using berawgn. However, beragwn depends on the modulation type ----
        % ---- So we need to check the modulation type and calculate the theoretical BER accordingly ----
        EbN0_theory_dB = linspace(min(EbN0_domain_dB), max(EbN0_domain_dB), 300);   % vector representing Eb/N0 range in dB with 300 points
        if strcmpi(ModType, 'qam')                                                  % if modulation type is QAM
            BER_theory = berawgn(EbN0_theory_dB, 'qam', ModOrder, 'nondiff');       % calculate theoretical BER for QAM
            label_base = sprintf('%d-QAM (Nbps=%d)', ModOrder, nbps_values(i));     % label for QAM
        elseif strcmpi(ModType, 'pam')                                              % if modulation type is PAM
             BER_theory = berawgn(EbN0_theory_dB, 'pam', ModOrder, 'nondiff');      % calculate theoretical BER for PAM
             if ModOrder == 2                                                       % if modulation order is 2 (BPSK), special case of PAM
                 label_base = sprintf('2-PAM/BPSK (Nbps=%d)', nbps_values(i));   % label for BPSK/PAM
             else
                 label_base = sprintf('%d-PAM (Nbps=%d)', ModOrder, nbps_values(i)); % label for PAM
             end
        end

        % --- Plot ---
        hTheo = semilogy(axBER, EbN0_theory_dB, BER_theory, 'HandleVisibility', 'off', 'LineWidth', 2);     % Theory line
        hSim = semilogy(axBER, EbN0_domain_dB, ber_results_cell{i}, ...                     % Simulation markers
         'LineStyle', 'none', 'Marker', 'o', 'DisplayName', label_base, ...
         'MarkerSize', 7);   

         % --- Match color ---
        hSim.Color = hTheo.Color;
        hSim.MarkerFaceColor = hTheo.Color; % Fill the markers with the same color as the line
    end

    % --- Final Formatting ---
    xlabel(axBER, '$E_b/N_0$ (dB)', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel(axBER, 'Bit Error Rate (BER)', 'FontSize', 14, 'Interpreter', 'latex');
    title(axBER, title_str, 'FontSize', 16, 'Interpreter', 'latex');
    legend(axBER, 'show', 'Location', 'northeast', 'FontSize', 16, 'Interpreter', 'latex')
    xlim(axBER, [min(params_cell{1}.simulation.EbN0_domain_dB) max(params_cell{1}.simulation.EbN0_domain_dB)]);  % set x-axis limits as the same as the first BER data
    ylim(axBER, [1e-3 1]);
    grid(axBER, 'on');     % Turn on grid
    axBER.FontSize = 12;   % Axes tick labels size
    hold(axBER, 'off');    % Release hold on axes

end