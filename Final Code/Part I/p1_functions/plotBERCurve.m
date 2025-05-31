function hFig = plotBERCurve(ber_data, params)
    % simulated BER results against theoretical BER values for AWGN channels.

    % Plotting Parameters
    figureName = 'BER Curve';
    figureNumberTitle = 'off';
    plotTitleFormat = 'BER Performance for %d-%s (OSF=%d)';
    xAxisLabel = 'Eb/N0 (dB)';
    yAxisLabel = 'Bit Error Rate (BER)';

    % Input Parameters
    EbN0_domain_dB = params.simulation.EbN0_domain_dB;
    ModType = params.modulation.ModulationType;
    ModOrder = params.modulation.ModulationOrder;
    OSF = params.sampling.OversamplingFactor;

    % Theoretical BER Calculation
    if strcmpi(ModType, 'qam')
        BER_theoretical_smooth = berawgn(EbN0_domain_dB, 'qam', ModOrder, 'nondiff');
        theory_label_base = sprintf('%d-QAM', ModOrder);
    elseif strcmpi(ModType, 'pam')
        BER_theoretical_smooth = berawgn(EbN0_domain_dB, 'pam', ModOrder, 'nondiff');
        if ModOrder == 2
            theory_label_base = 'BPSK (2-PAM)';
        else
            theory_label_base = sprintf('%d-PAM', ModOrder);
        end
    end
    theoryDisplayName = sprintf('Theoretical %s', theory_label_base);

    % Plotting
    hFig = figure('Name', figureName, 'NumberTitle', figureNumberTitle);
    
    semilogy(EbN0_domain_dB, BER_theoretical_smooth, 'LineStyle', '-', 'Marker', 'o', 'DisplayName', theoryDisplayName);
    hold on;
    semilogy(EbN0_domain_dB, ber_data, 'LineStyle', '-', 'Marker', '^', 'DisplayName', 'Simulated');

    % Plot Formatting
    grid('on');
    xlabel(xAxisLabel);
    ylabel(yAxisLabel);
    title(sprintf(plotTitleFormat, ModOrder, upper(ModType), OSF));
    legend('show', 'Location', 'northeast');

    ylim([1e-5 1]);
    xlim([min(EbN0_domain_dB) max(EbN0_domain_dB)]);
    hold off;
end