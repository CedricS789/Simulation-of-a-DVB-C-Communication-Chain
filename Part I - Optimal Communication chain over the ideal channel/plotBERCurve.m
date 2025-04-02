function hFig = plotBERCurve(EbN0_domain_dB, BER_simulated, ModType, ModOrder, OSF, Beta)
    %   Generates a Bit Error Rate (BER) curve plot, comparing simulated BER
    %   results against theoretical BER values for AWGN channels.
    %
    %   Inputs:
    %       EbN0_domain_dB - Vector of Eb/N0 values (dB) (x-axis).
    %       BER_simulated  - Vector of corresponding simulated BER values.
    %       ModType        - Modulation type string ('qam' or 'pam').
    %       ModOrder       - Modulation order (M).
    %       OSF            - Oversampling Factor (for annotation).
    %       Beta           - RRC Roll-off Factor (for annotation).
    %
    %   Output:
    %       hFig           - Handle to the generated figure object.

        % =====================================================================
        % == Plotting Parameters (Internal Configuration) ==
        % =====================================================================
        figureName          = 'BER Curve';
        figureNumberTitle   = 'off';
        plotTitleFormat     = 'BER Performance for %d-%s (OSF=%d, β=%.2f)';

        xAxisLabel          = 'Eb/N0 (dB)';
        yAxisLabel          = 'Bit Error Rate (BER)';

        simLineStyle        = '-';
        simMarker           = 'o';
        simColor            = 'b';
        simLineWidth        = 1.5;
        simMarkerSize       = 6;
        simDisplayName      = 'Simulated';

        theoryLineStyle     = '--';
        theoryMarker        = 'x';
        theoryColor         = 'r';
        theoryLineWidth     = 1.5;
        theoryDisplayNameFormat = 'Theoretical %s'; % e.g., "Theoretical 16-QAM"

        legendLocation      = 'southwest';
        gridVisible         = 'on';
        yAxisMinBERFloor    = 1e-7;        % Minimum floor for y-axis lower limit.

        % =====================================================================
        % == Theoretical BER Calculation ==
        % =====================================================================
        % Use MATLAB's 'berawgn' for theoretical AWGN performance.
        if strcmpi(ModType, 'qam')
            BER_theoretical = berawgn(EbN0_domain_dB, 'qam', ModOrder, 'nondiff');
            theory_label_base = sprintf('%d-QAM', ModOrder);
        elseif strcmpi(ModType, 'pam')
             BER_theoretical = berawgn(EbN0_domain_dB, 'pam', ModOrder, 'nondiff');
             if ModOrder == 2 % Special label for BPSK
                theory_label_base = 'BPSK (2-PAM)';
             else
                theory_label_base = sprintf('%d-PAM', ModOrder);
             end
        else
            warning('plotBERCurve: Unrecognized ModType "%s". Theoretical curve may be incorrect.', ModType);
            BER_theoretical = nan(size(EbN0_domain_dB));
            theory_label_base = 'Unknown';
        end
        theoryDisplayName = sprintf(theoryDisplayNameFormat, theory_label_base);

        % =====================================================================
        % == Plotting ==
        % =====================================================================
        hFig = figure('Name', figureName, 'NumberTitle', figureNumberTitle);

        semilogy(EbN0_domain_dB, BER_simulated,...
                 'LineStyle', simLineStyle,...
                 'Marker', simMarker,...
                 'Color', simColor,...
                 'LineWidth', simLineWidth,...
                 'MarkerSize', simMarkerSize,...
                 'DisplayName', simDisplayName);
        hold on;

        semilogy(EbN0_domain_dB, BER_theoretical,...
                 'LineStyle', theoryLineStyle,...
                 'Marker', theoryMarker,...
                 'Color', theoryColor,...
                 'LineWidth', theoryLineWidth,...
                 'DisplayName', theoryDisplayName);

        % =====================================================================
        % == Plot Formatting ==
        % =====================================================================
        grid(gridVisible);
        xlabel(xAxisLabel);
        ylabel(yAxisLabel);
        title(sprintf(plotTitleFormat, ModOrder, upper(ModType), OSF, Beta));
        legend('show', 'Location', legendLocation);

        % Adjust y-axis limits for better visualization
        min_ber_plot = min(BER_simulated(BER_simulated > 0));
        if isempty(min_ber_plot)
            plot_lower_limit = yAxisMinBERFloor;
        else
             % Set lower limit slightly below min simulated BER, but not below floor.
            plot_lower_limit = max(min_ber_plot / 10, yAxisMinBERFloor);
        end
        ylim([plot_lower_limit 1]);

        xlim([min(EbN0_domain_dB) max(EbN0_domain_dB)]);
        hold off;
end