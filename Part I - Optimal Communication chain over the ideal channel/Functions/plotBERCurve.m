function hFig = plotBERCurve(ber_averages_datas, params)
    %   Generates a simplified Bit Error Rate (BER) curve plot, comparing
    %   simulated BER results against theoretical BER values for AWGN channels.
    %   Assumes valid, non-empty inputs for plotting.
    %
    %   Inputs - params
    %       EbN0_domain_dB - Vector of Eb/N0 values (dB) for SIMULATION (x-axis).
    %       BER_simulated  - Vector of corresponding simulated BER values.
    %       ModType        - Modulation type string ('qam' or 'pam').
    %       ModOrder       - Modulation order (M).
    %       OSF            - Oversampling Factor (for annotation).
    %       Beta           - RRC Roll-off Factor (for annotation).
    %
    %   Output:
    %       hFig           - Handle to the generated figure object.

        
        % =====================================================================
        % == Plotting Parameters ==
        % =====================================================================
        figureName          = 'BER Curve';
        figureNumberTitle   = 'off';
        plotTitleFormat     = 'BER Performance for %d-%s (OSF=%d, β=%.2f)';
        xAxisLabel          = 'Eb/N0 (dB)';
        yAxisLabel          = 'Bit Error Rate (BER)';

        % --- SIMULATION STYLING ---
        simLineStyle        = 'none';
        simMarker           = 'o';
        simColor            = 'b';
        simLineWidth        = 1.5;
        simMarkerSize       = 6;
        simDisplayName      = 'Simulated';

        % --- THEORY STYLING ---
        theoryLineStyle     = '-';
        theoryMarker        = 'none';
        theoryColor         = 'r';
        theoryLineWidth     = 1;
        theoryDisplayNameFormat = 'Theoretical %s';
        nPointsTheory       = 300;                      % Number of points for smooth curve

        legendLocation      = 'southwest';
        gridVisible         = 'on';
        yAxisMinBERFloor    = 1e-4;                     % Reasonable floor for BER plots
        
        % =====================================================================
        % == Input Parameters ==
        % =====================================================================
        EbN0_domain_dB = params.simulation.EbN0_domain_dB;
        ModType        = params.modulation.ModulationType;
        ModOrder       = params.modulation.ModulationOrder;
        OSF            = params.sampling.OversamplingFactor;
        Beta           = params.filter.RolloffFactor;

        % =====================================================================
        % == Theoretical BER Calculation for SMOOTH CURVE ==
        % =====================================================================
        % Create a denser Eb/N0 vector for the theoretical plot
        % Assumes EbN0_domain_dB has at least two points
        EbN0_theory_dB = linspace(min(EbN0_domain_dB), max(EbN0_domain_dB), nPointsTheory);

        % Calculate theoretical BER over the DENSE vector
        % Assumes ModType is 'qam' or 'pam'
        if strcmpi(ModType, 'qam')
            BER_theoretical_smooth = berawgn(EbN0_theory_dB, 'qam', ModOrder, 'nondiff');
            theory_label_base = sprintf('%d-QAM', ModOrder);
        elseif strcmpi(ModType, 'pam')
             BER_theoretical_smooth = berawgn(EbN0_theory_dB, 'pam', ModOrder, 'nondiff');
             if ModOrder == 2 % Special label for BPSK
                theory_label_base = 'BPSK (2-PAM)';
             else
                theory_label_base = sprintf('%d-PAM', ModOrder);
             end
        end
        theoryDisplayName = sprintf(theoryDisplayNameFormat, theory_label_base);

        % =====================================================================
        % == Plotting ==
        % =====================================================================
        hFig = figure('Name', figureName, 'NumberTitle', figureNumberTitle);

        % Plot Theoretical Data
        semilogy(EbN0_theory_dB, BER_theoretical_smooth,... % Use dense vectors
                 'LineStyle', theoryLineStyle,...
                 'Marker', theoryMarker,...
                 'Color', theoryColor,...
                 'LineWidth', theoryLineWidth,...
                 'DisplayName', theoryDisplayName);

        hold on;

        % Plot Simulated Data
        semilogy(EbN0_domain_dB, ber_averages_datas,...
                 'LineStyle', simLineStyle,...
                 'Marker', simMarker,...
                 'Color', simColor,...
                 'LineWidth', simLineWidth,...
                 'MarkerSize', simMarkerSize,...
                 'MarkerFaceColor', simColor,...
                 'DisplayName', simDisplayName);

        % =====================================================================
        % == Plot Formatting ==
        % =====================================================================
        grid(gridVisible);
        xlabel(xAxisLabel);
        ylabel(yAxisLabel);
        title(sprintf(plotTitleFormat, ModOrder, upper(ModType), OSF, Beta));
        legend('show', 'Location', legendLocation);

        % Adjust y-axis limits based on SIMULATED data and floor
        valid_ber_sim = ber_averages_datas(ber_averages_datas > 0 & ~isnan(ber_averages_datas));
        if isempty(valid_ber_sim)
             min_ber_plot = yAxisMinBERFloor; % Default if no valid simulated points > 0
        else
             min_ber_plot = min(valid_ber_sim);
        end
        % Set lower limit slightly below min simulated BER, but not below floor.
        plot_lower_limit = max(min_ber_plot / 10, yAxisMinBERFloor);
        ylim([plot_lower_limit 1]); % Upper limit often 1 for BER

        % Set x-axis limits based on the original simulation Eb/N0 range
        % Assumes EbN0_domain_dB is not empty
        xlim([min(EbN0_domain_dB) max(EbN0_domain_dB)]);

        hold off;
end