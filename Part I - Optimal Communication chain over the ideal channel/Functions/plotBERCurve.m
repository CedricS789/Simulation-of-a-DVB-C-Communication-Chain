function hFig = plotBERCurve(EbN0_domain_dB, BER_simulated, ModType, ModOrder, OSF, Beta)
    %   Generates a Bit Error Rate (BER) curve plot, comparing simulated BER
    %   results against theoretical BER values for AWGN channels.
    %
    %   Inputs:
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

        % --- SIMULATION ---
        simLineStyle        = 'none';
        simMarker           = 'o';
        simColor            = 'b';
        simLineWidth        = 1.5;
        simMarkerSize       = 6;
        simDisplayName      = 'Simulated';

        % --- THEORY ---
        theoryLineStyle     = '-';
        theoryMarker        = 'none';
        theoryColor         = 'r';
        theoryLineWidth     = 1;
        theoryDisplayNameFormat = 'Theoretical %s';
        nPointsTheory       = 300;                      % Number of points for smooth  curve

        legendLocation      = 'southwest';
        gridVisible         = 'on';
        yAxisMinBERFloor    = 1e-4;                     % Reasonable floor for BER plots

        % =====================================================================
        % == Theoretical BER Calculation for SMOOTH CURVE ==
        % =====================================================================
        % Create a denser Eb/N0 vector just for the theoretical plot
        % Ensure EbN0_domain_dB is not empty and has at least two points for linspace
        EbN0_theory_dB = linspace(min(EbN0_domain_dB), max(EbN0_domain_dB), nPointsTheory);

        % Calculate theoretical BER over the DENSE vector if possible
        if ~isempty(EbN0_theory_dB)
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
            else
                warning('plotBERCurve: Unrecognized ModType "%s". Theoretical curve may be incorrect.', ModType);
                BER_theoretical_smooth = nan(size(EbN0_theory_dB)); % Use the dense vector size
                theory_label_base = 'Unknown';
            end
            theoryDisplayName = sprintf(theoryDisplayNameFormat, theory_label_base);
        else
            % Handle case where theoretical curve can't be generated
            BER_theoretical_smooth = [];
            theoryDisplayName = 'Theoretical (N/A)';
        end

        % =====================================================================
        % == Plotting ==
        % =====================================================================
        hFig = figure('Name', figureName, 'NumberTitle', figureNumberTitle);

        % Plot Theoretical Data (Smooth Line Only - using dense EbN0) if available
        if ~isempty(BER_theoretical_smooth)
            semilogy(EbN0_theory_dB, BER_theoretical_smooth,... % Use dense vectors
                     'LineStyle', theoryLineStyle,...
                     'Marker', theoryMarker,...
                     'Color', theoryColor,...
                     'LineWidth', theoryLineWidth,...
                     'DisplayName', theoryDisplayName);
        end
        hold on;
        
        % Plot Simulated Data (Dots Only - using original sparse EbN0)
        semilogy(EbN0_domain_dB, BER_simulated,...
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
        valid_ber_sim = BER_simulated(BER_simulated > 0 & ~isnan(BER_simulated));
        if isempty(valid_ber_sim)
             min_ber_plot = yAxisMinBERFloor; % Default if no valid points
        else
             min_ber_plot = min(valid_ber_sim);
        end

        if isempty(min_ber_plot) || min_ber_plot <= 0 % Check against zero or less
            plot_lower_limit = yAxisMinBERFloor;
        else
             % Set lower limit slightly below min simulated BER, but not below floor.
            plot_lower_limit = max(min_ber_plot / 10, yAxisMinBERFloor);
        end
        ylim([plot_lower_limit 1]); % Upper limit often 1 for BER

        % Set x-axis limits based on the original simulation Eb/N0 range
        if ~isempty(EbN0_domain_dB)
           xlim([min(EbN0_domain_dB) max(EbN0_domain_dB)]);
        end
        hold off;
end