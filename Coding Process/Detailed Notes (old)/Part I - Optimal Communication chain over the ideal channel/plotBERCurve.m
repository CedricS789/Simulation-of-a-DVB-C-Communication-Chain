function hFig = plotBERCurve(EbN0_domain_dB, BER_simulated, ModType, ModOrder, OSF, Beta)
    %   Generates a Bit Error Rate (BER) curve plot, comparing simulated BER
    %   results against theoretical BER values for Additive White Gaussian Noise (AWGN)
    %   channels over a range of Eb/N0 values.
    %
    %   This function visualizes the performance of a digital communication system
    %   by plotting BER as a function of the signal-to-noise ratio per bit (Eb/N0).
    %   It requires the simulation results (BER vs. Eb/N0) and system parameters
    %   to calculate and plot the corresponding theoretical curve for comparison.
    %
    %   The function uses customizable parameters (defined internally) for plot
    %   appearance (colors, line styles, labels, etc.). It assumes valid inputs;
    %   errors might occur if inputs are incompatible (e.g., invalid ModType/ModOrder
    %   for 'berawgn', or non-positive BER values in BER_simulated).
    %
    %   Inputs:
    %       EbN0_domain_dB - A vector containing the Eb/N0 values (in dB) for which
    %                        the simulation was run. This defines the x-axis of the plot.
    %       BER_simulated  - A vector of the corresponding simulated BER values obtained
    %                        from running the communication system simulation at each
    %                        EbN0_domain_dB point. Must be the same length as EbN0_domain_dB.
    %       ModType        - A string specifying the modulation type used in the simulation.
    %                        Expected values are 'qam' (Quadrature Amplitude Modulation)
    %                        or 'pam' (Pulse Amplitude Modulation). Case-insensitive.
    %       ModOrder       - An integer indicating the modulation order (M), representing
    %                        the number of points in the constellation (e.g., 2 for BPSK,
    %                        4 for QPSK/4-PAM, 16 for 16-QAM/16-PAM, 64 for 64-QAM).
    %       OSF            - The integer Oversampling Factor used in the simulation's
    %                        pulse shaping and matched filtering stages. Included for plot annotation.
    %       Beta           - The Roll-off Factor (0 <= Beta <= 1) of the RRC filter
    %                        used in the simulation. Included for plot annotation.
    %
    %   Output:
    %       hFig           - A handle to the generated figure object. This allows
    %                        further modification of the plot if needed after the
    %                        function call.

        % =====================================================================
        % == Plotting Parameters (Internal Configuration) ==
        % =====================================================================
        % --- Figure Properties ---
        figureName          = 'BER Curve'; % Name appearing in the figure window title bar
        figureNumberTitle   = 'off';       % Controls if 'Figure X:' prefix is shown ('on' or 'off')
        plotTitleFormat     = 'BER Performance for %d-%s (OSF=%d, β=%.2f)'; % Format string for the main plot title

        % --- Axis Labels ---
        xAxisLabel          = 'Eb/N0 (dB)'; % Label for the x-axis
        yAxisLabel          = 'Bit Error Rate (BER)'; % Label for the y-axis

        % --- Simulated Curve Appearance ---
        simLineStyle        = '-';         % Line style (e.g., '-', '--', ':', '-.')
        simMarker           = 'o';         % Marker style (e.g., 'o', 'x', '+', '*', 's')
        simColor            = 'b';         % Color (e.g., 'b', 'r', 'g', [0 0.4470 0.7410])
        simLineWidth        = 1.5;         % Line width
        simMarkerSize       = 6;           % Marker size
        simDisplayName      = 'Simulated'; % Text label for the simulated curve in the legend

        % --- Theoretical Curve Appearance ---
        theoryLineStyle     = '--';        % Line style
        theoryMarker        = 'x';         % Marker style
        theoryColor         = 'r';         % Color (e.g., 'r', [0.8500 0.3250 0.0980])
        theoryLineWidth     = 1.5;         % Line width
        theoryDisplayNameFormat = 'Theoretical %s'; % Format string for the theoretical curve label in the legend (e.g., "Theoretical 16-QAM")

        % --- General Plot Appearance ---
        legendLocation      = 'southwest'; % Position of the legend box (e.g., 'best', 'northeast')
        gridVisible         = 'on';        % Show grid lines ('on' or 'off')
        yAxisMinBERFloor    = 1e-7;        % A minimum floor for the y-axis lower limit to prevent excessively low limits when BER is very small.

        % =====================================================================
        % == Theoretical BER Calculation ==
        % =====================================================================
        % Calculate the theoretical BER values for an AWGN channel using MATLAB's
        % 'berawgn' function. The function requires Eb/N0 (dB), modulation type,
        % modulation order, and optionally coding information (here 'nondiff' for
        % non-differentially coded gray-mapped constellations).

        if strcmpi(ModType, 'qam')
            % Calculate theoretical BER for QAM.
            BER_theoretical = berawgn(EbN0_domain_dB, 'qam', ModOrder, 'nondiff');
            % Base label for the legend entry.
            theory_label_base = sprintf('%d-QAM', ModOrder);
        elseif strcmpi(ModType, 'pam')
             % Calculate theoretical BER for PAM (including BPSK as 2-PAM).
             BER_theoretical = berawgn(EbN0_domain_dB, 'pam', ModOrder, 'nondiff');
             % Special label for BPSK for clarity.
             if ModOrder == 2
                theory_label_base = 'BPSK (2-PAM)';
             else
                theory_label_base = sprintf('%d-PAM', ModOrder);
             end
        else
            % Handle unrecognized modulation types if necessary, though 'berawgn' will error.
            % Assign NaN to prevent plotting errors if ModType was invalid.
            warning('plotBERCurve: Unrecognized ModType "%s". Theoretical curve may be incorrect.', ModType);
            BER_theoretical = nan(size(EbN0_domain_dB));
            theory_label_base = 'Unknown';
        end
        % Format the final display name for the theoretical curve using the base label.
        theoryDisplayName = sprintf(theoryDisplayNameFormat, theory_label_base);

        % =====================================================================
        % == Plotting ==
        % =====================================================================
        % Create a new figure window with the specified properties.
        hFig = figure('Name', figureName, 'NumberTitle', figureNumberTitle);

        % Plot the simulated BER data using a logarithmic scale for the y-axis (BER).
        % Apply the defined appearance settings.
        semilogy(EbN0_domain_dB, BER_simulated,...
                 'LineStyle', simLineStyle,...
                 'Marker', simMarker,...
                 'Color', simColor,...
                 'LineWidth', simLineWidth,...
                 'MarkerSize', simMarkerSize,...
                 'DisplayName', simDisplayName); % Assign display name for legend
        hold on; % Keep the plot active to add the theoretical curve

        % Plot the theoretical BER data on the same axes.
        % 'semilogy' automatically handles NaN values by not plotting them.
        semilogy(EbN0_domain_dB, BER_theoretical,...
                 'LineStyle', theoryLineStyle,...
                 'Marker', theoryMarker,...
                 'Color', theoryColor,...
                 'LineWidth', theoryLineWidth,...
                 'DisplayName', theoryDisplayName); % Assign display name for legend

        % =====================================================================
        % == Plot Formatting ==
        % =====================================================================
        % Configure the appearance of the plot axes, title, legend, etc.
        grid(gridVisible);                     % Turn grid on or off
        xlabel(xAxisLabel);                    % Set x-axis label
        ylabel(yAxisLabel);                    % Set y-axis label
        % Set the main title using the format string and simulation parameters.
        title(sprintf(plotTitleFormat, ModOrder, upper(ModType), OSF, Beta));
        % Display the legend with the specified location.
        legend('show', 'Location', legendLocation);

        % Adjust y-axis limits for better visualization of BER data.
        % Find the minimum non-zero simulated BER to set a reasonable lower limit.
        min_ber_plot = min(BER_simulated(BER_simulated > 0));
        % Set the lower y-axis limit slightly below the minimum simulated BER,
        % but not lower than the defined floor (yAxisMinBERFloor).
        % Use max to ensure the limit isn't excessively low or zero. Use eps if min_ber_plot is empty.
        if isempty(min_ber_plot)
            plot_lower_limit = yAxisMinBERFloor;
        else
            plot_lower_limit = max(min_ber_plot / 10, yAxisMinBERFloor);
        end
        % Set the y-axis limits. The upper limit is set to 1 (max possible BER).
        ylim([plot_lower_limit 1]);

        % Set x-axis limits based on the provided range of Eb/N0 values.
        xlim([min(EbN0_domain_dB) max(EbN0_domain_dB)]);
        hold off; % Release the plot hold
    end