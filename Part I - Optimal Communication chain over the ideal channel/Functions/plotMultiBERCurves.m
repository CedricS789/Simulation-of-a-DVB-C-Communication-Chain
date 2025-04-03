function hFig = plotMultiBERCurves(results_cell_array)
    %PLOTMULTIBERCURVES Plots multiple BER curves on a single figure.
    %   Generates a Bit Error Rate (BER) curve plot, comparing simulated BER
    %   results against theoretical BER values for AWGN channels for multiple
    %   modulation schemes provided in the input cell array.
    %
    %   Inputs:
    %       results_cell_array - Cell array where each cell contains a struct
    %                            with fields:
    %                               - EbN0_dB     : Vector of Eb/N0 values (dB) for SIMULATION.
    %                               - BER_sim     : Vector of corresponding simulated BER values.
    %                               - params      : The parameter struct used for that simulation
    %                                               (needed for ModType, ModOrder, OSF, Beta).
    %
    %   Output:
    %       hFig           - Handle to the generated figure object.

    % =====================================================================
    % == Plotting Parameters ==
    % =====================================================================
    figureName          = 'BER Curve Comparison';
    figureNumberTitle   = 'off';
    plotTitle           = 'BER Performance Comparison'; % Title is now generic
    xAxisLabel          = 'Eb/N0 (dB)';
    yAxisLabel          = 'Bit Error Rate (BER)';

    % --- Define Colors and Markers to cycle through ---
    simColors   = {'b', 'g', 'm', 'c', 'k'};
    simMarkers  = {'o', 's', '^', 'd', 'p'}; % Circle, Square, Triangle, Diamond, Pentagon
    theoryColors = {'r', [0 0.6 0], [0.6 0 0.6], [0 0.7 0.7], [0.3 0.3 0.3]}; % Red, Dark Green, Purple, Teal, Gray
    theoryLineStyles = {'-', '--', ':', '-.', '-'}; % Solid, Dashed, Dotted, Dash-Dot, Solid

    simLineStyle        = 'none'; % Keep simulated as markers only
    simLineWidth        = 1.5;
    simMarkerSize       = 6;

    theoryMarker        = 'none'; % Keep theory as lines only
    theoryLineWidth     = 1;
    nPointsTheory       = 300;                      % Number of points for smooth curve

    legendLocation      = 'southwest';
    gridVisible         = 'on';
    yAxisMinBERFloor    = 1e-6;                     % Lower floor for better comparison range

    num_results = length(results_cell_array);
    if num_results == 0
        warning('plotMultiBERCurves: No results provided to plot.');
        hFig = figure; % Return an empty figure handle
        return;
    end

    % =====================================================================
    % == Plotting ==
    % =====================================================================
    hFig = figure('Name', figureName, 'NumberTitle', figureNumberTitle);
    hold on; % Hold on to add all curves

    legendEntries = {}; % Cell array to store legend text
    min_sim_ber_overall = 1; % Initialize to max possible BER
    max_ebn0_overall = -inf; % Initialize for x-axis limit calculation
    min_ebn0_overall = inf;  % Initialize for x-axis limit calculation

    for i = 1:num_results
        % --- Extract data for this curve ---
        current_result = results_cell_array{i};
        EbN0_domain_dB = current_result.EbN0_dB;
        BER_simulated  = current_result.BER_sim;
        params         = current_result.params;
        ModType        = params.modulation.ModulationType;
        ModOrder       = params.modulation.ModulationOrder;
        % OSF            = params.sampling.OversamplingFactor; % For annotation if needed
        % Beta           = params.filter.RolloffFactor;      % For annotation if needed

        % --- Assign colors/markers based on index ---
        color_idx = mod(i-1, length(simColors)) + 1;
        marker_idx = mod(i-1, length(simMarkers)) + 1;
        theory_color_idx = mod(i-1, length(theoryColors)) + 1;
        theory_style_idx = mod(i-1, length(theoryLineStyles)) + 1;

        % --- Format modulation label ---
         if strcmpi(ModType, 'qam')
             mod_label_base = sprintf('%d-QAM', ModOrder);
         elseif strcmpi(ModType, 'pam')
             if ModOrder == 2 % Special label for BPSK
                 mod_label_base = 'BPSK (2-PAM)';
             else
                 mod_label_base = sprintf('%d-PAM', ModOrder);
             end
         else
             mod_label_base = sprintf('%d-Unknown (%s)', ModOrder, upper(ModType));
         end

        simDisplayName = sprintf('Sim: %s', mod_label_base);
        theoryDisplayName = sprintf('Theory: %s', mod_label_base);

        % --- Plot Simulated Data ---
        semilogy(EbN0_domain_dB, BER_simulated,...
                 'LineStyle', simLineStyle,...
                 'Marker', simMarkers{marker_idx},...
                 'Color', simColors{color_idx},...
                 'LineWidth', simLineWidth,...
                 'MarkerSize', simMarkerSize,...
                 'MarkerFaceColor', simColors{color_idx},...
                 'DisplayName', simDisplayName);
        legendEntries{end+1} = simDisplayName; % Add to legend list

        % Update overall min BER and EbN0 range
        valid_ber_sim = BER_simulated(BER_simulated > 0 & ~isnan(BER_simulated));
        if ~isempty(valid_ber_sim)
            min_sim_ber_overall = min(min_sim_ber_overall, min(valid_ber_sim));
        end
        if ~isempty(EbN0_domain_dB)
            min_ebn0_overall = min(min_ebn0_overall, min(EbN0_domain_dB));
            max_ebn0_overall = max(max_ebn0_overall, max(EbN0_domain_dB));
        end

        % --- Calculate and Plot Theoretical Data ---
        % Create a dense Eb/N0 vector spanning the range relevant FOR THIS CURVE
        % or use the overall range found so far. Using overall is better for consistent theory lines.
        if isinf(min_ebn0_overall) || isinf(max_ebn0_overall) % Handle case of first run or no valid EbN0
           ebn0_theory_start = 0; ebn0_theory_end = 1; % Default placeholder
        else
           ebn0_theory_start = min_ebn0_overall;
           ebn0_theory_end   = max_ebn0_overall;
        end
        EbN0_theory_dB = linspace(ebn0_theory_start, ebn0_theory_end, nPointsTheory);

        % Calculate theoretical BER over the DENSE vector
        if strcmpi(ModType, 'qam') || strcmpi(ModType, 'pam')
            BER_theoretical_smooth = berawgn(EbN0_theory_dB, ModType, ModOrder, 'nondiff');

            semilogy(EbN0_theory_dB, BER_theoretical_smooth,...
                     'LineStyle', theoryLineStyles{theory_style_idx},...
                     'Marker', theoryMarker,...
                     'Color', theoryColors{theory_color_idx},...
                     'LineWidth', theoryLineWidth,...
                     'DisplayName', theoryDisplayName);
             legendEntries{end+1} = theoryDisplayName; % Add to legend list
        else
            warning('plotMultiBERCurves: Unrecognized ModType "%s" for %d-order. Skipping theoretical curve.', ModType, ModOrder);
            legendEntries{end+1} = sprintf('Theory: %s (N/A)', mod_label_base); % Add placeholder to legend
        end

    end % End loop through results

    % =====================================================================
    % == Plot Formatting ==
    % =====================================================================
    grid(gridVisible);
    xlabel(xAxisLabel);
    ylabel(yAxisLabel);
    title(plotTitle);
    % Use the collected legend entries
    legend(legendEntries, 'Location', legendLocation, 'Interpreter', 'none');

    % Adjust y-axis limits based on OVERALL simulated data and floor
    if min_sim_ber_overall >= 1 || isinf(min_sim_ber_overall) % If no valid BER points found below 1
         plot_lower_limit = yAxisMinBERFloor;
    else
         % Set lower limit slightly below min simulated BER, but not below floor.
         plot_lower_limit = max(min_sim_ber_overall / 10, yAxisMinBERFloor);
    end
    ylim([plot_lower_limit 1]); % Upper limit often 1 for BER

    % Set x-axis limits based on the overall simulation Eb/N0 range
    if ~isinf(min_ebn0_overall) && ~isinf(max_ebn0_overall)
       xlim([min_ebn0_overall max_ebn0_overall]);
    end

    hold off;
end