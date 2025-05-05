function hFig = plotConstellation_Tx_Rx(ModulationOrder, ModulationType, symb_tx, symb_rx)
    %PLOTCONSTELLATION_TX_RX Plots transmitted and received signal constellations.
    %   HFIG = PLOTCONSTELLATION_TX_RX(MODORDER, MODTYPE, SYMB_TX, SYMB_RX)
    %   creates a figure with two subplots showing the constellation diagrams
    %   of the transmitted symbols (SYMB_TX) and the received symbols (SYMB_RX)
    %   side-by-side.
    %
    %   This visualization helps in assessing the quality of the received signal
    %   by comparing the scattered received points to the ideal transmitted
    %   constellation points. Noise and distortion will cause the received points
    %   to deviate from the ideal locations.
    %
    %   Plotting appearance (colors, markers, titles, axes) is controlled by
    %   internal parameters. Assumes inputs are valid complex vectors.
    %
    %   Inputs:
    %       ModulationOrder - An integer indicating the modulation order (M), e.g.,
    %                         4 for QPSK, 16 for 16-QAM, 64 for 64-QAM. Used for title.
    %       ModulationType  - A string specifying the modulation type ('qam' or 'pam').
    %                         Used for title formatting.
    %       symb_tx         - A complex vector containing the transmitted symbols.
    %                         Each element represents a point in the ideal constellation.
    %       symb_rx         - A complex vector containing the received symbols,
    %                         typically after matched filtering and downsampling, but
    %                         before slicing/decision. These points show the effect
    %                         of channel noise and impairments.
    %
    %   Output:
    %       hFig            - A handle to the created figure object.

    % =====================================================================
    % == Plotting Parameters (Internal Configuration) ==
    % =====================================================================
    % --- Colors ---
    txColor = [0, 0.4470, 0.7410];      % Standard MATLAB blue for transmitted data
    rxColor = [0.8500, 0.3250, 0.0980]; % Standard MATLAB red for received data
    % --- Line/Marker Styles ---
    lineWidth = 1.5;                    % Line width (primarily for axes/markers if applicable)
    baseMarkerSize = 6;                 % Base size for scatter plot markers
    txMarkerStyle       = 'o';          % Marker style for transmitted points ('o', '+', '*', etc.)
    rxMarkerStyle       = 'x';          % Marker style for received points ('x', '.', 's', etc.)
    markerSizeMultiplier= 1;            % Multiplier for baseMarkerSize (optional scaling)
    % --- General Figure Settings ---
    figureName          = 'Constellation Diagrams (Tx vs Rx)'; % Figure window title
    figureNumberTitle   = 'off';        % Hide 'Figure X:' prefix ('on' or 'off')
    overallTitle        = 'Transmitted vs. Received Constellations'; % Main title for the entire figure
    % --- Subplot Titles ---
    subplotTitleFormat  = '%s %d-%s';   % Format string for subplot titles (e.g., "Transmitted 16-QAM")
    txTitlePrefix       = 'Transmitted';% Prefix for the transmitted constellation title
    rxTitlePrefix       = 'Received';   % Prefix for the received constellation title
    % --- Scatter Plot Settings ---
    % (Marker styles defined above)
    % --- Axes Properties ---
    xAxisLabel          = 'In-Phase (I)'; % Label for the x-axis (Real part)
    yAxisLabel          = 'Quadrature (Q)';% Label for the y-axis (Imaginary part)
    axisLimitPaddingFactor = 0.1;       % Factor (proportion of max abs value) for padding axis limits
    axisLimitPaddingConst = 0.1;        % Constant value added for padding axis limits (ensures some padding even if max is small)
    axisLineStyle       = 'k--';        % Style for the zero axes lines ('-', '--', ':', etc.)
    axisLineWidth       = 0.5;          % Line width for the zero axes lines
    enableGridMinor     = true;         % Enable minor grid lines ('true' or 'false')
    useAxisEqual        = true;         % Use 'axis equal' to ensure aspect ratio is 1:1 ('true' or 'false')
    linkConstAxes       = true;         % Link axes of the two subplots ('xy', 'x', 'y', or false/true for 'xy')
    % =====================================================================

    % --- Figure and Axes Creation ---
    % Create the main figure window.
    hFig = figure('Name', figureName, 'NumberTitle', figureNumberTitle);
    % Add an overall title to the figure.
    sgtitle(hFig, overallTitle);

    % --- Transmitted Constellation Subplot ---
    % Create the first subplot (left) for the transmitted symbols.
    ax_tx = subplot(1, 2, 1, 'Parent', hFig);
    % Format the title string for the Tx subplot.
    title_str_tx = sprintf(subplotTitleFormat, txTitlePrefix, ModulationOrder, upper(ModulationType));
    % Extract real (In-Phase) and imaginary (Quadrature) parts of Tx symbols.
    xData_tx = real(symb_tx);
    yData_tx = imag(symb_tx);
    % Create the scatter plot for transmitted symbols.
    scatter(ax_tx, xData_tx, yData_tx, baseMarkerSize * markerSizeMultiplier, ...
            txColor, txMarkerStyle, ...
            'MarkerFaceColor', txColor, 'LineWidth', lineWidth); % Use MarkerFaceColor for filled markers if style supports it
    hold(ax_tx, 'on'); % Hold plot for adding axes lines

    % Calculate axis limits with padding for Tx plot.
    minX_tx = min(xData_tx); maxX_tx = max(xData_tx);
    minY_tx = min(yData_tx); maxY_tx = max(yData_tx);
    % Determine padding based on both factor and constant, take the max to ensure sufficient space.
    padX_tx = max(abs([minX_tx, maxX_tx]) * axisLimitPaddingFactor + axisLimitPaddingConst);
    padY_tx = max(abs([minY_tx, maxY_tx]) * axisLimitPaddingFactor + axisLimitPaddingConst);
    % Calculate final axis limits.
    xLims_tx = [minX_tx - padX_tx, maxX_tx + padX_tx];
    yLims_tx = [minY_tx - padY_tx, maxY_tx + padY_tx];

    % Add horizontal and vertical lines at zero.
    plot(ax_tx, xLims_tx, [0 0], axisLineStyle, 'LineWidth', axisLineWidth, 'HandleVisibility','off'); % Hide from legend
    plot(ax_tx, [0 0], yLims_tx, axisLineStyle, 'LineWidth', axisLineWidth, 'HandleVisibility','off'); % Hide from legend
    % Configure grid, aspect ratio, limits, labels, title.
    grid(ax_tx, 'on');
    if enableGridMinor
        grid(ax_tx, 'minor');
    end
    if useAxisEqual
        axis(ax_tx, 'equal'); % Ensure I and Q axes have the same scale per unit
    end
    xlim(ax_tx, xLims_tx);
    ylim(ax_tx, yLims_tx);
    xlabel(ax_tx, xAxisLabel);
    ylabel(ax_tx, yAxisLabel);
    title(ax_tx, title_str_tx);
    box(ax_tx, 'on'); % Draw a box around the plot area
    hold(ax_tx, 'off'); % Release plot hold

    % --- Received Constellation Subplot ---
    % Create the second subplot (right) for the received symbols.
    ax_rx = subplot(1, 2, 2, 'Parent', hFig);
    % Format the title string for the Rx subplot.
    title_str_rx = sprintf(subplotTitleFormat, rxTitlePrefix, ModulationOrder, upper(ModulationType));
    % Extract real (In-Phase) and imaginary (Quadrature) parts of Rx symbols.
    xData_rx = real(symb_rx);
    yData_rx = imag(symb_rx);
    % Create the scatter plot for received symbols.
    scatter(ax_rx, xData_rx, yData_rx, baseMarkerSize * markerSizeMultiplier, ...
            rxColor, rxMarkerStyle, ...
            'MarkerFaceColor', rxColor, 'LineWidth', lineWidth); % Use MarkerFaceColor if applicable
    hold(ax_rx, 'on'); % Hold plot for adding axes lines

    % Calculate axis limits with padding for Rx plot (can be different from Tx due to noise).
    minX_rx = min(xData_rx); maxX_rx = max(xData_rx);
    minY_rx = min(yData_rx); maxY_rx = max(yData_rx);
    padX_rx = max(abs([minX_rx, maxX_rx]) * axisLimitPaddingFactor + axisLimitPaddingConst);
    padY_rx = max(abs([minY_rx, maxY_rx]) * axisLimitPaddingFactor + axisLimitPaddingConst);
    xLims_rx = [minX_rx - padX_rx, maxX_rx + padX_rx];
    yLims_rx = [minY_rx - padY_rx, maxY_rx + padY_rx];

    % Add horizontal and vertical lines at zero.
    plot(ax_rx, xLims_rx, [0 0], axisLineStyle, 'LineWidth', axisLineWidth, 'HandleVisibility','off'); % Hide from legend
    plot(ax_rx, [0 0], yLims_rx, axisLineStyle, 'LineWidth', axisLineWidth, 'HandleVisibility','off'); % Hide from legend
    % Configure grid, aspect ratio, limits, labels, title.
    grid(ax_rx, 'on');
    if enableGridMinor
        grid(ax_rx, 'minor');
    end
    if useAxisEqual
        axis(ax_rx, 'equal');
    end
    xlim(ax_rx, xLims_rx);
    ylim(ax_rx, yLims_rx);
    xlabel(ax_rx, xAxisLabel);
    ylabel(ax_rx, yAxisLabel);
    title(ax_rx, title_str_rx);
    box(ax_rx, 'on'); % Draw a box around the plot area
    hold(ax_rx, 'off'); % Release plot hold

    % --- Link Axes ---
    % Link the axes of the two subplots if requested (for synchronized zoom/pan).
    if linkConstAxes
        if ischar(linkConstAxes) || isstring(linkConstAxes)
             % Link specified axes ('x', 'y', or 'xy')
             linkaxes([ax_tx, ax_rx], linkConstAxes);
        elseif linkConstAxes % Treat boolean true as 'xy' linking
             % Link both x and y axes
             linkaxes([ax_tx, ax_rx], 'xy');
        end
        % If linkConstAxes is false, do nothing.
    end

end