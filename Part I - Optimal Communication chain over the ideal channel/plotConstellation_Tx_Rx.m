function hFig = plotConstellation_Tx_Rx(ModulationOrder, ModulationType, symb_tx, symb_rx)
    %PLOTCONSTELLATION_TX_RX Plots transmitted and received signal constellations.
    %
    %   Inputs:
    %       ModulationOrder - Modulation order (M).
    %       ModulationType  - Modulation type string ('qam' or 'pam').
    %       symb_tx         - Transmitted complex symbols vector.
    %       symb_rx         - Received complex symbols vector.
    %
    %   Output:
    %       hFig            - Handle to the created figure object.

    % =====================================================================
    % == Plotting Parameters ==
    % =====================================================================
    txColor = [0, 0.4470, 0.7410];      % Blue
    rxColor = [0.8500, 0.3250, 0.0980]; % Red
    lineWidth = 1.5;
    baseMarkerSize = 60;
    txMarkerStyle       = 'o';
    rxMarkerStyle       = 'x';
    markerSizeMultiplier= 1;
    figureName          = 'Constellation Diagrams (Tx vs Rx)';
    figureNumberTitle   = 'off';
    overallTitle        = 'Transmitted vs. Received Constellations';
    subplotTitleFormat  = '%s %d-%s';   % e.g., "Transmitted 16-QAM"
    txTitlePrefix       = 'Transmitted';
    rxTitlePrefix       = 'Received';
    xAxisLabel          = 'In-Phase (I)';
    yAxisLabel          = 'Quadrature (Q)';
    axisLimitPaddingFactor = 0.1;       % Proportional padding for axis limits
    axisLimitPaddingConst = 0.1;        % Constant padding for axis limits
    axisLineStyle       = 'k--';        % Style for zero axes lines
    axisLineWidth       = 0.5;
    enableGridMinor     = true;
    useAxisEqual        = true;         % Ensure 1:1 aspect ratio
    linkConstAxes       = true;         % Link axes of subplots ('xy', 'x', 'y', or false/true for 'xy')

    % --- Figure and Axes Creation ---
    hFig = figure('Name', figureName, 'NumberTitle', figureNumberTitle);
    sgtitle(hFig, overallTitle);

    % --- Transmitted Constellation Subplot ---
    ax_tx = subplot(1, 2, 1, 'Parent', hFig);
    title_str_tx = sprintf(subplotTitleFormat, txTitlePrefix, ModulationOrder, upper(ModulationType));
    xData_tx = real(symb_tx);
    yData_tx = imag(symb_tx);
    scatter(ax_tx, xData_tx, yData_tx, baseMarkerSize * markerSizeMultiplier, ...
            txColor, txMarkerStyle, 'MarkerFaceColor', txColor, 'LineWidth', lineWidth);
    hold(ax_tx, 'on');

    % Calculate axis limits with padding for Tx plot
    minX_tx = min(xData_tx); maxX_tx = max(xData_tx);
    minY_tx = min(yData_tx); maxY_tx = max(yData_tx);
    padX_tx = max(abs([minX_tx, maxX_tx]) * axisLimitPaddingFactor + axisLimitPaddingConst);
    padY_tx = max(abs([minY_tx, maxY_tx]) * axisLimitPaddingFactor + axisLimitPaddingConst);
    xLims_tx = [minX_tx - padX_tx, maxX_tx + padX_tx];
    yLims_tx = [minY_tx - padY_tx, maxY_tx + padY_tx];

    plot(ax_tx, xLims_tx, [0 0], axisLineStyle, 'LineWidth', axisLineWidth, 'HandleVisibility','off'); % Zero axes
    plot(ax_tx, [0 0], yLims_tx, axisLineStyle, 'LineWidth', axisLineWidth, 'HandleVisibility','off');
    grid(ax_tx, 'on');
    if enableGridMinor, grid(ax_tx, 'minor'); end
    if useAxisEqual, axis(ax_tx, 'equal'); end
    xlim(ax_tx, xLims_tx);
    ylim(ax_tx, yLims_tx);
    xlabel(ax_tx, xAxisLabel);
    ylabel(ax_tx, yAxisLabel);
    title(ax_tx, title_str_tx);
    box(ax_tx, 'on');
    hold(ax_tx, 'off');

    % --- Received Constellation Subplot ---
    ax_rx = subplot(1, 2, 2, 'Parent', hFig);
    title_str_rx = sprintf(subplotTitleFormat, rxTitlePrefix, ModulationOrder, upper(ModulationType));
    xData_rx = real(symb_rx);
    yData_rx = imag(symb_rx);
    scatter(ax_rx, xData_rx, yData_rx, baseMarkerSize * markerSizeMultiplier, ...
            rxColor, rxMarkerStyle, 'MarkerFaceColor', rxColor, 'LineWidth', lineWidth);
    hold(ax_rx, 'on');

    % Calculate axis limits with padding for Rx plot
    minX_rx = min(xData_rx); maxX_rx = max(xData_rx);
    minY_rx = min(yData_rx); maxY_rx = max(yData_rx);
    padX_rx = max(abs([minX_rx, maxX_rx]) * axisLimitPaddingFactor + axisLimitPaddingConst);
    padY_rx = max(abs([minY_rx, maxY_rx]) * axisLimitPaddingFactor + axisLimitPaddingConst);
    xLims_rx = [minX_rx - padX_rx, maxX_rx + padX_rx];
    yLims_rx = [minY_rx - padY_rx, maxY_rx + padY_rx];

    plot(ax_rx, xLims_rx, [0 0], axisLineStyle, 'LineWidth', axisLineWidth, 'HandleVisibility','off'); % Zero axes
    plot(ax_rx, [0 0], yLims_rx, axisLineStyle, 'LineWidth', axisLineWidth, 'HandleVisibility','off');
    grid(ax_rx, 'on');
    if enableGridMinor, grid(ax_rx, 'minor'); end
    if useAxisEqual, axis(ax_rx, 'equal'); end
    xlim(ax_rx, xLims_rx);
    ylim(ax_rx, yLims_rx);
    xlabel(ax_rx, xAxisLabel);
    ylabel(ax_rx, yAxisLabel);
    title(ax_rx, title_str_rx);
    box(ax_rx, 'on');
    hold(ax_rx, 'off');

    % --- Link Axes ---
    if linkConstAxes
        if ischar(linkConstAxes) || isstring(linkConstAxes) || linkConstAxes
            linkaxes([ax_tx, ax_rx], linkConstAxes); % Link specified axes ('x', 'y', or 'xy'), or both if true
        end
    end

end