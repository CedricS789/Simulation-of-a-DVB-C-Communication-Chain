function hFig = plotBitstream_Tx_Rx(bit_tx, bit_rx, bits_to_plot)
    %PLOTBITSTREAM_TX_RX Plots transmitted and received bit streams for comparison.
    %
    %   Inputs:
    %       bit_tx          - Transmitted binary vector.
    %       bit_rx          - Received binary vector.
    %       bits_to_plot    - Number of initial bits to display.
    %
    %   Output:
    %       hFig            - Handle to the created figure object.

    % =====================================================================
    % == Plotting Parameters (Internal Configuration) ==
    % =====================================================================
    txColor = [0, 0.4470, 0.7410];      % Blue
    rxColor = [0.8500, 0.3250, 0.0980]; % Red
    lineWidth = 1.5;
    figureName          = 'Bit Stream Comparison (Tx vs Rx)';
    figureNumberTitle   = 'off';
    overallTitle        = 'Transmitted vs. Received Bit Streams';
    txSubplotTitleFormat = 'Transmitted Bit Stream (First %d Bits)';
    rxSubplotTitleFormat = 'Received Bit Stream (First %d Bits)';
    xAxisLabel          = 'Bit Index';
    yAxisLabel          = 'Bit Value';
    xAxisLimitPadding   = 0.5;
    yAxisLimits         = [-0.1, 1.1];  % Fixed limits for 0/1 values
    yAxisTicks          = [0, 1];       % Explicit ticks at 0 and 1
    linkXAxes           = true;         % Link x-axes of subplots
    % =====================================================================

    % --- Input Validation ---
    actual_len_tx = length(bit_tx);
    actual_len_rx = length(bit_rx);
    bits_to_plot = min([bits_to_plot, actual_len_tx, actual_len_rx]);

    if bits_to_plot < 1
       warning('plotBitstream_Tx_Rx: No bits available to plot.');
       hFig = figure('Visible','off');
       return;
    end

    % --- Figure and Axes Creation ---
    hFig = figure('Name', figureName, 'NumberTitle', figureNumberTitle);
    sgtitle(hFig, overallTitle);
    plot_indices = 1:bits_to_plot;

    % --- Transmitted Bit Stream Subplot ---
    ax_tx = subplot(2, 1, 1, 'Parent', hFig);
    title_str_tx = sprintf(txSubplotTitleFormat, bits_to_plot);
    stairs(ax_tx, plot_indices, bit_tx(plot_indices), 'Color', txColor, 'LineWidth', lineWidth);
    title(ax_tx, title_str_tx);
    xlim(ax_tx, [plot_indices(1) - xAxisLimitPadding, plot_indices(end) + xAxisLimitPadding]);
    ylim(ax_tx, yAxisLimits);
    yticks(ax_tx, yAxisTicks);
    grid(ax_tx, 'on');
    xlabel(ax_tx, xAxisLabel);
    ylabel(ax_tx, yAxisLabel);
    box(ax_tx, 'on');

    % --- Received Bit Stream Subplot ---
    ax_rx = subplot(2, 1, 2, 'Parent', hFig);
    title_str_rx = sprintf(rxSubplotTitleFormat, bits_to_plot);
    stairs(ax_rx, plot_indices, bit_rx(plot_indices), 'Color', rxColor, 'LineWidth', lineWidth);
    title(ax_rx, title_str_rx);
    xlim(ax_rx, [plot_indices(1) - xAxisLimitPadding, plot_indices(end) + xAxisLimitPadding]);
    ylim(ax_rx, yAxisLimits);
    yticks(ax_rx, yAxisTicks);
    grid(ax_rx, 'on');
    xlabel(ax_rx, xAxisLabel);
    ylabel(ax_rx, yAxisLabel);
    box(ax_rx, 'on');

    % --- Link Axes ---
    if linkXAxes
        linkaxes([ax_tx, ax_rx], 'x');
    end

end