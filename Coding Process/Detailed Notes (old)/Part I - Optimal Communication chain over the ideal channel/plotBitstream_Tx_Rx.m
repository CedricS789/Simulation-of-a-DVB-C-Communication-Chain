function hFig = plotBitstream_Tx_Rx(bit_tx, bit_rx, bits_to_plot)
    %PLOTBITSTREAM_TX_RX Plots transmitted and received bit streams for comparison.
    %   HFIG = PLOTBITSTREAM_TX_RX(BIT_TX, BIT_RX, BITS_TO_PLOT) creates a
    %   figure with two subplots. The top subplot shows the first BITS_TO_PLOT
    %   bits of the transmitted bit stream (BIT_TX), and the bottom subplot
    %   shows the first BITS_TO_PLOT bits of the received bit stream (BIT_RX).
    %   This visualization helps in qualitatively assessing the signal
    %   transmission by comparing the input and output bits side-by-side.
    %
    %   The plots use stairs visualization, suitable for binary data.
    %   Plotting appearance (colors, titles, labels) is controlled by internal parameters.
    %
    %   Inputs:
    %       bit_tx          - A binary vector (containing 0s and 1s) representing
    %                         the transmitted bit stream.
    %       bit_rx          - A binary vector representing the received bit stream
    %                         after demodulation and detection. Should ideally have
    %                         the same length as bit_tx, but comparison length is limited.
    %       bits_to_plot    - A positive integer specifying the number of initial bits
    %                         from both streams to display in the plots. If this value
    %                         exceeds the length of either bit stream, it will be
    %                         clipped to the maximum available length.
    %
    %   Output:
    %       hFig            - A handle to the created figure object, allowing further
    %                         manipulation if needed.

    % =====================================================================
    % == Plotting Parameters (Internal Configuration) ==
    % =====================================================================
    % --- Colors ---
    txColor = [0, 0.4470, 0.7410];      % Standard MATLAB blue for transmitted data
    rxColor = [0.8500, 0.3250, 0.0980]; % Standard MATLAB red for received data
    % --- Line ---
    lineWidth = 1.5;                    % Line width for the stairs plots
    % --- General Figure ---
    figureName          = 'Bit Stream Comparison (Tx vs Rx)'; % Figure window title
    figureNumberTitle   = 'off';        % Hide 'Figure X:' prefix in title ('on' or 'off')
    overallTitle        = 'Transmitted vs. Received Bit Streams'; % Main title for the entire figure
    % --- Subplot Titles ---
    txSubplotTitleFormat = 'Transmitted Bit Stream (First %d Bits)'; % Format for top subplot title
    rxSubplotTitleFormat = 'Received Bit Stream (First %d Bits)';   % Format for bottom subplot title
    % --- Axes Properties ---
    xAxisLabel          = 'Bit Index';  % Label for the common x-axis
    yAxisLabel          = 'Bit Value';  % Label for the common y-axis
    xAxisLimitPadding   = 0.5;          % Padding added to x-axis limits
    yAxisLimits         = [-0.1, 1.1];  % Fixed y-axis limits suitable for 0/1 values
    yAxisTicks          = [0, 1];       % Set y-axis ticks explicitly at 0 and 1
    linkXAxes           = true;         % Link the x-axes of the two subplots for synchronized zooming/panning
    % =====================================================================

    % --- Input Validation ---
    % Ensure the number of bits to plot does not exceed the length of the input vectors.
    actual_len_tx = length(bit_tx);
    actual_len_rx = length(bit_rx);
    bits_to_plot = min(bits_to_plot, actual_len_tx);
    bits_to_plot = min(bits_to_plot, actual_len_rx); % Use the shorter length if they differ

    % Check if there are any bits to plot after validation.
    if bits_to_plot < 1
       warning('plotBitstream_Tx_Rx: No bits available to plot or bits_to_plot was non-positive.');
       % Return an empty or invisible figure handle if nothing can be plotted.
       hFig = figure('Visible','off'); % Create an invisible figure to return a valid handle
       return;
    end

    % --- Figure and Axes Creation ---
    % Create the main figure window.
    hFig = figure('Name', figureName, 'NumberTitle', figureNumberTitle);
    % Add an overall title to the figure.
    sgtitle(hFig, overallTitle);

    % Define the range of indices to plot.
    plot_indices = 1:bits_to_plot;

    % --- Transmitted Bit Stream Subplot ---
    % Create the first subplot (top) for the transmitted bits.
    ax_tx = subplot(2, 1, 1, 'Parent', hFig);
    % Format the title string for the Tx subplot.
    title_str_tx = sprintf(txSubplotTitleFormat, bits_to_plot);
    % Plot the transmitted bits using stairs for clear visualization of levels.
    stairs(ax_tx, plot_indices, bit_tx(plot_indices), 'Color', txColor, 'LineWidth', lineWidth);
    % Set the title for the Tx subplot.
    title(ax_tx, title_str_tx);
    % Set x-axis limits with padding.
    xlim(ax_tx, [plot_indices(1) - xAxisLimitPadding, plot_indices(end) + xAxisLimitPadding]);
    % Set y-axis limits and ticks for binary data.
    ylim(ax_tx, yAxisLimits);
    yticks(ax_tx, yAxisTicks);
    % Enable grid lines.
    grid(ax_tx, 'on');
    % Set axis labels.
    xlabel(ax_tx, xAxisLabel);
    ylabel(ax_tx, yAxisLabel);
    % Draw a box around the plot area.
    box(ax_tx, 'on');

    % --- Received Bit Stream Subplot ---
    % Create the second subplot (bottom) for the received bits.
    ax_rx = subplot(2, 1, 2, 'Parent', hFig);
    % Format the title string for the Rx subplot.
    title_str_rx = sprintf(rxSubplotTitleFormat, bits_to_plot);
    % Plot the received bits using stairs.
    stairs(ax_rx, plot_indices, bit_rx(plot_indices), 'Color', rxColor, 'LineWidth', lineWidth);
    % Set the title for the Rx subplot.
    title(ax_rx, title_str_rx);
    % Set x-axis limits with padding.
    xlim(ax_rx, [plot_indices(1) - xAxisLimitPadding, plot_indices(end) + xAxisLimitPadding]);
    % Set y-axis limits and ticks.
    ylim(ax_rx, yAxisLimits);
    yticks(ax_rx, yAxisTicks);
    % Enable grid lines.
    grid(ax_rx, 'on');
    % Set axis labels.
    xlabel(ax_rx, xAxisLabel);
    ylabel(ax_rx, yAxisLabel);
    % Draw a box around the plot area.
    box(ax_rx, 'on');

    % --- Link Axes ---
    % Link the x-axes of the two subplots if requested.
    if linkXAxes
        linkaxes([ax_tx, ax_rx], 'x');
    end

end