% Raised Cosine (RC) Filter ISI Cancellation Demonstration
% This script shows the time-domain response of an RC filter 
% and demonstrates its ISI cancellation property

% Parameters
rolloff = 0.35;  % Rolloff factor (0 <= α <= 1)
symbol_rate = 1; % Symbol rate
samples_per_symbol = 8; % Oversampling factor
span = 6;        % Filter span in symbol periods

% Time vector
t = linspace(-span/2, span/2, span * samples_per_symbol + 1);

% Raised Cosine Filter Impulse Response
h = zeros(size(t));
for i = 1:length(t)
    if t(i) == 0
        h(i) = 1;
    elseif abs(t(i)) == 1 / (2 * symbol_rate)
        h(i) = pi / (4 * symbol_rate);
    else
        h(i) = sin(pi * t(i) * symbol_rate) ./ (pi * t(i) * symbol_rate) .* ...
               cos(pi * rolloff * t(i) * symbol_rate) ./ ...
               (1 - (2 * rolloff * t(i) * symbol_rate).^2);
    end
end

% Normalize the filter
h = h / sum(h);

% Simulate multiple symbol pulses
symbols = [1, -1, 1, 1, -1];  % Example symbol sequence
pulse_train = zeros(1, length(symbols) * samples_per_symbol);

% Create pulse train
for i = 1:length(symbols)
    pulse_train((i-1)*samples_per_symbol + 1) = symbols(i);
end

% Convolve pulse train with RC filter
filtered_signal = conv(pulse_train, h, 'same');

% Plotting
figure('Position', [100, 100, 1000, 600]);

% Subplot 1: Impulse Response
subplot(2,2,1);
plot(t, h, 'LineWidth', 2);
title('RC Filter Impulse Response');
xlabel('Time');
ylabel('Amplitude');
grid on;

% Subplot 2: Original Pulse Train
subplot(2,2,2);
stem(pulse_train, 'filled', 'LineWidth', 2);
title('Original Pulse Train');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;

% Subplot 3: Filtered Signal
subplot(2,2,3);
plot(filtered_signal, 'LineWidth', 2);
title('Filtered Signal');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;

% Subplot 4: Sampling Points to Show ISI Cancellation
subplot(2,2,4);
sample_indices = (0:length(symbols)-1) * samples_per_symbol + samples_per_symbol/2;
sampled_values = filtered_signal(sample_indices);
stem(sample_indices, sampled_values, 'filled', 'LineWidth', 2);
title('Sampled Points (ISI Cancellation)');
xlabel('Sample Index');
ylabel('Amplitude');
ylim([-1.5, 1.5]);
grid on;

% Adjust overall figure
sgtitle('Raised Cosine (RC) Filter - ISI Cancellation', 'FontSize', 14);