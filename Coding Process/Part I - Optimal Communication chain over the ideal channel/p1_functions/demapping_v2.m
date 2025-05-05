function [bit_rx] = demapping_v2(symb_rx, Nbps, modulation)
% Demaps received symbols back to bits.
%   Handles both PAM and QAM modulations with Gray coding.
%   Includes fixes for common issues like de2bi input size and
%   handling complex inputs for PAM after noise addition.
%
% INPUTS:
% - symb_rx : Column vector of received complex symbols (potentially noisy).
% - Nbps : Number of bits per symbol (integer).
% - modulation : String specifying modulation type ('pam' or 'qam').
%
% OUTPUTS:
% - bit_rx : Column vector of estimated output bits.

    % --- Input Validation and Setup ---
    % Ensure symb_rx is a column vector for consistent processing
    symb_rx = symb_rx(:);
    Nsymb = size(symb_rx, 1); % Number of symbols

    % Handle empty input case gracefully
    if Nsymb == 0
        bit_rx = [];
        return;
    end

    % --- Demapping based on Modulation Type ---
    switch lower(modulation) % Use lower case for robustness

        case 'pam'
            % PAM demodulation relies only on the real part of the symbol.
            M = 2^Nbps; % Modulation order

            % --- Calculate Scaling Factor ---
            % Define ideal integer levels and find the standard deviation
            % used for normalization during mapping.
            pam_levels = 0:(M-1);
            mean_level = mean(pam_levels); % (M-1)/2
            sigma = sqrt(mean((pam_levels - mean_level).^2));
            % sigma = sqrt(sum(([0:M-1]-(M-1)/2).^2)/M); % Original equivalent calculation

            % --- Symbol to Integer Conversion ---
            % Invert the mapping: symbol = (1/sigma) * (integer - mean_level)
            % => integer = sigma * symbol + mean_level
            % *** FIX for AWGN ***: Use only the REAL part of the received symbol
            % because AWGN makes symb_rx complex, but PAM decision is 1D (real).
            int_rx = sigma * real(symb_rx) + mean_level;

            % --- Hard Decision Integer Detection ---
            % Round the calculated real value to the nearest valid integer.
            int_det = round(int_rx);

            % Clip to the valid integer range [0, M-1] to handle noise driving
            % values outside the expected range.
            int_det(int_det < 0) = 0;
            int_det(int_det > M-1) = M-1;

            % --- Integer to Binary Conversion (Natural Binary) ---
            % Convert detected integers to their binary representations.
            % *** FIX ***: Specify the number of bits (Nbps) required.
            % Omitting this causes de2bi to use the minimum bits for the max
            % value *present* in int_det, which can be less than Nbps.
            % 'left-msb' ensures standard binary order (e.g., 6 -> [1 1 0]).
            try
                mapp_rx = de2bi(int_det, Nbps, 'left-msb');
            catch ME % Catch potential errors from de2bi if int_det is bad
                 fprintf('Error in de2bi for PAM: %s\n', ME.message);
                 fprintf('Size of int_det: %d x %d\n', size(int_det,1), size(int_det,2));
                 fprintf('Min/Max of int_det: %f / %f\n', min(int_det), max(int_det));
                 fprintf('Nbps: %d\n', Nbps);
                 rethrow(ME); % Stop execution
            end

            % --- Sanity Check for de2bi Output Size ---
            % This check should technically not be needed after specifying Nbps
            % and clipping int_det, but serves as a safeguard.
            if size(mapp_rx, 1) ~= Nsymb || size(mapp_rx, 2) ~= Nbps
                % Handle potential edge case where input was valid but de2bi failed.
                % If int_det was empty, mapp_rx would be 0x0, need 0xNbps.
                 if Nsymb == 0 && isempty(mapp_rx)
                     mapp_rx = zeros(0, Nbps); % Correct size for empty input
                 elseif size(mapp_rx, 2) < Nbps % Pad if somehow too few columns
                    warning('Padding de2bi output in PAM - check inputs.');
                    mapp_rx = [zeros(Nsymb, Nbps - size(mapp_rx, 2)), mapp_rx];
                 elseif size(mapp_rx, 1) ~= Nsymb
                     error('Demapping error: Row count mismatch after de2bi (PAM).');
                 end
            end

            % --- Binary to Gray Conversion ---
            % Convert the natural binary mapping back to Gray coded bits.
            bit_rx2 = zeros(Nsymb, Nbps); % Pre-allocate for speed
            bit_rx2(:, 1) = mapp_rx(:, 1); % First bit is the same
            for ii = 2:Nbps
                % Gray bit = current binary XOR previous binary
                bit_rx2(:, ii) = xor(mapp_rx(:, ii - 1), mapp_rx(:, ii));
            end

            % --- Reshape to Output Bit Vector ---
            % Convert the Nsymb x Nbps matrix into a single column vector.
            % *** FIX ***: Use `[]` for the second dimension argument in reshape.
            % This makes it robust, automatically calculating the dimension (1)
            % based on the total number of elements. Explicit '1' failed for Nbps=1.
            bit_rx = reshape(bit_rx2', Nsymb * Nbps, []); % Transpose first!

        case 'qam'
            % QAM demodulation involves separate processing for Real (I) and Imaginary (Q) parts.
            if mod(Nbps, 2) ~= 0
                error('Nbps must be even for QAM modulation.');
            end

            % --- Split Nbps for I and Q ---
            Nbps_per_dim = Nbps / 2; % Bits per dimension (I or Q)
            M_per_dim = 2^Nbps_per_dim; % PAM order per dimension

            % --- Calculate Scaling Factor (Same for I and Q) ---
            pam_levels_1D = 0:(M_per_dim - 1);
            mean_level_1D = mean(pam_levels_1D); % (M_per_dim - 1) / 2
            sigma_1D = sqrt(mean((pam_levels_1D - mean_level_1D).^2));
            % Power is split equally between I and Q for unit variance QAM,
            % hence the sqrt(2) factor during mapping inversion.

            % = = = = = = = = = = = REAL PART (I) = = = = = = = = = = =
            symb_rxI = real(symb_rx); % Extract real part

            % --- Symbol to Integer (I) ---
            int_rxI = sigma_1D * sqrt(2) * symb_rxI + mean_level_1D;

            % --- Integer Detection (I) ---
            int_detI = round(int_rxI);
            int_detI(int_detI < 0) = 0;
            int_detI(int_detI > M_per_dim - 1) = M_per_dim - 1;

            % --- Integer to Binary (I) ---
            % *** FIX ***: Specify Nbps_per_dim for de2bi.
            try
                mapp_rxI = de2bi(int_detI, Nbps_per_dim, 'left-msb');
            catch ME
                 fprintf('Error in de2bi for QAM-I: %s\n', ME.message);
                 rethrow(ME);
            end

            % --- Sanity Check (I) ---
            if size(mapp_rxI, 1) ~= Nsymb || size(mapp_rxI, 2) ~= Nbps_per_dim
                 if Nsymb == 0 && isempty(mapp_rxI)
                     mapp_rxI = zeros(0, Nbps_per_dim);
                 elseif size(mapp_rxI, 2) < Nbps_per_dim
                    warning('Padding de2bi output in QAM-I - check inputs.');
                    mapp_rxI = [zeros(Nsymb, Nbps_per_dim - size(mapp_rxI, 2)), mapp_rxI];
                 elseif size(mapp_rxI, 1) ~= Nsymb
                     error('Demapping error: Row count mismatch after de2bi (QAM-I).');
                 end
            end

            % --- Binary to Gray (I) ---
            bit_rx2I = zeros(Nsymb, Nbps_per_dim);
            bit_rx2I(:, 1) = mapp_rxI(:, 1);
            for ii = 2:Nbps_per_dim
                bit_rx2I(:, ii) = xor(mapp_rxI(:, ii - 1), mapp_rxI(:, ii));
            end

            % = = = = = = = = = = IMAGINARY PART (Q) = = = = = = = = = =
            symb_rxQ = imag(symb_rx); % Extract imaginary part

            % --- Symbol to Integer (Q) ---
            int_rxQ = sigma_1D * sqrt(2) * symb_rxQ + mean_level_1D;

            % --- Integer Detection (Q) ---
            int_detQ = round(int_rxQ);
            int_detQ(int_detQ < 0) = 0;
            int_detQ(int_detQ > M_per_dim - 1) = M_per_dim - 1;

            % --- Integer to Binary (Q) ---
             % *** FIX ***: Specify Nbps_per_dim for de2bi.
            try
                mapp_rxQ = de2bi(int_detQ, Nbps_per_dim, 'left-msb');
            catch ME
                 fprintf('Error in de2bi for QAM-Q: %s\n', ME.message);
                 rethrow(ME);
            end


            % --- Sanity Check (Q) ---
            if size(mapp_rxQ, 1) ~= Nsymb || size(mapp_rxQ, 2) ~= Nbps_per_dim
                 if Nsymb == 0 && isempty(mapp_rxQ)
                     mapp_rxQ = zeros(0, Nbps_per_dim);
                 elseif size(mapp_rxQ, 2) < Nbps_per_dim
                    warning('Padding de2bi output in QAM-Q - check inputs.');
                    mapp_rxQ = [zeros(Nsymb, Nbps_per_dim - size(mapp_rxQ, 2)), mapp_rxQ];
                 elseif size(mapp_rxQ, 1) ~= Nsymb
                     error('Demapping error: Row count mismatch after de2bi (QAM-Q).');
                 end
            end

            % --- Binary to Gray (Q) ---
            bit_rx2Q = zeros(Nsymb, Nbps_per_dim);
            bit_rx2Q(:, 1) = mapp_rxQ(:, 1);
            for ii = 2:Nbps_per_dim
                bit_rx2Q(:, ii) = xor(mapp_rxQ(:, ii - 1), mapp_rxQ(:, ii));
            end

            % --- Bit Concatenation & Reshape ---
            % Combine Gray bits from I and Q paths side-by-side, then reshape.
            bit_rx_combined = [bit_rx2I, bit_rx2Q]; % Nsymb x Nbps matrix
            bit_rx = reshape(bit_rx_combined', Nsymb * Nbps, []); % Transpose and reshape to column vector

        otherwise
            error('Unknown modulation type: %s. Use ''pam'' or ''qam''.', modulation);

    end

    % --- Final Output ---
    % Ensure output is always a column vector
    bit_rx = bit_rx(:);

end