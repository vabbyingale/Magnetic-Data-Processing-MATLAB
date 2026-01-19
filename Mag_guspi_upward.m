function upward_field = Mag_guspi_upward(input_field, fish_depth, dx, dy, wl, ws, zlev, pad_factor, OPT)
% GUSPI_UPWARD_CONTINUATION - Upward continuation using Guspi's 1987 method
% with a Tukey-window inspired cosine-tapered bandpass filter and controlled convergence
%
% Inputs:
%   input_field - 2D matrix of magnetic field measurements
%   fish_depth  - depth of measurements (scalar or matrix, positive downward)
%   dx, dy      - grid spacing (x and y directions)
%   wl          - long wavelength cutoff (low frequency bound)
%   ws          - short wavelength cutoff (high frequency bound)
%   zlev        - target level for upward continuation (positive downward)
%   OPT         - [max_iter, max_terms, gtol, tol, alpha] (default: [20, 20, 0.01, 0.01, 0.5])
%
% Output:
%   upward_field - upward continued magnetic field

    % Defaults
    if nargin < 9
        OPT = [300, 300, 0.01, 0.01];
    end
    [max_iter, max_terms, gtol, tol] = deal(OPT(1), OPT(2), OPT(3), OPT(4));
    [ny, nx] = size(input_field);

    % Symmetric padding
    pad_size = round(pad_factor * min([ny, nx]));  % symmetric padding
    input_padded = padarray(input_field, [pad_size pad_size], 'symmetric');
    depth_padded = padarray(fish_depth, [pad_size pad_size], 'symmetric');
    [ny_pad, nx_pad] = size(input_padded);

    % Reference depth
    zref = mean(depth_padded(:));
    h = depth_padded - zref;
    dz = zlev - zref;

    % Mean center
    mean_input = mean(input_padded(:));
    field_centered = input_padded - mean_input;

    % Wavelength limits
    if ws == 0, ws = max([dx, dy]); end
    if wl == 0, wl = min(nx * dx, ny * dy); end
    if ws >= wl
        error('Short wavelength cutoff (%.2f) must be less than long wavelength cutoff (%.2f)', ws, wl);
    end

    % Wavenumber grid
    kx = ifftshift( (-floor(nx_pad/2):ceil(nx_pad/2)-1) * (2*pi / (nx_pad*dx)) );
    ky = ifftshift( (-floor(ny_pad/2):ceil(ny_pad/2)-1) * (2*pi / (ny_pad*dy)) );
    [KX, KY] = meshgrid(kx, ky);
    K = sqrt(KX.^2 + KY.^2);

    % Bandpass filter (cosine-tapered)
    f_low = 1 / wl;
    f_high = 1 / ws;
    k_low = 2*pi*f_low;
    k_high = 2*pi*f_high;
    dk = k_high - k_low;
    dk_taper = 0.3 * dk;
    k1 = k_low - dk_taper;
    k2 = k_low;
    k3 = k_high;
    k4 = k_high + dk_taper;

    W = zeros(size(K));
    idx1 = (K >= k1) & (K < k2);
    W(idx1) = 0.5 * (1 + cos(pi * (K(idx1) - k2) / dk_taper));
    idx2 = (K >= k2) & (K <= k3);
    W(idx2) = 1;
    idx3 = (K > k3) & (K <= k4);
    W(idx3) = 0.5 * (1 + cos(pi * (K(idx3) - k3) / dk_taper));

    % Precompute factorial terms
    factorial_terms = arrayfun(@factorial, 1:max_terms);

    % Guspi iteration
    m_ref = field_centered;
    convergence_history = NaN(max_iter, 1);

    % For rollback logic
    best_m_ref = m_ref;
    best_rms_err = inf;
    max_bad_iters = 3;
    bad_iter_count = 0;

    fprintf('Starting Guspi upward continuation\n');
    fprintf('  Depth range: ref %.2f m to target %.2f m (dz = %.2f m)\n', zref, zlev, dz);
    fprintf('  Wavelength band: %.2f m - %.2f m\n', ws, wl);

    if any(isnan(field_centered(:))) || any(isinf(field_centered(:)))
        error('Input field contains NaN or Inf');
    end
    if any(isnan(h(:))) || any(isinf(h(:)))
        error('Depths (fish_depth) contain NaN or Inf');
    end
    if any(isnan(W(:))) || any(isinf(W(:)))
        error('Bandpass filter contains NaN or Inf');
    end
    if any(isnan(K(:))) || any(isinf(K(:)))
        error('Wavenumber matrix K contains NaN or Inf');
    end


    for iter = 1:max_iter
        M_fft = fft2(m_ref);
        taylor_sum = zeros(size(m_ref));

        for n = 1:max_terms
            h_term = ((-h).^n) / factorial_terms(n);
            freq_term = (K.^n) .* M_fft .* W;
            spatial_term = real(ifft2(freq_term));
            contribution = h_term .* spatial_term;
            taylor_sum = taylor_sum + contribution;

            if max(abs(contribution(:))) < gtol
                fprintf('  Iter %d: Taylor converged at term %d\n', iter, n);
                break;
            end
        end

        m_new = field_centered - taylor_sum;

        % RMS error relative to previous
        rms_err = rms(m_new(:) - m_ref(:));
        convergence_history(iter) = rms_err;

        fprintf('  Iter %2d: RMS error = %.6e\n', iter, rms_err);

        if any(isnan(m_new(:))) || any(isinf(m_new(:)))
            warning('NaN or Inf detected — aborting');
            upward_field = best_m_ref + mean_input;
            return;
        end

        % Track best and tolerate a few bad steps
        if rms_err < best_rms_err
            best_rms_err = rms_err;
            best_m_ref = m_new;
            bad_iter_count = 0;
        else
            bad_iter_count = bad_iter_count + 1;
            fprintf('  Iteration %d RMS increased (bad step %d of %d)\n', iter, bad_iter_count, max_bad_iters);
            if bad_iter_count >= max_bad_iters
                fprintf('Too many bad steps — reverting to best result and stopping\n');
                upward_field = best_m_ref + mean_input;
                return;
            end
        end

        % Absolute convergence check
        if rms_err < tol
            fprintf('Converged at iteration %d (RMS < %.2e)\n', iter, tol);
            break;
        end

        % Update for next iteration
        m_ref = m_new;
    end

    % Final upward continuation in frequency domain
    upward_filter = exp(-K * dz);
    m_continued = real(ifft2(fft2(best_m_ref) .* upward_filter));
    upward_field_padded = m_continued + mean_input;

    % Crop to original size
    upward_field = upward_field_padded(pad_size+1:end-pad_size, pad_size+1:end-pad_size);
end
