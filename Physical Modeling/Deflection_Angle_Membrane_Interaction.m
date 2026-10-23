% Ben G. 2026
% N protein - condensate-membrane bending model

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Below are all of the custom MATLAB files used to model the
% condensate-membrane free energy landscape. 

function [E_surface] = surfaceEnergy(delta_yDM_yMS, yDS, lambda, theta, phi, R, k_bend)

    E_surface = (delta_yDM_yMS * 2 * pi * ((R.^2) .* (sin(phi).^2))/(sin(theta - phi).^2)) * (1 - cos(theta - phi))...
                + (yDS * 2 * pi * (R.^2) .* (1 - cos(phi)))...
                + (lambda * 2 * pi * R .* sin(phi))...
                + (4 * pi * k_bend * ((R.^2) .* (sin(phi).^2))/(sin(theta - phi).^2)) * (1 - cos(theta - phi))*(((R .* (sin(phi))/(sin(theta - phi))).^(-1)))^(2);

end

function [R_prime, phi, theta_intrinsic, theta_intrinsic_vector] = doubleSphericalCapConstraint(delta_yDM_yMS)

    theta_intrinsic = acos(-1*delta_yDM_yMS);
    phi = linspace(0.01,(theta_intrinsic-0.01),50); % we could let phi go to zero - but this would be hue
    R_prime = zeros(1,length(phi));
    theta_intrinsic_vector = zeros(1,length(phi));
    
    for i = 1:length(phi)
    
        psi = theta_intrinsic - phi(i);

        if psi == 0 % no membrane being condition - volume is completely contained in the top spherical cap. 

            R_prime(i) = (4/(2 - 3*cos(phi(i)) + cos(phi(i))^3))^(1/3);
            theta_intrinsic_vector(i) = theta_intrinsic;
    
        elseif psi <= (pi/12) % sin small angle approximation to avoid issue
    
            R_prime(i) = (4/((2 - 3*cos(phi(i)) + (cos(phi(i))^3)) + ((sin(phi(i))^3)/((psi - ((psi^3)/6) + ((psi^5)/120)))^3)*(2 - 3*cos(psi) + (cos(psi)^3))))^(1/3);
            theta_intrinsic_vector(i) = theta_intrinsic;

        else 
    
            R_prime(i) = (4/((2 - 3*cos(phi(i)) + (cos(phi(i))^3)) + ((sin(phi(i))^3)/(sin(psi)^3))*(2 - 3*cos(psi) + (cos(psi)^3))))^(1/3);
            theta_intrinsic_vector(i) = theta_intrinsic;
        
        end
    
    end
end

function [E_surface, phi, theta_intrinsic, R_prime] = energyForGivenContactAngle(delta_yDM_yMS, k_bend, lambda)

    % obtain contrainted geometries for a given contact angle
    [R_prime, phi, theta_intrinsic] = doubleSphericalCapConstraint(delta_yDM_yMS);

    %k_bend = 1; % normalized to yDS*R_0^2
    %lamba = 3; % normalized to yDS*R_0^2 CAN BE POSITIVE OF NEGATIVE!!

    [E_surface] = surfaceEnergy(delta_yDM_yMS, 1, lambda, theta_intrinsic, phi, R_prime, k_bend);

end

function [dE_dphi, phi_mid, k_bend] = energyFittingFromDeflectionAngle(theta_exp, phi_exp)

    theta = theta_exp*pi;
    lambda = 0.5;
    k_bend = 0.0;
    delta_y = -1*cos(theta);

    % Get the gradient of the energy landscape and the radius for a given contact angle
    [E_surface, phi, theta_intrinsic, R] = energyForGivenContactAngle(delta_y, k_bend, lambda);
    dE_dphi = diff(E_surface) ./ diff(phi);  % numeric derivative
    phi_mid  = (phi(1:end-1) + phi(2:end)) / 2;

    Bend_term_phi_coef = (4 * pi  * ((R.^2) .* (sin(phi).^2))/(sin(theta - phi).^2)) * (1 - cos(theta - phi))*(((R .* (sin(phi))/(sin(theta - phi))).^(-1)))^(2);
    dEbend_dphi = diff(Bend_term_phi_coef)./diff(phi);

    % Get the index that is closest to the experimentally measured membrane
    % deflection phi.
    
    differences = abs(phi - (phi_exp*pi));
    %indexes = zeros(length(differences));
    [minDiff, indexes] = min(differences);

    %phi_equilibrium = phi(indexes);
   
    % Get the k_bend asssuming energy equilibirium
    k_bend = -dE_dphi(indexes)/(dEbend_dphi(indexes));

end

function [energyMatrix, phiMatrix, thetaMatrix] = energyLandscape(theta_range, k_bend, lambda)
    % Inputs:
    %   theta_range: array of intrinsic contact angles
    %   k_bend, lambda: physical parameters
    %
    % Outputs:
    %   energyMatrix: 2D matrix of free energy values
    %   phiMatrix: 2D matrix of deflection angles (φ)
    %   thetaMatrix: 2D matrix of contact angles (θ), matched in size

    n = length(theta_range);
    % Preallocate lists (use cell first since phi length might vary slightly)
    energyCell = cell(1, n);
    phiCell = cell(1, n);
    thetaCell = cell(1, n);

    for i = 1:n
        theta = theta_range(i);

        % Evaluate energy and phi vector for given theta
        try
            [E_surface, phi_vec, theta_intrinsic_vec] = energyForGivenContactAngle(theta, k_bend, lambda);

            energyCell{i} = E_surface;
            phiCell{i} = phi_vec;
            thetaCell{i} = theta_intrinsic_vec;
        catch
            % In case of failure (e.g., NaN), fill with NaNs
            energyCell{i} = nan(size(phi_vec));
            phiCell{i} = phi_vec;
            thetaCell{i} = theta * ones(size(phi_vec));
        end
    end

    % Convert cells to full matrices using padding
    maxLength = max(cellfun(@length, phiCell));
    energyMatrix = nan(maxLength, n);
    phiMatrix = nan(maxLength, n);
    thetaMatrix = nan(maxLength, n);

    for i = 1:n
        len = length(phiCell{i});
        energyMatrix(1:len, i) = energyCell{i};
        phiMatrix(1:len, i) = phiCell{i};
        thetaMatrix(1:len, i) = thetaCell{i};
    end
end
 
function [dE_dPhi] = gradientCountourPlot(energyMatrix, phi_fraction, theta_normalized, delta_range, k_bend, lambda)

    % Already computed:
    dTheta = 1;
    dPhi = mean(diff(phi_fraction(1,:)));
    [~, dE_dPhi] = gradient(energyMatrix, dTheta, dPhi);
    [CP, x, y, pp, xspline, yspline] = fitSplineToCriticalPoints(delta_range, k_bend, lambda);

    % === Determine symmetric color axis for consistent zero color ===
    maxAbsGrad = max(abs(dE_dPhi(:)));
    clim = [-maxAbsGrad, +maxAbsGrad];
    
    % Plot gradient contour
    figure(1);
    contourf(theta_normalized, phi_fraction, dE_dPhi, 30, 'LineColor', 'none');

    % Set diverging colormap (blue ↔ white ↔ red)
    colormap(redWhiteBlue());
    caxis(clim);  % Fix color axis so 0 is always centered
    colorbar;

    xlabel('\theta / \pi (Normalized Contact Angle)');
    ylabel('\phi / \theta (Fractional Deflection)');
    zlabel('Free Energy Gradient');
    title('Contour of \partial E / \partial \phi (Gradient of Energy Landscape)');
    axis tight;
    hold on;

    % === Plot critical points ===
    % critPoints = getCriticalPoints(delta_range, k_bend, lambda);
    % for i = 1:length(critPoints)
    %     theta = critPoints(i).theta;
    %     phi_crit = critPoints(i).phi_crit;
    %     is_stable = critPoints(i).is_stable;
    % 
    %     if isempty(phi_crit)
    %         continue;
    %     end
    % 
    %     if isscalar(theta)
    %         theta = repmat(theta, size(phi_crit));
    %     end
    % 
    %     theta_norm = theta / pi;
    %     phi_frac = phi_crit ./ theta;
    % 
    %     % Plot stable points (minima) as filled circles
    %     plot(theta_norm(is_stable), phi_frac(is_stable), ...
    %         'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 6);
    % 
    %     % Plot unstable points (maxima) as empty circles
    %     plot(theta_norm(~is_stable), phi_frac(~is_stable), ...
    %         'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 6);
    %         %'ko', 'MarkerSize', 6);
    % end
    % hold on

    % == Plot spline == 
    plot(xspline, yspline, 'b-', 'LineWidth',1.4);


    hold on
    
end

function [dE_dPhi] = gradientCountourPlotNoFitting(energyMatrix, phi_fraction, theta_normalized, delta_range, k_bend, lambda)

    % Already computed:
    dTheta = 1;
    dPhi = mean(diff(phi_fraction(1,:)));
    [~, dE_dPhi] = gradient(energyMatrix, dTheta, dPhi);

    % === Determine symmetric color axis for consistent zero color ===
    maxAbsGrad = max(abs(dE_dPhi(:)));
    clim = [-maxAbsGrad, +maxAbsGrad];
    
    % Plot gradient contour
    figure(1);
    contourf(theta_normalized, phi_fraction, dE_dPhi, 30, 'LineColor', 'none');

    % Set diverging colormap (blue ↔ white ↔ red)
    colormap(redWhiteBlue());
    caxis(clim);  % Fix color axis so 0 is always centered
    colorbar;

    xlabel('\theta / \pi (Normalized Contact Angle)');
    ylabel('\phi / \theta (Fractional Deflection)');
    zlabel('Free Energy Gradient');
    title('Contour of \partial E / \partial \phi (Gradient of Energy Landscape)');
    axis tight;
    hold on;

     %=== Plot critical points ===
    critPoints = getCriticalPoints(delta_range, k_bend, lambda);
    for i = 1:length(critPoints)
        theta = critPoints(i).theta;
        phi_crit = critPoints(i).phi_crit;
        is_stable = critPoints(i).is_stable;

        if isempty(phi_crit)
            continue;
        end

        if isscalar(theta)
            theta = repmat(theta, size(phi_crit));
        end

        theta_norm = theta / pi;
        phi_frac = phi_crit ./ theta;

        % Plot stable points (minima) as filled circles
        plot(theta_norm(is_stable), phi_frac(is_stable), ...
            'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 6);

        % Plot unstable points (maxima) as empty circles
        plot(theta_norm(~is_stable), phi_frac(~is_stable), ...
            'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 6);
            %'ko', 'MarkerSize', 6);
    end
    hold on

end

function [criticalPoints] = getCriticalPoints(delta_range, k_bend, lambda)
% Extended to include stability classification (min vs max)

    criticalPoints = struct('delta', {}, 'phi_crit', {}, 'E_crit', {}, ...
                            'theta', {}, 'is_stable', {});
    
    for i = 1:length(delta_range)
        delta = delta_range(i);
        
        [E, phi, theta] = energyForGivenContactAngle(delta, k_bend, lambda);
        
        dphi = diff(phi);
        dE = diff(E);
        phi_mid = (phi(1:end-1) + phi(2:end)) / 2;
        dE_dphi = dE ./ dphi;
    
        sign_change = dE_dphi(1:end-1) .* dE_dphi(2:end) < 0;
        idx_zero_cross = find(sign_change);
        
        phi_crit = zeros(size(idx_zero_cross));
        E_crit = zeros(size(idx_zero_cross));
        is_stable = false(size(idx_zero_cross));
    
        for j = 1:length(idx_zero_cross)
            idx = idx_zero_cross(j);
    
            % Interpolate critical phi (zero crossing)
            x1 = phi_mid(idx); x2 = phi_mid(idx+1);
            y1 = dE_dphi(idx); y2 = dE_dphi(idx+1);
            phi_zero = x1 - y1 * (x2 - x1) / (y2 - y1);
            phi_crit(j) = phi_zero;
    
            % Interpolated energy at critical point
            E_crit(j) = interp1(phi, E, phi_zero, 'spline');
    
            % Estimate second derivative: central difference on gradient
            if idx > 1 && idx < length(dE_dphi)-1
                dd = (dE_dphi(idx+1) - dE_dphi(idx-1)) / (phi_mid(idx+1) - phi_mid(idx-1));
            else
                dd = NaN;
            end
            is_stable(j) = dd > 0;
        end
    
        criticalPoints(i).delta = delta;
        criticalPoints(i).phi_crit = phi_crit;
        criticalPoints(i).E_crit = E_crit;
        criticalPoints(i).theta = theta;
        criticalPoints(i).is_stable = is_stable;
    end
end

function [critPoints, x, y, pp, xq, yq] = fitSplineToCriticalPoints(delta_range, k_bend, lambda)

    [critPoints] = getCriticalPoints(delta_range, k_bend, lambda);

    x = [];
    y = [];

    for i = 1:length(critPoints)

        theta = critPoints(i).theta;
        phi_crit = critPoints(i).phi_crit;
        is_stable = critPoints(i).is_stable;
    
        if isempty(phi_crit)
            continue;
        end
    
        if isscalar(theta)
            theta = repmat(theta, size(phi_crit));
        end
    
        x = [x, theta / pi];
        y = [y, (phi_crit ./ theta)];

    end

    pp = spline(x, y);

    % Evaluate on fine grid
    xq = linspace(min(x), max(x), 200);
    yq = ppval(pp, xq);

    % % Plot
    % figure
    % plot(x, y, 'ko', 'MarkerFaceColor','k'); 
    % hold on
    % plot(xq, yq, 'b-', 'LineWidth',1.4);
    % xlabel('\theta/\pi');
    % ylabel('\phi/\theta');
    % title('Cubic Spline Fit to Critical Points');

end

function [k_fit, lambda_fit] = experimentalSplineFit(expTheta, expPhi, delta_range, initial_guess)

    lb = [0.0001, -10];
    ub = [20, 10];

    % --- ERROR FUNCTION ---
    % This function will be minimized
    errFun = @(params) splineResiduals(params, delta_range, expTheta, expPhi);
    
    % --- OPTIMIZATION OPTIONS ---
    opts = optimoptions('lsqnonlin', ...
        'Display','iter', ...
        'MaxFunctionEvaluations', 1000, ...
        'FunctionTolerance', 1e-6);
    
    % --- SOLVE FOR k_bend ---
    [fit_params, ~] = lsqnonlin(errFun, initial_guess, lb, ub, opts);

    % Unpack results
    k_fit = fit_params(1);
    lambda_fit = fit_params(2);
    
    fprintf('Best-fit k_bend = %.4f, lambda = %.4f\n', k_fit, lambda_fit);
   
end

function [residuals] = splineResiduals(params, delta_range, expTheta, expPhi)

    % Unpack parameter vector
    k_bend = params(1);
    lambda  = params(2);

    [critPoints, x, y, pp, xq, yq] = fitSplineToCriticalPoints(delta_range, k_bend, lambda);

    y_pred = ppval(pp, expTheta);

    residuals = y_pred - expPhi;
    
end

function crit = getCriticalPointsOverLambda(delta_range, k_bend, lambda_range)
%GETCRITICALPOINTSOVERLAMBDA  Critical points & stability over delta and lambda.
%
%   crit = GETCRITICALPOINTSOVERLAMBDA(delta_range, k_bend, lambda_range)
%   uses energyForGivenContactAngle(delta, k_bend, lambda), assumed to return:
%       [E, phi, theta]
%
%   Output 'crit' is a struct array with fields:
%       delta      - intrinsic contact angle parameter from delta_range
%       lambda     - lambda value
%       phi_crit   - phi at critical point (where dE/dphi = 0)
%       theta_crit - theta at the critical point (or scalar theta if no phi-dependence)
%       E_crit     - energy at critical point
%       is_stable  - logical, true if d2E/dphi2 > 0 (local minimum)

    crit = struct('delta', {}, 'lambda', {}, ...
                  'phi_crit', {}, 'theta_crit', {}, ...
                  'E_crit', {}, 'is_stable', {});

    next = 1;  % index for filling struct array

    for iDelta = 1:numel(delta_range)
        delta = delta_range(iDelta);

        for iLam = 1:numel(lambda_range)
            lambda = lambda_range(iLam);

            % --- 1) Call your existing energy function
            [E, phi, theta] = energyForGivenContactAngle(delta, k_bend, lambda);

            % Ensure column vectors for phi and E
            phi = phi(:);
            E   = E(:);

            % --- 2) First derivative dE/dphi using finite differences
            dphi    = diff(phi);
            dE      = diff(E);
            phi_mid = (phi(1:end-1) + phi(2:end)) / 2;
            dE_dphi = dE ./ dphi;

            % --- 3) Find zero crossings of dE/dphi (critical points)
            sign_change    = dE_dphi(1:end-1) .* dE_dphi(2:end) < 0;
            idx_zero_cross = find(sign_change);

            if isempty(idx_zero_cross)
                continue
            end

            for j = 1:numel(idx_zero_cross)
                idx = idx_zero_cross(j);

                % --- 3a) Interpolate phi where dE/dphi = 0 (linear in gradient)
                x1 = phi_mid(idx);     x2 = phi_mid(idx+1);
                y1 = dE_dphi(idx);     y2 = dE_dphi(idx+1);

                phi_zero = x1 - y1*(x2 - x1)/(y2 - y1);
                phi_crit = phi_zero;

                % --- 3b) Interpolate energy at phi_crit
                E_crit = interp1(phi, E, phi_zero, 'spline');

                % --- 3c) Handle theta (may be scalar or vector)
                if isscalar(theta)
                    % theta is independent of phi
                    theta_crit = theta;
                elseif numel(theta) == numel(phi)
                    theta = theta(:);
                    theta_crit = interp1(phi, theta, phi_zero, 'spline');
                else
                    % Size mismatch → cannot safely interpolate
                    theta_crit = NaN;
                end

                % --- 4) Estimate second derivative (no Curve Fitting Toolbox)
                if idx > 1 && idx < numel(dE_dphi)-1
                    % Central difference of gradient (on midpoints)
                    num2 = dE_dphi(idx+1) - dE_dphi(idx-1);
                    den2 = phi_mid(idx+1)  - phi_mid(idx-1);
                    d2E  = num2 / den2;
                elseif idx > 1 && idx < numel(phi)-1
                    % Fallback: second difference on E(phi) at grid points
                    h = phi(idx+1) - phi(idx);    % assume roughly uniform
                    d2E = (E(idx+1) - 2*E(idx) + E(idx-1)) / (h^2);
                else
                    % At boundaries we don't have enough neighbors
                    d2E = NaN;
                end

                is_stable = (d2E > 0);

                % --- 5) Store in struct array
                crit(next).delta      = delta;
                crit(next).lambda     = lambda;
                crit(next).phi_crit   = phi_crit;
                crit(next).theta_crit = theta_crit;
                crit(next).E_crit     = E_crit;
                crit(next).is_stable  = is_stable;

                next = next + 1;
            end
        end
    end
end

function plotStabilityTransition(crit)
%PLOTSTABILITYTRANSITION  Visualize stability vs lambda.
%
%   plotStabilityTransition(crit)
%   where crit is the struct array from getCriticalPointsOverLambda.

    if isempty(crit)
        warning('plotStabilityTransition:Empty', ...
                'No critical points in input struct.');
        return
    end

    theta = [crit.theta_crit];
    lambda = [crit.lambda];
    isStable = [crit.is_stable];

    figure;
    hold on;

    % Stable: blue filled circles
    plot(theta(isStable)/pi, lambda(isStable), ...
        'bo', 'MarkerFaceColor','b', 'DisplayName','Stable');

    % Unstable: red crosses
    plot(theta(~isStable)/pi, lambda(~isStable), ...
        'rx', 'LineWidth',1.5, 'DisplayName','Unstable');

    xlabel('\theta_c / \pi');
    ylabel('\lambda');
    legend('Location','best');
    grid on;
    box on;
end

function plotPhaseBoundary(crit, nbins)
%PLOTPHASEBOUNDARY  Overlay a stability (phase) boundary on your scatter plot.
%
%   plotPhaseBoundary(crit, nbins)
%   crit  : struct array from getCriticalPointsOverLambda
%   nbins : number of theta bins (default 60)

    if nargin < 2
        nbins = 60;
    end

    T = struct2table(crit);

    theta = T.theta_crit / pi;   % normalized x-axis
    lambda = T.lambda;
    isStable = T.is_stable;

    % ---- 1) Make your scatter plot (if not already plotted)
    figure; hold on;
    plot(theta(isStable),   lambda(isStable),   'bo', ...
        'MarkerFaceColor','b', 'DisplayName','Stable');
    plot(theta(~isStable),  lambda(~isStable),  'rx', ...
        'LineWidth',1.0, 'DisplayName','Unstable');
    xlabel('\theta_c / \pi');
    ylabel('\lambda');
    grid on; box on;
    legend('Location','southwest');

    % ---- 2) Bin in theta and find boundary in each bin
    % theta edges from min to max:
    edges = linspace(min(theta), max(theta), nbins+1);
    centers = 0.5*(edges(1:end-1) + edges(2:end));

    lamBoundary = NaN(size(centers));

    for k = 1:numel(centers)
        inBin = theta >= edges(k) & theta < edges(k+1);

        lamStable = lambda(inBin & isStable);

        if ~isempty(lamStable)
            % upper envelope of stable region
            lamBoundary(k) = max(lamStable);
        end
    end

    % Remove empty bins
    mask = ~isnan(lamBoundary);

    % ---- 3) Plot the phase boundary
    plot(centers(mask), lamBoundary(mask), 'k-', 'LineWidth',2, ...
        'DisplayName','Phase boundary');

    legend('Location','southwest');
end

function plotPhaseBoundariesPoly(crit, nbins, degStable, degUnstable)

%PLOTPHASEBOUNDARIESPOLY  Smooth polynomial boundaries for stable & unstable regions.
%
%   plotPhaseBoundariesPoly(crit)
%   plotPhaseBoundariesPoly(crit, nbins)
%   plotPhaseBoundariesPoly(crit, nbins, degStable, degUnstable)
%
%   crit       : struct array from getCriticalPointsOverLambda
%   nbins      : number of theta bins for envelope construction (default 150)
%   degStable  : poly degree for stable boundary lambda(theta) (default 5)
%   degUnstable: poly degree for unstable boundary lambda(theta) (default 5)

    if nargin < 2 || isempty(nbins)
        nbins = 150;
    end
    if nargin < 3 || isempty(degStable)
        degStable = 5;
    end
    if nargin < 4 || isempty(degUnstable)
        degUnstable = 5;
    end

    T = struct2table(crit);

    theta    = T.theta_crit / pi;  % normalized x-axis
    lambda   = T.lambda;
    isStable = T.is_stable;

    %--------------------------------------------------------------
    % 1) Scatter plot
    %--------------------------------------------------------------
    figure; hold on;

    plot(theta(isStable),  lambda(isStable),  'bo', ...
        'MarkerFaceColor','b', 'DisplayName','Stable');
    plot(theta(~isStable), lambda(~isStable), 'rx', ...
        'LineWidth',1.0, 'DisplayName','Unstable');

    xlabel('\theta_c / \pi');
    ylabel('\lambda');
    grid on; box on;

    %--------------------------------------------------------------
    % 2) Bin in theta
    %--------------------------------------------------------------
    edges   = linspace(min(theta), max(theta), nbins+1);
    centers = 0.5*(edges(1:end-1) + edges(2:end));

    lamStableMax   = NaN(size(centers));   % upper envelope of stable region
    lamUnstableMin = NaN(size(centers));   % lower envelope of unstable region

    for k = 1:numel(centers)
        inBinStable   = theta >= edges(k) & theta < edges(k+1) & isStable;
        inBinUnstable = theta >= edges(k) & theta < edges(k+1) & ~isStable;

        if any(inBinStable)
            lamStableMax(k) = max(lambda(inBinStable));
        end
        if any(inBinUnstable)
            lamUnstableMin(k) = min(lambda(inBinUnstable));
        end
    end

    %--------------------------------------------------------------
    % 3) Stable boundary: clean, fit, and plot
    %--------------------------------------------------------------
    ms = ~isnan(lamStableMax);
    theta_s = centers(ms).';
    lam_s   = lamStableMax(ms).';

    [theta_s, idxS] = sort(theta_s);
    lam_s = lam_s(idxS);

    % raw points (optional, dotted)
    plot(theta_s, lam_s, 'k:', 'LineWidth',1.0, ...
        'DisplayName','Stable boundary (raw)');

    % polynomial fit & smooth curve
    pS = polyfit(theta_s, lam_s, degStable);
    theta_s_fit = linspace(min(theta_s), max(theta_s), 400);
    lam_s_fit   = polyval(pS, theta_s_fit);

    plot(theta_s_fit, lam_s_fit, 'k-', 'LineWidth',2.0, ...
        'DisplayName', sprintf('Stable boundary (deg %d)', degStable));

    %--------------------------------------------------------------
    % 4) Unstable boundary: clean, fit, and plot
    %--------------------------------------------------------------
    mu = ~isnan(lamUnstableMin);
    theta_u = centers(mu).';
    lam_u   = lamUnstableMin(mu).';

    [theta_u, idxU] = sort(theta_u);
    lam_u = lam_u(idxU);

    % raw unstable boundary (optional, dotted)
    plot(theta_u, lam_u, 'g:', 'LineWidth',1.0, ...
        'DisplayName','Unstable boundary (raw)');

    if numel(theta_u) > degUnstable
        pU = polyfit(theta_u, lam_u, degUnstable);
        theta_u_fit = linspace(min(theta_u), max(theta_u), 400);
        lam_u_fit   = polyval(pU, theta_u_fit);

        plot(theta_u_fit, lam_u_fit, 'g-', 'LineWidth',2.0, ...
            'DisplayName', sprintf('Unstable boundary (deg %d)', degUnstable));
    else
        warning('Not enough unstable boundary points to fit a degree-%d polynomial.', ...
                degUnstable);
    end

    legend('Location','southwest');
end

function cmap = redWhiteBlue(n)

    if nargin < 1
        n = 256;  % default colormap size
    end
    half = floor(n/2);
    
    % Blue to white
    blue = [linspace(0,1,half)', linspace(0,1,half)', ones(half,1)];
    
    % White to red
    red = [ones(half,1), linspace(1,0,half)', linspace(1,0,half)'];
    
    cmap = [blue; red];

end

function [R_prime, R_1_prime, phi, theta_intrinsic] = plotCapGeometries(theta_intrinsic)

    [R_prime, phi, theta_intrinsic, theta_intrinsic_vector] = doubleSphericalCapConstraint(theta_intrinsic);

    % compute R'_prime to get bottom half of spherical cap
    R_1_prime = (R_prime.*sin(phi))./(sin(theta_intrinsic - phi));

    plotDoubleCapConnected2D(R_prime,R_1_prime,phi,theta_intrinsic);

end

function plotDoubleCapConnected2D(R_top, R_bot, theta_above, theta_intrinsic, opts)
% plotDoubleCapConnected2D  Plot 2D "double spherical cap" (lens) shapes.
%
% Inputs (vectors must be same length):
%   R_top           top sphere radii (can be normalized)
%   R_bot           bottom sphere radii (can be normalized)
%   theta_above     angle ABOVE the baseline, measured from vertical radius (rad)
%   theta_intrinsic intrinsic contact angle theta (rad), material property
%
% Geometry:
%   theta_below = theta_intrinsic - theta_above
%   a_top = R_top*sin(theta_above)
%   a_bot = R_bot*sin(theta_below)
%   Use a = 0.5*(a_top + a_bot) (they should match if computed consistently)
%
% Arc formulas (baseline is y=0):
%   y_top(x) = -c_top + sqrt(R_top^2 - x^2),   c_top = sqrt(R_top^2 - a^2)
%   y_bot(x) =  c_bot - sqrt(R_bot^2 - x^2),   c_bot = sqrt(R_bot^2 - a^2)

    arguments
        R_top (:,1) double
        R_bot (:,1) double
        theta_above (:,1) double
        theta_intrinsic (1,1) double
        opts.NumPoints (1,1) double {mustBeInteger,mustBePositive} = 400
        opts.LineWidth (1,1) double {mustBePositive} = 1.5
        opts.ShowBaseLine (1,1) logical = true
        opts.ColorByIndex (1,1) logical = true
    end

    n = numel(R_top);
    assert(numel(R_bot)==n && numel(theta_above)==n, ...
        'R_top, R_bot, theta_above must have the same length.');

    theta_below = theta_intrinsic - theta_above;

    figure
    ax = axes; hold(ax,'on'); grid(ax,'on')
    axis(ax,'equal')
    xlabel(ax,'x'); ylabel(ax,'y')
    title(ax,'Double spherical caps connected at endpoints')
    
    N = 800;                        % more points helps for near-vertical sides
    cmap = turbo(numel(R_top));     % optional nice coloring :contentReference[oaicite:1]{index=1}
    
    for i = 1:numel(R_top)
        Rt = R_top(i);
        Rb = R_bot(i);
        aA = theta_above(i);
        bB = theta_below(i);
    
        if any(~isfinite([Rt Rb aA bB])) || any([Rt Rb] <= 0) || any([aA bB] <= 0)
            continue
        end
    
        % Base half-width implied by each cap (should match if your model is consistent)
        abase_top = Rt*sin(aA);
        abase_bot = Rb*sin(bB);
    
        % If they don't match numerically, pick a policy:
        % (1) enforce shared base exactly by adjusting bB:
        % bB = asin(min(max(abase_top/Rb, -1), 1));
        % abase = abase_top;
        %
        % (2) or average (keeps your precomputed angles):
        abase = 0.5*(abase_top + abase_bot);
    
        % Skip impossible geometry
        if abase <= 0 || abase > min(Rt,Rb)
            continue
        end
    
        % --- TOP ARC (center at y = -Rt*cos(aA), parameter t in [-aA, +aA]) ---
        yt_center = -Rt*cos(aA);
        tTop = linspace(-aA, +aA, N);
        xTop = Rt*sin(tTop);
        yTop = yt_center + Rt*cos(tTop);
    
        % --- BOTTOM ARC (center at y = +Rb*cos(bB), parameter t in [pi-bB, pi+bB]) ---
        yb_center =  Rb*cos(bB);
        tBot = linspace(pi-bB, pi+bB, N);
        xBot = Rb*sin(tBot);
        yBot = yb_center + Rb*cos(tBot);
    
        col = cmap(i,:);
        plot(ax, xTop, yTop, 'Color', col, 'LineWidth', 1.6)
        plot(ax, xBot, yBot, 'Color', col, 'LineWidth', 1.6)
    end
    
    % baseline
    xl = xlim(ax);
    plot(ax, xl, [0 0], 'k-', 'LineWidth', 1)
    hold(ax,'off')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set each if block to true to execute anlysis, data fitting, or plot
% generation. Each block will be described here (in the code) and in the
% Physical_Modeling README.txt file.

% Code used to produce Extended Data Fig. 7 
% Plots a given contact angle for different above and below membrane
% deflection angles

if false == true
    [R_prime, R_1_prime, phi, theta] = plotCapGeometries(0.75);
end

% Code used to produce Extended Data Fig. 8
% Calculates the energy landscape and energy landscape gradient for
% different normalized k_bend and lambda values for all possible intrinsic
% contact angles (-pi, pi) and above membrane deflection angles (0, 1). 

if false == true

    delta_range = linspace(-0.98, 0.98, 40);
    k_bend = 0.5;
    lambda = -0.5;

    [energyMatrix, phiMatrix, thetaMatrix] = energyLandscape(delta_range, k_bend, lambda);

    % Normalize axes
    theta_normalized = thetaMatrix / pi;
    phi_fraction = phiMatrix ./ thetaMatrix;
    
    % Plot surface
    figure(3);
    surf(theta_normalized, phi_fraction, energyMatrix);
    xlabel('\theta / \pi (Normalized Contact Angle)');
    ylabel('\phi / \theta (Fractional Deflection)');
    zlabel('Free Energy');
    title('Normalized Free Energy Landscape');
    shading interp;
    view(135, 30);
    colormap turbo;
    
    [dE_dPhi] = gradientCountourPlotNoFitting(energyMatrix,phi_fraction,theta_normalized, delta_range, k_bend, lambda);
    
end

% Code used to produce Fig. 4.
% Assuming equilibrium solutions, the CA and DA data point are fit to
% spline. k_bend and lambda fits are extracted using lsqnonlin fit.
% K_bend is assumed to be possitve.
% lambda is assumed to be negative. 

if false == true

    % Uncommment the following data vectors to reproduce gradient spline
    % fitting. The data corresponds to the measured contact angles (CA) and the
    % measured above membrane deflection angles (DA).

    % 0.1nM vgRNA
    CA = [140 130 137 133 142 114 126 136 135 114 124 106 113 125 135 141 118 139 131 133 125];
    DA = [96 67 51 55 71 49 79 65 87 64 64 29 59 88 74 101 40 88 50 77 37];

    % 0.5nM F3RNA
    % CA = [132 91 93 131 117 102 122 114 112];
    % DA = [112 65 41 71 49 40 88 85 89];
    
    % no RNA
    % CA = [47 47 40 51 47];
    % DA = [26 40 40 51 47];
    
    CA_norm = CA/180; % fraction of pi
    DA_norm = (DA./CA); % fraction of pi
    
    delta_range = linspace(-0.98, 0.98, 40); % provide theta range using the Young-Dupres relation - relate the difference in surface tension to an intrinsic contact angle. The difference in surface tensions is also an important term in the surface free energy
    initial_guess = [0.1, -1];
    
    [k_fit, lambda_fit] = experimentalSplineFit(CA_norm, DA_norm, delta_range, initial_guess);
    
    [energyMatrix, phiMatrix, thetaMatrix] = energyLandscape(delta_range, k_fit, lambda_fit);

    % Normalize axes
    theta_normalized = thetaMatrix / pi;
    phi_fraction = phiMatrix ./ thetaMatrix;
    
    % Plot surface
    figure(3);
    surf(theta_normalized, phi_fraction, energyMatrix);
    xlabel('\theta / \pi (Normalized Contact Angle)');
    ylabel('\phi / \theta (Fractional Deflection)');
    zlabel('Free Energy');
    title('Normalized Free Energy Landscape');
    shading interp;
    view(135, 30);
    colormap turbo;
    
    [dE_dPhi] = gradientCountourPlot(energyMatrix,phi_fraction,theta_normalized, delta_range, k_fit, lambda_fit);
    
    figure(1);
    scatter(CA_norm, DA_norm, "filled");

end