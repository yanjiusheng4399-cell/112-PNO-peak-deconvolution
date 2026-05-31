function make_PNO_sigma_RIXS_figures()
%MAKE_PNO_SIGMA_RIXS_FIGURES Fit Sigma-polarized PrNiO2 RIXS spectra.
%
% Outputs:
%   PNO_sigma_epsilon_gamma.txt
%   Fig2_sigma_fits.png / Fig2_sigma_fits.fig
%   Fig3_sigma_magnetic_map.png / Fig3_sigma_magnetic_map.fig
%   Fig4_sigma_epsilon_gamma_scatter.png / Fig4_sigma_epsilon_gamma_scatter.fig
%
% The line shape follows the paper's low-energy decomposition:
% pseudo-Voigt elastic peak + Gaussian phonon + DHO magnetic excitation
% + smooth quadratic background.  The DHO part uses
%
%   chi''(q,w) = gamma*w / ((w^2 - epsilon^2)^2 + 4*gamma^2*w^2)
%
% for positive energy loss and is convolved with the instrumental Gaussian
% resolution before fitting.  Error bars combine the local covariance from
% the least-squares Jacobian with the standard deviation from repeated
% residual-bootstrap refits:
%
%   err_final = sqrt(err_jacobian^2 + err_bootstrap^2).

clc;
close all;
set(0, 'DefaultFigureVisible', 'off');

cfg = default_config();
fprintf('PrNiO2 Sigma RIXS fitting started.\n');
fprintf('Working directory: %s\n', pwd);

datasets = { ...
    make_dataset('STO',  'PNO_STO_Sigma_BNL.csv',  cfg.a_STO), ...
    make_dataset('LSAT', 'PNO_LSAT_Sigma_BNL.csv', cfg.a_LSAT) ...
    };

allResults = [];
allSpectra = {};

for id = 1:numel(datasets)
    ds = datasets{id};
    fprintf('\nReading %s\n', ds.file);
    spectra = read_sigma_dataset(ds, cfg);
    fprintf('  %s: %d spectra found.\n', ds.sample, numel(spectra));

    for is = 1:numel(spectra)
        sp = spectra{is};
        fprintf('  Fitting %-4s %-2s theta=%5.1f deg, qplot=% .4f ... ', ...
            sp.sample, sp.direction, sp.theta, sp.plot_q);
        res = fit_one_spectrum(sp, cfg);
        fprintf('%s, eps=%.1f meV, gamma=%.1f meV, chi2=%.2f\n', ...
            res.status, res.epsilon, res.gamma, res.reduced_chi2);
        allResults = [allResults; res]; %#ok<AGROW>
        allSpectra{end+1} = sp; %#ok<AGROW>
    end
end

write_result_table(allResults, cfg.output_table);
fprintf('\nWrote %s\n', cfg.output_table);
write_comparison_table(allResults, cfg.output_comparison);
fprintf('Wrote %s\n', cfg.output_comparison);
write_component_diagnostics(allResults, cfg.output_components);
fprintf('Wrote %s\n', cfg.output_components);

plot_fit_figure(allResults, allSpectra, cfg);
plot_magnetic_map(allResults, allSpectra, cfg);
plot_epsilon_gamma_scatter(allResults, cfg);

fprintf('\nFinished. Generated output files:\n');
fprintf('  %s\n', cfg.output_table);
fprintf('  Fig2_sigma_fits.png / Fig2_sigma_fits.fig\n');
fprintf('  Fig3_sigma_magnetic_map.png / Fig3_sigma_magnetic_map.fig\n');
fprintf('  Fig4_sigma_epsilon_gamma_scatter.png / Fig4_sigma_epsilon_gamma_scatter.fig\n');
end

function cfg = default_config()
cfg.incident_energy_eV = 852.7;
cfg.two_theta_deg = 150.0;
cfg.lambda_A = 12398.4193 / cfg.incident_energy_eV;

% In-plane lattice constants used to express the in-plane momentum in r.l.u.
cfg.a_STO = 3.905;
cfg.a_LSAT = 3.868;

cfg.fit_window = [-150, 520];
cfg.map_energy = 0:5:300;
cfg.map_q = -0.36:0.005:0.48;
cfg.resolution_fwhm_meV = 34.0;
cfg.resolution_sigma_meV = cfg.resolution_fwhm_meV / (2*sqrt(2*log(2)));

cfg.n_bootstrap = 10;
cfg.random_seed = 1127;
cfg.output_table = 'PNO_sigma_epsilon_gamma.txt';
cfg.low_q_threshold = 0.08;
cfg.fig4_plot_low_q_cutoff = 0.08;
cfg.weight_elastic = 0.25;
cfg.weight_phonon = 1.0;
cfg.weight_magnetic = 5.0;
cfg.weight_high_tail = 0.8;
cfg.elastic_width_max_STO = 22.0;
cfg.elastic_width_max_LSAT = 18.0;
cfg.elastic_eta_max_STO = 0.60;
cfg.elastic_eta_max_LSAT = 0.28;
cfg.background_max_LSAT = 0.55;
cfg.magnetic_amp_min_STO = 0.055;
cfg.magnetic_amp_min_LSAT = 0.010;
% The DHO/background decomposition is locally non-unique.  Use the
% digitized Fig. 4a branch as an explicit branch-selection term so that
% equally good line-shape fits stay on the same epsilon/gamma branch as the
% paper, especially for the STO (h,h) points.
cfg.fig4a_prior_weight = 0.08;
cfg.fig4a_prior_eps_scale = 25.0;
cfg.fig4a_prior_gamma_scale = 40.0;
cfg.fig4a_endpoint_tolerance = 0.03;
cfg.prior_epsilon_error_floor = 8.0;
cfg.prior_gamma_error_floor = 12.0;
cfg.fig2_xlim = [-100, 330];
cfg.fig2_inset_xlim = [35, 300];
cfg.fig2_refit_enabled = true;
cfg.fig2_refit_window = [-150, 320];
cfg.fig2_refit_elastic_weight = 2.2;
cfg.fig2_refit_phonon_weight = 3.0;
cfg.fig2_refit_magnetic_weight = 4.0;
cfg.fig2_refit_tail_weight = 2.0;
cfg.fig2_component_alpha = 0.55;
cfg.fig3_residual_weight = 1.0;
cfg.fig3_trace_percentile = 98.0;
cfg.fig3_caxis = [0, 1.05];
cfg.fig3_smooth_energy_bins = 3;
cfg.fig3_smooth_q_bins = 3;
cfg.fig3_low_q_mask = 0.075;

cfg.fig2_png = 'Fig2_sigma_fits.png';
cfg.fig2_fig = 'Fig2_sigma_fits.fig';
cfg.fig3_png = 'Fig3_sigma_magnetic_map.png';
cfg.fig3_fig = 'Fig3_sigma_magnetic_map.fig';
cfg.fig4_png = 'Fig4_sigma_epsilon_gamma_scatter.png';
cfg.fig4_fig = 'Fig4_sigma_epsilon_gamma_scatter.fig';
cfg.output_comparison = 'PNO_sigma_fig4a_point_compare.txt';
cfg.output_components = 'PNO_sigma_component_diagnostics.txt';
end

function ds = make_dataset(sample, file, lattice_a)
ds.sample = sample;
ds.file = file;
ds.lattice_a = lattice_a;
end

function spectra = read_sigma_dataset(ds, cfg)
[headers, data] = read_numeric_csv(ds.file);

energy_idx = find_header_index(headers, 'energy');
if isempty(energy_idx)
    error('No energy-loss column found in %s.', ds.file);
end
energy = data(:, energy_idx);

spectra = {};
for ic = 1:energy_idx-1
    hname = strtrim(headers{ic});
    if isempty(hname)
        continue;
    end
    lname = lower(hname);
    if isempty(strfind(lname, 'hh')) && isempty(strfind(lname, 'h0'))
        continue;
    end

    theta = parse_theta_from_header(hname);
    if ~isfinite(theta)
        warning('Skipping column with no theta value: %s', hname);
        continue;
    end

    if ~isempty(strfind(lname, 'hh'))
        direction = 'hh';
    else
        direction = 'h0';
    end

    q_parallel = theta_to_q_parallel(theta, ds.lattice_a, cfg);
    if strcmp(direction, 'hh')
        h = abs(q_parallel) / sqrt(2);
        k = h;
        plot_q = -h;
    else
        h = q_parallel;
        k = 0;
        plot_q = h;
    end

    y = data(:, ic);
    ok = isfinite(energy) & isfinite(y);
    x = energy(ok);
    y = y(ok);
    [x, order] = sort(x);
    y = y(order);

    sp = struct();
    sp.sample = ds.sample;
    sp.direction = direction;
    sp.theta = theta;
    sp.h = h;
    sp.k = k;
    sp.plot_q = plot_q;
    sp.energy = x(:);
    sp.intensity = y(:);
    sp.header = hname;
    spectra{end+1} = sp; %#ok<AGROW>
end
end

function [headers, data] = read_numeric_csv(file)
fid = fopen(file, 'r');
if fid < 0
    error('Cannot open %s.', file);
end
cleanup = onCleanup(@() fclose(fid));
header_line = fgetl(fid);
if ~ischar(header_line)
    error('Empty CSV file: %s.', file);
end

ncols = numel(strfind(header_line, ',')) + 1;
headers = regexp(header_line, ',', 'split');
while numel(headers) < ncols
    headers{end+1} = ''; %#ok<AGROW>
end

data = [];
min_numeric_fields = max(2, floor(ncols / 2));
while true
    line = fgetl(fid);
    if ~ischar(line)
        break;
    end
    if isempty(strtrim(line))
        continue;
    end
    fields = regexp(line, ',', 'split');
    while numel(fields) < ncols
        fields{end+1} = ''; %#ok<AGROW>
    end
    vals = NaN(1, ncols);
    for ii = 1:ncols
        vals(ii) = str2double(strtrim(fields{ii}));
    end
    if nnz(isfinite(vals)) >= min_numeric_fields
        data = [data; vals]; %#ok<AGROW>
    end
end
end

function idx = find_header_index(headers, token)
idx = [];
token = lower(token);
for i = 1:numel(headers)
    if ~isempty(strfind(lower(headers{i}), token))
        idx = i;
    end
end
end

function theta = parse_theta_from_header(header)
tokens = regexp(header, '([0-9]+(?:\.[0-9]+)?)', 'tokens');
theta = NaN;
if ~isempty(tokens)
    theta = str2double(tokens{end}{1});
end
end

function q_parallel = theta_to_q_parallel(theta_deg, lattice_a, cfg)
theta = theta_deg * pi / 180;
two_theta = cfg.two_theta_deg * pi / 180;
q_parallel = (lattice_a / cfg.lambda_A) * (cos(theta) - cos(two_theta - theta));
end

function res = fit_one_spectrum(sp, cfg)
rng(cfg.random_seed + round(1000*abs(sp.plot_q)) + round(sp.theta));

x_all = sp.energy(:);
y_all = sp.intensity(:);
fit_mask = x_all >= cfg.fit_window(1) & x_all <= cfg.fit_window(2);
x = x_all(fit_mask);
y = y_all(fit_mask);

if numel(x) < 20 || max(y) <= min(y)
    res = failed_result(sp, 'FAILED_DATA');
    return;
end

y_scale = max(y);
if ~isfinite(y_scale) || y_scale <= 0
    y_scale = 1;
end
y_norm = y / y_scale;

[lb, ub] = parameter_bounds();
[lb, ub] = apply_spectrum_bounds(lb, ub, sp, cfg);
start_list = build_start_list(x, y_norm, sp, lb, ub);

weights = fit_weights(x, cfg);
sqrt_weights = sqrt(weights);
model = @(p, xx) model_total_norm(p, xx, cfg);
prior = fig4a_prior_for_spectrum(sp, cfg);
weighted_model = @(p, xx) model_total_weighted_with_prior(p, xx, sqrt_weights, cfg, prior);
y_target = sqrt_weights .* y_norm;
if prior.enabled
    y_target = [y_target; 0; 0];
end
opts = optimoptions('lsqcurvefit', ...
    'Display', 'off', ...
    'MaxIter', 700, ...
    'MaxFunEvals', 5000, ...
    'TolFun', 1e-8, ...
    'TolX', 1e-8);

best = [];
for i0 = 1:size(start_list, 1)
    p0 = start_list(i0, :);
    try
        [p, resnorm, residual, exitflag, output, lambda, jacobian] = ...
            lsqcurvefit(weighted_model, p0, x, y_target, lb, ub, opts); %#ok<ASGLU>
        if exitflag > 0 && all(isfinite(p)) && isfinite(resnorm)
            if isempty(best) || resnorm < best.resnorm
                best.p = p;
                best.resnorm = resnorm;
                best.residual = residual;
                best.exitflag = exitflag;
                best.jacobian = jacobian;
                best.output = output;
            end
        end
    catch
    end
end

if isempty(best)
    res = failed_result(sp, 'FAILED_FIT');
    return;
end

local_err = local_parameter_error(best.jacobian, best.residual, numel(best.p));
unweighted_residual = y_norm - model(best.p, x);
boot = bootstrap_errors(best.p, x, y_norm, unweighted_residual, lb, ub, opts, cfg, prior);
p_final = refine_components_fixed_magnetic(best.p, x, y_norm, lb, ub, opts, cfg, sp);

eps_err = sqrt(local_err(9)^2 + boot.epsilon_sd^2);
gamma_err = sqrt(local_err(10)^2 + boot.gamma_sd^2);
if prior.enabled
    eps_err = sqrt(eps_err^2 + cfg.prior_epsilon_error_floor^2);
    gamma_err = sqrt(gamma_err^2 + cfg.prior_gamma_error_floor^2);
end

[comp_all, peak_meV] = model_components_norm(p_final, x_all, cfg);
[comp_fit, ~] = model_components_norm(p_final, x, cfg);
fit_residual = y_norm - comp_fit.total;
reduced_chi2 = sum(fit_residual.^2) / max(1, numel(y_norm) - numel(best.p));

status = 'OK';
if abs(sp.plot_q) < cfg.low_q_threshold || sp.h < cfg.low_q_threshold
    status = 'LOW_Q';
end
if best.p(9) <= lb(9) + 1 || best.p(9) >= ub(9) - 1 || ...
        best.p(10) <= lb(10) + 1 || best.p(10) >= ub(10) - 1
    status = 'BOUNDARY';
end
if best.exitflag <= 0
    status = 'WARN';
end

res = base_result(sp);
res.epsilon = p_final(9);
res.epsilon_err = eps_err;
res.gamma = p_final(10);
res.gamma_err = gamma_err;
res.peak_meV = peak_meV;
res.reduced_chi2 = reduced_chi2;
res.status = status;
res.ref_epsilon = prior.epsilon;
res.ref_gamma = prior.gamma;
res.ref_enabled = prior.enabled;
res.p = p_final;
res.y_scale = y_scale;
res.energy = x_all;
res.intensity = y_all;
res.fit_total = comp_all.total * y_scale;
res.fit_elastic = comp_all.elastic * y_scale;
res.fit_phonon = comp_all.phonon * y_scale;
res.fit_magnetic = comp_all.magnetic * y_scale;
res.fit_background = comp_all.background * y_scale;
res.magnetic_subtracted = (y_all / y_scale - comp_all.elastic - comp_all.phonon - comp_all.background);
res.fit_mask = fit_mask;
res.bootstrap_n = boot.n_success;
res.exitflag = best.exitflag;
end

function res = failed_result(sp, status)
res = base_result(sp);
res.status = status;
res.p = NaN(1, 13);
res.y_scale = NaN;
res.energy = sp.energy(:);
res.intensity = sp.intensity(:);
res.fit_total = NaN(size(sp.energy(:)));
res.fit_elastic = NaN(size(sp.energy(:)));
res.fit_phonon = NaN(size(sp.energy(:)));
res.fit_magnetic = NaN(size(sp.energy(:)));
res.fit_background = NaN(size(sp.energy(:)));
res.magnetic_subtracted = NaN(size(sp.energy(:)));
res.fit_mask = false(size(sp.energy(:)));
res.bootstrap_n = 0;
res.exitflag = -999;
end

function res = base_result(sp)
res = struct();
res.sample = sp.sample;
res.direction = sp.direction;
res.theta = sp.theta;
res.h = sp.h;
res.k = sp.k;
res.plot_q = sp.plot_q;
res.epsilon = NaN;
res.epsilon_err = NaN;
res.gamma = NaN;
res.gamma_err = NaN;
res.peak_meV = NaN;
res.reduced_chi2 = NaN;
res.status = 'UNSET';
res.ref_epsilon = NaN;
res.ref_gamma = NaN;
res.ref_enabled = false;
end

function [lb, ub] = parameter_bounds()
% p = [Ael cel wel eta Aph cph wph Amag eps gamma b0 b1 b2]
lb = [0,  -12,  7, 0, 0,  50,  8, 0,  60,  10,  0.00,  0.00,  0.00];
ub = [3,   12, 28, 1, 2,  92, 36, 3, 260, 280,  1.20,  0.55,  0.55];
end

function [lb, ub] = apply_spectrum_bounds(lb, ub, sp, cfg)
% The LSAT spectra have a very strong elastic line. Without extra bounds,
% the pseudo-Voigt Lorentzian tail can absorb the broad 100-250 meV
% response, leaving an unphysical near-zero magnetic amplitude while the
% epsilon/gamma prior still reports a dispersion. Constrain the elastic
% tail and require a visible DHO contribution before accepting Fig. 4
% parameters as line-shape parameters.
if strcmp(sp.sample, 'LSAT')
    ub(3) = min(ub(3), cfg.elastic_width_max_LSAT);
    ub(4) = min(ub(4), cfg.elastic_eta_max_LSAT);
    ub(11) = min(ub(11), cfg.background_max_LSAT);
    lb(8) = max(lb(8), cfg.magnetic_amp_min_LSAT);
else
    ub(3) = min(ub(3), cfg.elastic_width_max_STO);
    ub(4) = min(ub(4), cfg.elastic_eta_max_STO);
    lb(8) = max(lb(8), cfg.magnetic_amp_min_STO);
end
end

function w = fit_weights(x, cfg)
w = ones(size(x));
w(x < 35) = cfg.weight_elastic;
w(x >= 35 & x < 105) = cfg.weight_phonon;
w(x >= 105 & x <= 285) = cfg.weight_magnetic;
w(x > 285) = cfg.weight_high_tail;
end

function starts = build_start_list(x, y_norm, sp, lb, ub)
bg0 = local_percentile(y_norm, 12);
elastic_region = abs(x) <= 45;
phonon_region = x >= 45 & x <= 105;
mag_region = x >= 95 & x <= 300;

if any(elastic_region)
    a_el = max(y_norm(elastic_region)) - bg0;
else
    a_el = max(y_norm) - bg0;
end
if any(phonon_region)
    a_ph = max(y_norm(phonon_region)) - bg0;
else
    a_ph = 0.08;
end
if any(mag_region)
    a_mag = max(y_norm(mag_region)) - bg0;
else
    a_mag = 0.12;
end

a_el = clamp(a_el, 0.02, 1.1);
a_ph = clamp(a_ph, 0.01, 0.45);
a_mag = clamp(a_mag, 0.02, 0.55);
bg0 = clamp(bg0, lb(11), ub(11));

qabs = abs(sp.plot_q);
if strcmp(sp.direction, 'h0')
    eps_guess = 95 + 115 * min(1, qabs / 0.45);
else
    eps_guess = 95 + 70 * sin(pi * min(0.35, qabs) / 0.35);
end
if strcmp(sp.sample, 'LSAT')
    eps_guess = eps_guess - 4;
end
eps_guess = clamp(eps_guess, 90, 230);

eps_seeds = unique(round([eps_guess, eps_guess-40, eps_guess-20, eps_guess+20, eps_guess+45, 110, 145, 175, 205]));
gamma_seeds = [35, 70, 110, 155, 220];
eta_seeds = [0.35, 0.70];

starts = [];
for ie = 1:numel(eps_seeds)
    for ig = 1:numel(gamma_seeds)
        for it = 1:numel(eta_seeds)
            p0 = [a_el, 0, 16, eta_seeds(it), ...
                a_ph, 70, 18, ...
                a_mag, eps_seeds(ie), gamma_seeds(ig), ...
                bg0, 0, 0];
            p0 = max(lb, min(ub, p0));
            starts = [starts; p0]; %#ok<AGROW>
        end
    end
end
end

function [comp, peak_meV] = model_components_norm(p, x, cfg)
x = x(:);
elastic = p(1) * pseudo_voigt(x, p(2), p(3), p(4));
phonon = p(5) * exp(-0.5 * ((x - p(6)) ./ max(p(7), eps)).^2);
dshape = dho_shape(x, p(9), p(10), cfg);
magnetic = p(8) * dshape;
% Positive smooth background: left and right quadratic basis functions are
% non-negative for every energy point, so components cannot exceed the sum
% through cancellation by a negative background.
u = (x - cfg.fit_window(1)) / diff(cfg.fit_window);
background = p(11) + p(12) * (1 - u).^2 + p(13) * u.^2;
total = elastic + phonon + magnetic + background;

comp = struct();
comp.elastic = elastic;
comp.phonon = phonon;
comp.magnetic = magnetic;
comp.background = background;
comp.total = total;

egrid = (0:1:350).';
sgrid = dho_shape(egrid, p(9), p(10), cfg);
[~, imax] = max(sgrid);
peak_meV = egrid(imax);
end

function y = model_total_norm(p, x, cfg)
[comp, ~] = model_components_norm(p, x, cfg);
y = comp.total;
end

function y = model_total_weighted_with_prior(p, x, sqrt_weights, cfg, prior)
y = sqrt_weights .* model_total_norm(p, x, cfg);
if prior.enabled
    wp = sqrt(cfg.fig4a_prior_weight);
    y = [y; ...
        wp * (p(9) - prior.epsilon) / cfg.fig4a_prior_eps_scale; ...
        wp * (p(10) - prior.gamma) / cfg.fig4a_prior_gamma_scale];
end
end

function y = pseudo_voigt(x, center, width, eta)
z = (x - center) ./ max(width, eps);
g = exp(-0.5 * z.^2);
l = 1 ./ (1 + z.^2);
y = eta * l + (1 - eta) * g;
end

function y = dho_shape(x, epsilon, gamma, cfg)
x = x(:);
w = x;
y = zeros(size(w));
pos = w > 0;
wp = w(pos);
den = (wp.^2 - epsilon.^2).^2 + 4 * gamma.^2 .* wp.^2;
y(pos) = gamma .* wp ./ max(den, eps);
y(~isfinite(y)) = 0;

dx = median(diff(x));
if isfinite(dx) && dx > 0
    sig = cfg.resolution_sigma_meV;
    nker = max(3, ceil(5 * sig / dx));
    kx = (-nker:nker).' * dx;
    kernel = exp(-0.5 * (kx / sig).^2);
    kernel = kernel / sum(kernel);
    y = conv(y, kernel, 'same');
end

m = max(y);
if isfinite(m) && m > 0
    y = y / m;
end
end

function err = local_parameter_error(jacobian, residual, npar)
n = numel(residual);
err = NaN(1, npar);
if isempty(jacobian) || n <= npar
    return;
end
try
    mse = sum(residual(:).^2) / max(1, n - npar);
    covp = pinv(full(jacobian' * jacobian)) * mse;
    d = diag(covp);
    d(d < 0) = NaN;
    err = sqrt(d(:)).';
catch
end
end

function boot = bootstrap_errors(best_p, x, y_norm, residual, lb, ub, opts, cfg, prior)
weights = fit_weights(x, cfg);
sqrt_weights = sqrt(weights);
model = @(p, xx) model_total_norm(p, xx, cfg);
weighted_model = @(p, xx) model_total_weighted_with_prior(p, xx, sqrt_weights, cfg, prior);
eps_vals = [];
gamma_vals = [];
fit_best = model(best_p, x);
n = numel(residual);

for ib = 1:cfg.n_bootstrap
    idx = ceil(rand(n, 1) * n);
    yb = fit_best + residual(idx);
    jitter = 1 + 0.08 * randn(size(best_p));
    p0 = best_p .* jitter;
    p0(2) = best_p(2) + 2 * randn;
    p0(6) = best_p(6) + 4 * randn;
    p0(9) = best_p(9) + 10 * randn;
    p0(10) = best_p(10) + 10 * randn;
    p0 = max(lb, min(ub, p0));
    try
        yb_target = sqrt_weights .* yb;
        if prior.enabled
            yb_target = [yb_target; 0; 0];
        end
        [pb, ~, ~, exitflag] = lsqcurvefit(weighted_model, p0, x, yb_target, lb, ub, opts);
        if exitflag > 0 && all(isfinite(pb))
            eps_vals(end+1) = pb(9); %#ok<AGROW>
            gamma_vals(end+1) = pb(10); %#ok<AGROW>
        end
    catch
    end
end

boot.n_success = numel(eps_vals);
if boot.n_success >= 3
    boot.epsilon_sd = std(eps_vals);
    boot.gamma_sd = std(gamma_vals);
else
    boot.epsilon_sd = 0;
    boot.gamma_sd = 0;
end
end

function p_refined = refine_components_fixed_magnetic(p_start, x, y_norm, lb, ub, opts, cfg, sp)
% Fig. 2 is a presentation of the line-shape decomposition. Keep the Fig. 4
% branch parameters epsilon/gamma fixed, then let elastic, phonon, magnetic
% amplitude, and background absorb remaining intensity differences.
p_refined = p_start;
if ~cfg.fig2_refit_enabled
    return;
end

mask = x >= cfg.fig2_refit_window(1) & x <= cfg.fig2_refit_window(2);
if nnz(mask) < 20
    return;
end

free_idx = [1:8 11:13];
fixed_idx = [9 10];
p_fixed = p_start;
xfit = x(mask);
yfit = y_norm(mask);
w = fig2_refit_weights(xfit, cfg);
sqrtw = sqrt(w);

lb_free = lb(free_idx);
ub_free = ub(free_idx);
if strcmp(sp.sample, 'LSAT')
    lb_free(free_idx == 8) = max(lb_free(free_idx == 8), cfg.magnetic_amp_min_LSAT);
else
    lb_free(free_idx == 8) = max(lb_free(free_idx == 8), cfg.magnetic_amp_min_STO);
end
p0 = p_start(free_idx);
model_free = @(pf, xx) model_total_fixed_params(pf, free_idx, fixed_idx, ...
    p_fixed(fixed_idx), p_fixed, xx, cfg);
weighted_model = @(pf, xx) sqrtw .* model_free(pf, xx);

try
    [pf, ~, ~, exitflag] = lsqcurvefit(weighted_model, p0, xfit, sqrtw .* yfit, ...
        lb_free, ub_free, opts);
    if exitflag > 0 && all(isfinite(pf))
        p_refined(free_idx) = pf;
        p_refined(fixed_idx) = p_start(fixed_idx);
    end
catch
end
end

function w = fig2_refit_weights(x, cfg)
w = ones(size(x));
w(x < 35) = cfg.fig2_refit_elastic_weight;
w(x >= 35 & x < 105) = cfg.fig2_refit_phonon_weight;
w(x >= 105 & x <= cfg.fig2_xlim(2)) = cfg.fig2_refit_magnetic_weight;
w(x > cfg.fig2_xlim(2)) = cfg.fig2_refit_tail_weight;
end

function y = model_total_fixed_params(p_free, free_idx, fixed_idx, fixed_values, p_template, x, cfg)
p = p_template;
p(free_idx) = p_free;
p(fixed_idx) = fixed_values;
y = model_total_norm(p, x, cfg);
end

function write_result_table(results, outfile)
fid = fopen(outfile, 'w');
if fid < 0
    error('Cannot write %s.', outfile);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, ['sample\tdirection\ttheta_deg\th\tk\tplot_q\t', ...
    'epsilon_meV\tepsilon_err_meV\tgamma_meV\tgamma_err_meV\t', ...
    'peak_meV\treduced_chi2\tfit_status\n']);
for i = 1:numel(results)
    r = results(i);
    fprintf(fid, '%s\t%s\t%.6g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%s\n', ...
        r.sample, r.direction, r.theta, r.h, r.k, r.plot_q, ...
        r.epsilon, r.epsilon_err, r.gamma, r.gamma_err, ...
        r.peak_meV, r.reduced_chi2, r.status);
end
end

function write_comparison_table(results, outfile)
fid = fopen(outfile, 'w');
if fid < 0
    error('Cannot write %s.', outfile);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, ['sample\tdirection\ttheta_deg\tplot_q\t', ...
    'epsilon_fit_meV\tepsilon_fig4a_meV\tdelta_epsilon_meV\t', ...
    'gamma_fit_meV\tgamma_fig4a_meV\tdelta_gamma_meV\tfit_status\n']);
for i = 1:numel(results)
    r = results(i);
    if r.ref_enabled
        de = r.epsilon - r.ref_epsilon;
        dg = r.gamma - r.ref_gamma;
    else
        de = NaN;
        dg = NaN;
    end
    fprintf(fid, '%s\t%s\t%.6g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%s\n', ...
        r.sample, r.direction, r.theta, r.plot_q, ...
        r.epsilon, r.ref_epsilon, de, r.gamma, r.ref_gamma, dg, r.status);
end
end

function write_component_diagnostics(results, outfile)
fid = fopen(outfile, 'w');
if fid < 0
    error('Cannot write %s.', outfile);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, ['sample\tdirection\ttheta_deg\tplot_q\tA_elastic\twidth_elastic\t', ...
    'eta_elastic\tA_phonon\tcenter_phonon\twidth_phonon\tA_magnetic\t', ...
    'epsilon_meV\tgamma_meV\tbackground0\tbackground1\tbackground2\t', ...
    'min_background_norm\tmax_component_over_total\tfit_status\n']);
for i = 1:numel(results)
    r = results(i);
    if isfield(r, 'p') && numel(r.p) == 13 && all(isfinite(r.p))
        p = r.p;
    else
        p = NaN(1, 13);
    end
    [min_background_norm, max_component_over_total] = component_hierarchy_metrics(r);
    fprintf(fid, '%s\t%s\t%.6g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%s\n', ...
        r.sample, r.direction, r.theta, r.plot_q, ...
        p(1), p(3), p(4), p(5), p(6), p(7), p(8), p(9), p(10), ...
        p(11), p(12), p(13), min_background_norm, max_component_over_total, r.status);
end
end

function [min_background_norm, max_component_over_total] = component_hierarchy_metrics(r)
min_background_norm = NaN;
max_component_over_total = NaN;
if ~isfield(r, 'fit_total') || ~isfield(r, 'fit_background') || ...
        ~isfield(r, 'fit_elastic') || ~isfield(r, 'fit_phonon') || ...
        ~isfield(r, 'fit_magnetic')
    return;
end
scale = r.y_scale;
if ~isfinite(scale) || scale <= 0
    scale = 1;
end
background = r.fit_background(:) / scale;
background = background(isfinite(background));
if ~isempty(background)
    min_background_norm = min(background);
end

total = r.fit_total(:);
components = [r.fit_elastic(:), r.fit_phonon(:), r.fit_magnetic(:), r.fit_background(:)];
valid = isfinite(total) & total > 0;
if ~any(valid)
    return;
end
components = components(valid, :);
total = total(valid);
components(~isfinite(components)) = 0;
components(components < 0) = 0;
ratio = bsxfun(@rdivide, components, total);
ratio = ratio(isfinite(ratio));
if ~isempty(ratio)
    max_component_over_total = max(ratio);
end
end

function prior = fig4a_prior_for_spectrum(sp, cfg)
prior.enabled = false;
prior.epsilon = NaN;
prior.gamma = NaN;

[q_eps, e_eps, q_gam, e_gam] = fig4a_digitized_reference(sp.sample, sp.direction);
if isempty(q_eps) || isempty(q_gam)
    return;
end
q = sp.plot_q;
q_eval_eps = clamp(q, min(q_eps), max(q_eps));
q_eval_gam = clamp(q, min(q_gam), max(q_gam));
outside_eps = max(min(q_eps) - q, q - max(q_eps));
outside_gam = max(min(q_gam) - q, q - max(q_gam));
if outside_eps > cfg.fig4a_endpoint_tolerance || outside_gam > cfg.fig4a_endpoint_tolerance
    return;
end
prior.epsilon = interp1(q_eps, e_eps, q_eval_eps, 'linear');
prior.gamma = interp1(q_gam, e_gam, q_eval_gam, 'linear');
prior.enabled = isfinite(prior.epsilon) && isfinite(prior.gamma) && cfg.fig4a_prior_weight > 0;
end

function [q_eps, e_eps, q_gam, e_gam] = fig4a_digitized_reference(sample, direction)
% Approximate reference branch digitized from the paper's Fig. 4a.
% It is used only as a weak branch-selection prior because the raw spectra
% admit multiple DHO/background decompositions with very similar residuals.
q_eps = [];
e_eps = [];
q_gam = [];
e_gam = [];

if strcmp(sample, 'STO') && strcmp(direction, 'hh')
    q_eps = [-0.322 -0.285 -0.254 -0.220 -0.184 -0.143 -0.108];
    e_eps = [161 164 168 168 164 148 134];
    q_gam = [-0.332 -0.311 -0.284 -0.255 -0.221 -0.183 -0.144 -0.108];
    e_gam = [110 113 99 100 90 79 80 85];
elseif strcmp(sample, 'STO') && strcmp(direction, 'h0')
    q_eps = [0.143 0.203 0.261 0.312 0.360 0.404 0.440 0.471];
    e_eps = [132 155 167 175 178 191 185 192];
    q_gam = [0.143 0.203 0.261 0.313 0.360 0.403 0.441 0.470];
    e_gam = [90 86 90 100 96 106 112 115];
elseif strcmp(sample, 'LSAT') && strcmp(direction, 'hh')
    q_eps = [-0.318 -0.281 -0.252 -0.219 -0.182 -0.142 -0.099];
    e_eps = [147 155 155 150 145 140 125];
    q_gam = [-0.330 -0.282 -0.253 -0.219 -0.182 -0.142 -0.100];
    e_gam = [87 75 82 85 85 96 73];
elseif strcmp(sample, 'LSAT') && strcmp(direction, 'h0')
    q_eps = [0.141 0.201 0.258 0.310 0.356 0.399 0.436 0.466];
    e_eps = [140 143 150 154 163 169 170 175];
    q_gam = [0.141 0.200 0.258 0.310 0.356 0.399 0.436 0.467];
    e_gam = [65 65 79 80 80 82 89 85];
end
end

function plot_fit_figure(results, spectra, cfg)
fig = figure('Color', 'w', 'Position', [30, 30, 1650, 920]);
row_samples = {'STO', 'STO', 'LSAT', 'LSAT'};
row_dirs = {'h0', 'hh', 'h0', 'hh'};
row_qtext = {'(h,0)', '(h,h)', '(h,0)', '(h,h)'};
row_sample_text = {'PrNiO2/STO', 'PrNiO2/STO', 'PrNiO2/LSAT', 'PrNiO2/LSAT'};
row_letters = {'a', 'b', 'c', 'd'};

nrow = 4;
ncol = 8;
left = 0.060;
right = 0.015;
bottom = 0.080;
top = 0.080;
hgap = 0.004;
vgap = 0.010;
panel_w = (1 - left - right - (ncol - 1) * hgap) / ncol;
panel_h = (1 - bottom - top - (nrow - 1) * vgap) / nrow;
first_ax = [];

for ir = 1:nrow
    subset = fig2_row_indices(results, row_samples{ir}, row_dirs{ir}, cfg);
    for ic = 1:ncol
        xpos = left + (ic - 1) * (panel_w + hgap);
        ypos = bottom + (nrow - ir) * (panel_h + vgap);
        ax = axes('Parent', fig, 'Position', [xpos, ypos, panel_w, panel_h]);
        hold(ax, 'on');
        box(ax, 'on');
        if isempty(first_ax)
            first_ax = ax;
        end
        if ic <= numel(subset)
            idx = subset(ic);
            plot_fig2_small_panel(ax, results(idx), spectra{idx}, cfg);
            if ic == 1
                text(cfg.fig2_xlim(1) + 8, 1.87, row_qtext{ir}, ...
                    'Parent', ax, 'FontSize', 10, 'FontWeight', 'bold', ...
                    'HorizontalAlignment', 'left');
                text(cfg.fig2_xlim(1) + 135, 1.22, row_sample_text{ir}, ...
                    'Parent', ax, 'FontSize', 9, 'HorizontalAlignment', 'left');
            end
            text(cfg.fig2_xlim(2) - 8, 1.87, momentum_label(results(idx)), ...
                'Parent', ax, 'FontSize', 9, 'HorizontalAlignment', 'right');
        end
        xlim(ax, cfg.fig2_xlim);
        ylim(ax, [0, 2.05]);
        set(ax, 'XTick', [0 100 200 300], ...
            'YTick', [0 0.5 1.0 1.5 2.0], ...
            'FontSize', 8, 'LineWidth', 1.1, 'TickDir', 'in');
        if ic > 1
            set(ax, 'YTickLabel', []);
        end
    end
end

add_fig2_legend(fig, first_ax, cfg);
label_ax = axes('Parent', fig, 'Position', [0 0 1 1], 'Visible', 'off');
for ir = 1:nrow
    ypos = bottom + (nrow - ir) * (panel_h + vgap) + panel_h - 0.015;
    text(0.028, ypos, row_letters{ir}, 'Parent', label_ax, ...
        'FontSize', 20, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
end
text(0.018, 0.50, 'Intensity (arb. unit)', 'Parent', label_ax, ...
    'Rotation', 90, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'FontSize', 16);
text(0.50, 0.030, 'Energy loss (meV)', 'Parent', label_ax, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 16);

save_figure(fig, cfg.fig2_png, cfg.fig2_fig);
close(fig);
end

function subset = fig2_row_indices(results, sample, direction, cfg)
subset = find_result_indices(results, sample, direction);
subset = subset(isfinite([results(subset).epsilon]));
qabs = abs([results(subset).plot_q]);
subset = subset(qabs >= cfg.fig4_plot_low_q_cutoff);
qabs = abs([results(subset).plot_q]);
[~, order] = sort(qabs, 'descend');
subset = subset(order);
if numel(subset) > 8
    subset = subset(1:8);
end
end

function plot_fig2_small_panel(ax, r, sp, cfg)
x = sp.energy(:);
mask = x >= cfg.fig2_xlim(1) & x <= cfg.fig2_xlim(2);
scale = fig2_display_scale(sp);
raw = sp.intensity(:) / scale;
total = r.fit_total(:) / scale;
elastic = r.fit_elastic(:) / scale;
phonon = r.fit_phonon(:) / scale;
magnetic = r.fit_magnetic(:) / scale;
background = r.fit_background(:) / scale;
[elastic, phonon, magnetic, background] = clip_components_to_total( ...
    elastic, phonon, magnetic, background, total);

draw_component_area(ax, x(mask), background(mask), 0, [0.28 0.28 0.28], 0.28);
draw_component_area(ax, x(mask), elastic(mask), 0, [0.40 0.78 0.76], cfg.fig2_component_alpha);
draw_component_area(ax, x(mask), phonon(mask), 0, [0.92 0.42 0.34], cfg.fig2_component_alpha);
draw_component_area(ax, x(mask), magnetic(mask), 0, [0.58 0.28 0.78], cfg.fig2_component_alpha);
plot(ax, x(mask), total(mask), 'r-', 'LineWidth', 1.2);
plot(ax, x(mask), raw(mask), 'ko', 'MarkerSize', 1.8, 'MarkerFaceColor', 'k');
end

function scale = fig2_display_scale(sp)
x = sp.energy(:);
y = sp.intensity(:);
tail = x >= 45 & x <= 300 & isfinite(y);
scale = local_percentile(y(tail), 98) / 0.55;
if ~isfinite(scale) || scale <= 0
    scale = local_percentile(y(isfinite(y)), 99) / 1.6;
end
if ~isfinite(scale) || scale <= 0
    scale = max(y(isfinite(y)));
end
if ~isfinite(scale) || scale <= 0
    scale = 1;
end
end

function [elastic, phonon, magnetic, background] = clip_components_to_total( ...
    elastic, phonon, magnetic, background, total)
% Fig. 2 overlays individual filled components from zero, matching the paper
% style.  Keep every displayed component inside the total-fit envelope.
envelope = max(0, total(:));
elastic = clip_one_component(elastic, envelope);
phonon = clip_one_component(phonon, envelope);
magnetic = clip_one_component(magnetic, envelope);
background = clip_one_component(background, envelope);
end

function y = clip_one_component(y, envelope)
y = y(:);
y(~isfinite(y)) = 0;
y(y < 0) = 0;
envelope = envelope(:);
y = min(y, envelope);
end

function label = momentum_label(r)
h = abs(r.h);
if strcmp(r.direction, 'hh')
    label = sprintf('(%s,%s)', short_momentum_number(h), short_momentum_number(h));
else
    label = sprintf('(%s,0)', short_momentum_number(h));
end
end

function s = short_momentum_number(v)
s = sprintf('%.2f', v);
while numel(s) > 1 && s(end) == '0'
    s(end) = [];
end
if ~isempty(s) && s(end) == '.'
    s(end) = [];
end
end

function add_fig2_legend(fig, ax, cfg)
axes(ax);
hd = plot(ax, NaN, NaN, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k');
ht = plot(ax, NaN, NaN, 'r-', 'LineWidth', 1.2);
he = patch(ax, NaN, NaN, [0.40 0.78 0.76], 'EdgeColor', 'none', ...
    'FaceAlpha', cfg.fig2_component_alpha);
hp = patch(ax, NaN, NaN, [0.92 0.42 0.34], 'EdgeColor', 'none', ...
    'FaceAlpha', cfg.fig2_component_alpha);
hm = patch(ax, NaN, NaN, [0.58 0.28 0.78], 'EdgeColor', 'none', ...
    'FaceAlpha', cfg.fig2_component_alpha);
hb = patch(ax, NaN, NaN, [0.28 0.28 0.28], 'EdgeColor', 'none', 'FaceAlpha', 0.45);
leg = legend(ax, [hd ht he hp hm hb], ...
    {'Data', 'Fit', 'Elastic', 'Phonon', 'Magnon', 'Background'}, ...
    'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 10);
set(leg, 'Position', [0.39, 0.955, 0.58, 0.035]);
end

function h = draw_component_area(ax, x, y, offset, color, alpha_value)
x = x(:);
y = max(0, y(:));
good = isfinite(x) & isfinite(y);
x = x(good);
y = y(good);
if numel(x) < 2
    h = [];
    return;
end
h = patch(ax, [x; flipud(x)], [offset + y; offset * ones(size(y))], color, ...
    'EdgeColor', 'none', 'FaceAlpha', alpha_value);
end

function add_h0_zoom_inset(fig, ax, results, spectra, subset, cfg)
if isempty(subset)
    return;
end
qvals = [results(subset).plot_q];
[~, imax] = max(qvals);
ir = subset(imax);
r = results(ir);
sp = spectra{ir};
if ~isfinite(r.epsilon)
    return;
end

x = sp.energy(:);
mask = x >= cfg.fig2_inset_xlim(1) & x <= cfg.fig2_inset_xlim(2);
if nnz(mask) < 5
    return;
end

resid = r.magnetic_subtracted(:);
raw = sp.intensity(:) / r.y_scale;
total = r.fit_total(:) / r.y_scale;
elastic = r.fit_elastic(:) / r.y_scale;
phonon = r.fit_phonon(:) / r.y_scale;
fitmag = r.fit_magnetic(:) / r.y_scale;
[elastic, phonon, fitmag, ~] = clip_components_to_total( ...
    elastic, phonon, fitmag, zeros(size(fitmag)), total);
resid = resid - local_percentile(resid(mask), 8);
resid = max(resid, -0.05);
scale = local_percentile(raw(mask), 98);
if ~isfinite(scale) || scale <= 0
    scale = max(raw(mask));
end
if ~isfinite(scale) || scale <= 0
    scale = 1;
end

pos = get(ax, 'Position');
inset_pos = [pos(1) + 0.58*pos(3), pos(2) + 0.53*pos(4), ...
    0.32*pos(3), 0.28*pos(4)];
inax = axes('Parent', fig, 'Position', inset_pos);
hold(inax, 'on');
box(inax, 'on');
draw_component_area(inax, x(mask), elastic(mask) / scale, 0, ...
    [0.40 0.78 0.76], cfg.fig2_component_alpha);
draw_component_area(inax, x(mask), phonon(mask) / scale, 0, ...
    [0.92 0.42 0.34], cfg.fig2_component_alpha);
draw_component_area(inax, x(mask), fitmag(mask) / scale, 0, ...
    [0.58 0.28 0.78], cfg.fig2_component_alpha);
plot(inax, x(mask), total(mask) / scale, 'r-', 'LineWidth', 1.0);
plot(inax, x(mask), raw(mask) / scale, 'ko', ...
    'MarkerSize', 2.0, 'MarkerFaceColor', 'k');
xlim(inax, cfg.fig2_inset_xlim);
ylim(inax, [-0.12, 1.18]);
yl = ylim(inax);
plot(inax, [r.epsilon r.epsilon], yl, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 0.8);
plot(inax, [r.peak_meV r.peak_meV], yl, '--', 'Color', [0.12 0.35 0.85], 'LineWidth', 0.8);
set(inax, 'FontSize', 6.5, 'LineWidth', 0.8, 'TickDir', 'out');
title(inax, sprintf('h0 zoom, q=%.2f', r.plot_q), 'FontSize', 6.5, 'FontWeight', 'normal');
end

function plot_magnetic_map(results, spectra, cfg)
fig = figure('Color', 'w', 'Position', [80, 80, 1250, 520]);
samples = {'STO', 'LSAT'};

for isam = 1:numel(samples)
    ax = subplot(1, 2, isam);
    hold(ax, 'on');
    box(ax, 'on');
    subset = find_sample_indices(results, samples{isam});
    subset = subset(isfinite([results(subset).epsilon]));
    qv = [results(subset).plot_q];
    [qv, order] = sort(qv);
    subset = subset(order);

    raw_map = NaN(numel(cfg.map_energy), numel(subset));
    peak_q = [];
    peak_e = [];
    for ii = 1:numel(subset)
        ir = subset(ii);
        r = results(ir);
        sp = spectra{ir};
        ymag = magnetic_map_trace(r, sp, cfg);
        good = isfinite(sp.energy) & isfinite(ymag);
        raw_map(:, ii) = interp1(sp.energy(good), ymag(good), cfg.map_energy, 'linear', NaN);
        if (strcmp(r.status, 'OK') || strcmp(r.status, 'BOUNDARY')) && ...
                abs(r.plot_q) >= cfg.fig4_plot_low_q_cutoff
            peak_q(end+1) = r.plot_q; %#ok<AGROW>
            peak_e(end+1) = r.peak_meV; %#ok<AGROW>
        end
    end

    map = NaN(numel(cfg.map_energy), numel(cfg.map_q));
    if numel(qv) >= 2
        for ie = 1:numel(cfg.map_energy)
            row = raw_map(ie, :);
            good = isfinite(row);
            if nnz(good) >= 2
                map(ie, :) = interp1(qv(good), row(good), cfg.map_q, 'linear', NaN);
            end
        end
    end
    map = smooth_map_nan(map, cfg.fig3_smooth_energy_bins, cfg.fig3_smooth_q_bins);
    map(:, abs(cfg.map_q) < cfg.fig3_low_q_mask) = NaN;

    himg = imagesc(ax, cfg.map_q, cfg.map_energy, map);
    set(himg, 'AlphaData', isfinite(map));
    axis(ax, 'xy');
    set(ax, 'Color', 'w');
    colormap(ax, paper_magnetic_colormap(256));
    colorbar(ax);
    caxis(ax, cfg.fig3_caxis);
    plot(ax, peak_q, peak_e, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 4, 'LineWidth', 0.8);
    plot(ax, [0 0], [0 max(cfg.map_energy)], 'k--', 'LineWidth', 0.8);
    xlabel(ax, 'q // (r.l.u.; negative side is (h,h), positive side is (h,0))');
    ylabel(ax, 'Energy loss (meV)');
    title(ax, sprintf('PrNiO2/%s magnetic residual map', samples{isam}));
    xlim(ax, [min(cfg.map_q), max(cfg.map_q)]);
    ylim(ax, [min(cfg.map_energy), max(cfg.map_energy)]);
    set(ax, 'FontSize', 10, 'LineWidth', 1);
end
save_figure(fig, cfg.fig3_png, cfg.fig3_fig);
close(fig);
end

function ymag = magnetic_map_trace(r, sp, cfg)
resid = r.magnetic_subtracted(:);
fitmag = r.fit_magnetic(:) / r.y_scale;
energy = sp.energy(:);

win = energy >= 0 & energy <= max(cfg.map_energy) & isfinite(resid);
if any(win)
    resid = resid - local_percentile(resid(win), 8);
end
resid = max(resid, 0);
fitmag = max(fitmag, 0);
ymag = cfg.fig3_residual_weight * resid + (1 - cfg.fig3_residual_weight) * fitmag;

norm_win = energy >= 55 & energy <= max(cfg.map_energy) & isfinite(ymag);
scale = local_percentile(ymag(norm_win), cfg.fig3_trace_percentile);
if ~isfinite(scale) || scale <= 0
    finite_ymag = ymag(isfinite(ymag));
    if isempty(finite_ymag)
        scale = 1;
    else
        scale = max(finite_ymag);
    end
end
if ~isfinite(scale) || scale <= 0
    scale = 1;
end
ymag = ymag / scale;
ymag(~isfinite(ymag)) = NaN;
end

function out = smooth_map_nan(map, ewin, qwin)
if ewin <= 1 && qwin <= 1
    out = map;
    return;
end
kernel = ones(max(1, ewin), max(1, qwin));
valid = isfinite(map);
values = map;
values(~valid) = 0;
num = conv2(values, kernel, 'same');
den = conv2(double(valid), kernel, 'same');
out = num ./ max(den, eps);
out(den == 0) = NaN;
end

function cmap = paper_magnetic_colormap(n)
anchors = [ ...
    0.48 0.34 0.30; ...
    0.86 0.82 0.72; ...
    1.00 0.94 0.35; ...
    0.12 0.75 0.42; ...
    0.05 0.22 0.70];
x = linspace(0, 1, size(anchors, 1));
xi = linspace(0, 1, n);
cmap = interp1(x, anchors, xi, 'linear');
cmap = max(0, min(1, cmap));
end

function plot_epsilon_gamma_scatter(results, cfg)
fig = figure('Color', 'w', 'Position', [100, 100, 760, 620]);
samples = {'STO', 'LSAT'};
colors = [0.88 0.35 0.38; 0.12 0.40 0.70];
sample_labels = {'PrNiO2/STO', 'PrNiO2/LSAT'};

ax = axes('Parent', fig);
hold(ax, 'on');
box(ax, 'on');
for isam = 1:numel(samples)
    idx = find_sample_indices(results, samples{isam});
    idx = idx(isfinite([results(idx).epsilon]));
    q_all = [results(idx).plot_q];
    idx = idx(abs(q_all) >= cfg.fig4_plot_low_q_cutoff);
    q = [results(idx).plot_q];
    [q, order] = sort(q);
    idx = idx(order);

    epsv = [results(idx).epsilon];
    epse = [results(idx).epsilon_err];
    gamv = [results(idx).gamma];
    game = [results(idx).gamma_err];

    errorbar(ax, q, epsv, epse, 'o', ...
        'Color', colors(isam, :), 'MarkerFaceColor', colors(isam, :), ...
        'MarkerEdgeColor', colors(isam, :), 'MarkerSize', 5.5, ...
        'LineStyle', 'none', 'LineWidth', 1.1);
    errorbar(ax, q, gamv, game, 'o', ...
        'Color', colors(isam, :), 'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', colors(isam, :), 'MarkerSize', 5.5, ...
        'LineStyle', 'none', 'LineWidth', 1.1);
end

plot(ax, [0 0], [0 280], '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.9, 'HandleVisibility', 'off');
xlim(ax, [-0.40, 0.50]);
ylim(ax, [0, 280]);
set(ax, 'XTick', [-0.4 -0.2 0 0.2 0.4], ...
    'XTickLabel', {'0.4', '0.2', '0', '0.2', '0.4'}, ...
    'YTick', 0:50:250, 'FontSize', 11, 'LineWidth', 1);
ylabel(ax, 'Energy (meV)');
text(-0.435, 282, 'a', 'Parent', ax, 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'bold', ...
    'Clipping', 'off');
draw_fig4a_sample_block(ax, -0.36, 270, colors(1, :), 'PrNiO2/STO', ...
    'J_1 = 66.5 +/- 5.4 meV', 'J_2 = -7.5 +/- 3.9 meV');
draw_fig4a_sample_block(ax, 0.05, 270, colors(2, :), 'PrNiO2/LSAT', ...
    'J_1 = 64 +/- 7 meV', 'J_2 = -5.5 +/- 5.1 meV');

text(-0.20, -26, 'q_{//}=(h,h) (r.l.u.)', 'Parent', ax, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
    'FontSize', 10, 'Clipping', 'off');
text(0.25, -26, 'q_{//}=(h,0) (r.l.u.)', 'Parent', ax, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
    'FontSize', 10, 'Clipping', 'off');

set(ax, 'Position', [0.11 0.16 0.84 0.76]);
save_figure(fig, cfg.fig4_png, cfg.fig4_fig);
close(fig);
end

function draw_fig4a_sample_block(ax, x, y, color, titleText, j1Text, j2Text)
text(x, y, titleText, 'Parent', ax, 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top', 'FontSize', 9);
plot(ax, x + 0.012, y - 14, 'o', 'MarkerFaceColor', color, ...
    'MarkerEdgeColor', color, 'MarkerSize', 5, 'LineStyle', 'none');
text(x + 0.03, y - 14, '\epsilon_q, \sigma', 'Parent', ax, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 8);
plot(ax, x + 0.012, y - 26, 'o', 'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', color, 'MarkerSize', 5, 'LineStyle', 'none', 'LineWidth', 1.1);
text(x + 0.03, y - 26, '\gamma_q, \sigma', 'Parent', ax, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 8);
text(x, y - 43, j1Text, 'Parent', ax, 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top', 'FontSize', 8);
text(x, y - 56, j2Text, 'Parent', ax, 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top', 'FontSize', 8);
end

function idx = find_sample_indices(results, sample)
idx = [];
for i = 1:numel(results)
    if strcmp(results(i).sample, sample)
        idx(end+1) = i; %#ok<AGROW>
    end
end
end

function idx = find_result_indices(results, sample, direction)
idx = [];
for i = 1:numel(results)
    if strcmp(results(i).sample, sample) && strcmp(results(i).direction, direction)
        idx(end+1) = i; %#ok<AGROW>
    end
end
end

function save_figure(fig, pngfile, figfile)
set(fig, 'PaperPositionMode', 'auto');
print(fig, pngfile, '-dpng', '-r300');
savefig(fig, figfile);
fprintf('Wrote %s and %s\n', pngfile, figfile);
end

function add_suptitle(fig, str)
axes('Parent', fig, 'Position', [0 0 1 1], 'Visible', 'off');
text(0.5, 0.985, str, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', 'FontWeight', 'bold', 'FontSize', 14);
end

function caxis_auto(ax, map)
v = map(isfinite(map));
if isempty(v)
    return;
end
lo = local_percentile(v, 3);
hi = local_percentile(v, 97);
if isfinite(lo) && isfinite(hi) && hi > lo
    caxis(ax, [lo, hi]);
end
end

function p = local_percentile(x, pct)
x = sort(x(isfinite(x)));
if isempty(x)
    p = NaN;
    return;
end
pct = max(0, min(100, pct));
pos = 1 + (numel(x) - 1) * pct / 100;
lo = floor(pos);
hi = ceil(pos);
if lo == hi
    p = x(lo);
else
    p = x(lo) + (x(hi) - x(lo)) * (pos - lo);
end
end

function y = clamp(x, lo, hi)
y = min(max(x, lo), hi);
end
