% Run ONE Newtonian simulation of Gabbard et al. (2025) and dump every state
% variable, so the film pressure at contact onset can be read off a working
% simulation instead of argued about.
%
% Parameters are the production sweep's own: sigma, rho, Ro, nu, Bo taken verbatim
% from sweeper_experiments.m, and the largest of its five Weber numbers.

here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, 'LWDR', 'matlab', '1_code')));

sigma = 20.5; rho = 0.96; Ro = 0.0203; Bo = 0.0189;
nu = 0.01 * 20;                       % cgs, as in the sweep
g  = Bo * sigma / (rho * Ro^2);
We  = 1.0;
v_e = -sqrt(sigma / (rho * Ro) * We);
X   = 2 * g * Ro * 0.02;
v0  = -sqrt(v_e^2 + X);

Oh = nu * sqrt(rho / (sigma * Ro));
fprintf('We = %.4f   Oh = %.6f   v0 = %.4f cm/s   g = %.4g\n', We, Oh, v0, g);

physical = struct('undisturbed_radius', Ro, 'initial_velocity', v0, ...
                  'rhoS', rho, 'sigmaS', sigma, 'g', g, 'nu', nu, ...
                  'initial_contact_points', 0);
numerical = struct('harmonics_qtt', 90, 'version', 3, 'order', 1, ...
                   'simulation_time', inf);
options = struct('live_plotting', false, 'debug_flag', false, ...
                 'folder', fullfile(here, 'lwdr_out'), 'prefix', 'probe', ...
                 'save_results', false, 'optimize_for_bounce', true, ...
                 'saving_frequency', 2e-7);

t_start = tic;
[rc, rt, PC] = solve_motion_v2(physical, numerical, options);
fprintf('solve_motion_v2 finished in %.1f s, %d recorded states\n', toc(t_start), numel(rc));

% ---- flatten to CSV -------------------------------------------------------------
n = numel(rc);
M = rc{end}.nb_harmonics;
scal = zeros(n, 5);
pres = zeros(n, M + 1);
defo = zeros(n, M);
for i = 1:n
    scal(i, :) = [rc{i}.current_time, rc{i}.center_of_mass, ...
                  rc{i}.center_of_mass_velocity, rc{i}.contact_points, rc{i}.dt];
    pres(i, :) = reshape(rc{i}.pressure_amplitudes, 1, M + 1);
    defo(i, :) = reshape(rc{i}.deformation_amplitudes, 1, M);
end
outdir = fullfile(here, 'lwdr_probe');
if ~exist(outdir, 'dir'); mkdir(outdir); end
writematrix(scal, fullfile(outdir, 'scalars.csv'));       % t, com, comvel, cp, dt
writematrix(pres, fullfile(outdir, 'pressure.csv'));      % l = 0..M
writematrix(defo, fullfile(outdir, 'deformation.csv'));   % l = 1..M
writematrix(PC.theta_vector(:), fullfile(outdir, 'theta.csv'));
fid = fopen(fullfile(outdir, 'meta.txt'), 'w');
fprintf(fid, 'We %.10g\nOh %.10g\nBo %.10g\nM %d\nnangles %d\nFr %.10g\n', ...
        We, Oh, Bo, M, numel(PC.theta_vector), PC.froude_nb);
fclose(fid);
fprintf('wrote %s\n', outdir);
