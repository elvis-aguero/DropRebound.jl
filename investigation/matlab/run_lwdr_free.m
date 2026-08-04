% FREE OSCILLATION, NO CONTACT, at Oh = 0 -- a cross-code check of the dynamics with
% the contact model switched off entirely.
%
% The drop starts high above the substrate with a seeded l = 2 deformation and a slow
% approach velocity, so it never touches within the simulated window. What is left is
% the free surface dynamics, which both codes must get right independently: theirs as
% a modal ODE with Lamb damping, ours as the variational assembly. If the contours
% agree, every remaining disagreement is in the contact model. If they do not, the
% integrator is the problem and the contact debugging has been chasing a symptom.

here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, 'LWDR', 'matlab', '1_code')));

sigma = 20.5; rho = 0.96; Ro = 0.0203;
nu = 0.0;                              % Oh = 0
Bo = 0.0189;
g  = Bo * sigma / (rho * Ro^2);
M  = 90;

time_unit = sqrt(rho * Ro^3 / sigma);

amps = zeros(1, M);
amps(2) = 0.05 * Ro;                   % zeta_2 = 0.05 once divided by length_unit

physical = struct('undisturbed_radius', Ro, 'initial_velocity', -1.0, ...
                  'initial_height', 2.0, ...
                  'initial_amplitudes', amps, ...
                  'rhoS', rho, 'sigmaS', sigma, 'g', g, 'nu', nu, ...
                  'initial_contact_points', 0);
numerical = struct('harmonics_qtt', M, 'version', 3, 'order', 1, ...
                   'simulation_time', 4e-3);
options = struct('live_plotting', false, 'debug_flag', false, ...
                 'folder', fullfile(here, 'lwdr_out'), 'prefix', 'free', ...
                 'save_results', false, 'optimize_for_bounce', false, ...
                 'saving_frequency', 2e-6);

fprintf('time_unit = %.6g s, simulating %.4g time units\n', time_unit, 4e-3/time_unit);
[rc, rt, PC] = solve_motion_v2(physical, numerical, options);
fprintf('done, %d states\n', numel(rc));

n = numel(rc);
scal = zeros(n, 4);
defo = zeros(n, M);
cpv  = zeros(n, 1);
for i = 1:n
    scal(i, :) = [rc{i}.current_time, rc{i}.center_of_mass, ...
                  rc{i}.center_of_mass_velocity, rc{i}.dt];
    defo(i, :) = reshape(rc{i}.deformation_amplitudes, 1, M);
    cpv(i)     = rc{i}.contact_points;
end
outdir = fullfile(here, 'lwdr_free');
if ~exist(outdir, 'dir'); mkdir(outdir); end
writematrix(scal, fullfile(outdir, 'scalars.csv'));
writematrix(defo, fullfile(outdir, 'deformation.csv'));
writematrix(cpv,  fullfile(outdir, 'contact.csv'));
fid = fopen(fullfile(outdir, 'meta.txt'), 'w');
fprintf(fid, 'Ro %.10g\ntime_unit %.10g\nOh %.10g\nBo %.10g\nM %d\nmax_cp %d\n', ...
        Ro, time_unit, nu * sqrt(rho/(sigma*Ro)), 1/PC.froude_nb, M, max(cpv));
fclose(fid);
fprintf('wrote %s ; max contact points = %d\n', outdir, max(cpv));
