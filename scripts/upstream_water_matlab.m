% Their solver, our water conditions, our Weber grid.
%
% cgs, matching the Julia sweep: rho = 1.0 g/cc, sigma = 72.8 dyn/cm,
% R = 0.03 cm, eta = 1e-3 Pa s -> nu = 0.01 St, g = 981 -> Bo = 0.0121.
% Oh = nu*sqrt(rho/(sigma*R)) = 0.00677, matching ours.
addpath(fullfile(pwd, 'simulation_code'), '-begin');

rho = 1.0; sig = 72.8; R = 0.03; nu = 0.01; g = 981.0;
Oh = nu*sqrt(rho/(sig*R));
Bo = rho*g*R^2/sig;
fprintf('Oh = %.5f   Bo = %.5f\n', Oh, Bo);

Wes = [0.0077 0.0161 0.0335 0.0698 0.1453 0.3028 0.6309 1.3144 2.7386];
fid = fopen('/Users/harrislab/Documents/GitHub/DropRebound.jl/outputs/csv/water_upstream.csv','w');
for k = 1:numel(Wes)
    We = Wes(k);
    v0 = -sqrt(We*sig/(rho*R));           % cgs, negative = towards impact
    phys = struct('undisturbed_radius', R, 'initial_velocity', v0, ...
                  'rhoS', rho, 'sigmaS', sig, 'g', g, 'nu', nu);
    num  = struct('harmonics_qtt', 90, 'version', 3, 'simulation_time', inf);
    opts = struct('live_plotting', false, 'save_results', false, ...
                  'optimize_for_bounce', true, 'debug_flag', false);
    t0 = tic;
    try
        [rc, rt, PC] = solve_motion_v2(phys, num, opts);
        % recorded_conditions is a CELL array of structs, not a struct array
        vel = cellfun(@(s) s.center_of_mass_velocity, rc);
        cps = cellfun(@(s) s.contact_points, rc);
        idx = find(cps > 0);
        if isempty(idx)
            cor = NaN; tc = NaN;
        else
            i1 = idx(1); i2 = idx(end);
            cor = abs(vel(i2)/vel(i1));
            tc  = rt(i2) - rt(i1);
        end
        fprintf('We=%.4f  cor=%.5f  tc=%.4f  %.0fs\n', We, cor, tc, toc(t0));
        fprintf(fid, '%g,%g,%g,%g\n', We, cor, tc, Oh);
    catch ME
        fprintf('We=%.4f  FAILED: %s\n', We, ME.message);
        fprintf(fid, '%g,NaN,NaN,%g\n', We, Oh);
    end
end
fclose(fid);
