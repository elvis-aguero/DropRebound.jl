% Do the ancestor's equations survive We = 5 and 10 at low Ohnesorge, where mine
% return |zeta| ~ 10 and never release? Same dimensionless groups, its own solver.
here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, 'LWDR', 'matlab', '1_code')));

sigma = 20.5; rho = 0.96; Ro = 0.0203; Bo = 0.0189;
g = Bo*sigma/(rho*Ro^2);
T_CAP = sqrt(rho*Ro^3/sigma);
fprintf('%-6s %-6s | %-8s %-8s %-8s %-9s %s\n','Oh','We','released','t_c','maxcp','max|zeta|','states');
for Oh = [0.03 0.30]
  nu = Oh*sqrt(sigma*Ro/rho);
  for We = [5.0 10.0]
    v0 = -sqrt(We*sigma/(rho*Ro));
    physical = struct('undisturbed_radius',Ro,'initial_velocity',v0,'rhoS',rho, ...
                      'sigmaS',sigma,'g',g,'nu',nu,'initial_contact_points',0);
    numerical = struct('harmonics_qtt',90,'version',3,'order',1,'simulation_time',20e-3);
    options = struct('live_plotting',false,'debug_flag',false, ...
                     'folder',fullfile(here,'lwdr_out'),'prefix','hw', ...
                     'save_results',false,'optimize_for_bounce',false, ...
                     'saving_frequency',2e-6);
    try
      [rc,~,PC] = solve_motion_v2(physical, numerical, options);
      n = numel(rc); cps = zeros(n,1); mz = 0;
      for i=1:n
        cps(i) = rc{i}.contact_points;
        mz = max(mz, max(abs(rc{i}.deformation_amplitudes))/Ro);
      end
      inc = find(cps>0,1); lastc = find(cps>0,1,'last');
      rel = ~isempty(lastc) && lastc < n && all(cps(lastc+1:end)==0);
      tc = 0;
      if ~isempty(inc) && ~isempty(lastc)
        tc = (rc{lastc}.current_time - rc{inc}.current_time)/T_CAP;
      end
      fprintf('%-6.2f %-6.1f | %-8s %-8.3f %-8d %-9.3f %d\n', Oh, We, ...
              string(rel), tc, max(cps), mz, n);
    catch ME
      fprintf('%-6.2f %-6.1f | ERROR: %s\n', Oh, We, ME.message);
    end
  end
end
