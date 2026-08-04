% Does the ancestor's contact count chatter on the same case where mine swings +-4?
here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here,'LWDR','matlab','1_code')));
sigma=20.5; rho=0.96; Ro=0.0203; Bo=0.0189; g=Bo*sigma/(rho*Ro^2);
Oh=0.03; We=10.0; nu=Oh*sqrt(sigma*Ro/rho); v0=-sqrt(We*sigma/(rho*Ro));
physical=struct('undisturbed_radius',Ro,'initial_velocity',v0,'rhoS',rho, ...
                'sigmaS',sigma,'g',g,'nu',nu,'initial_contact_points',0);
numerical=struct('harmonics_qtt',90,'version',3,'order',1,'simulation_time',6e-3);
options=struct('live_plotting',false,'debug_flag',false,'folder',fullfile(here,'lwdr_out'), ...
               'prefix','cq','save_results',false,'optimize_for_bounce',false, ...
               'saving_frequency',2e-6);
[rc,~,~]=solve_motion_v2(physical,numerical,options);
n=numel(rc); cps=zeros(n,1);
for i=1:n; cps(i)=rc{i}.contact_points; end
inc=sum(diff(cps)>0); dec=sum(diff(cps)<0);
big=sum(abs(diff(cps))>1);
fprintf('ANCESTOR cp: max=%d  increases=%d  decreases=%d  jumps>1 node=%d  states=%d\n', ...
        max(cps), inc, dec, big, n);
k=find(cps>0,1);
fprintf('first 30 from onset: %s\n', mat2str(cps(k:min(k+29,n))'));
