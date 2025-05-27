%% main code to run sims for ISB 2025 tutorials

% set sim number and date - to save results
parms.fp_custMyoSim = '/Users/surabhisimha/GitHub/ISB2025/SimhaModel/';
parms.simNo = 1;
parms.date  = 20250527;
parms.generateResults = 1;
parms.plotFigs = 1;
parms.time_step = 0.001;
%% set extrafusal fiber active parameters
sarcE.pCa=[];
sarcE.act=[];
sarcE.initial_act       = 0;
sarcE.activating_act    = 100; 
sarcE.power_stroke      = 2.5;
sarcE.hs_length         = 1300;
sarcE.compliance_factor = 1;
sarcE.bin_min = -20;      % min x value for myosin distributions in nm
sarcE.bin_max = 20;       % max x value for myosin distributions in nm
sarcE.bin_width = 0.5;
sarcE.k_cb = 0.001;
sarcE.thick_filament_length = 815;
sarcE.thin_filament_length = 1120;
sarcE.bare_zone_length = 80;
sarcE.k_falloff = 0.002; 
sarcE.max_rate = 5000;
sarcE.cb_number_density = 6.9e16;

%% extrafusal fiber passive parameteres
sarcE.isTendon           = 0;
sarcE.tendon_stiffness   = 5e3;
sarcE.passive_force_mode = 'linear'; % linear exponential
sarcE.hsl_slack          = 850;
sarcE.k_passive          = 0.0017; % only for linear parallel passive force
sarcE.passive_sigma      = 0.01; % only for exponential parallel passive force
sarcE.passive_L          = 45; % only for exponential parallel passive force

%% call relevant protocol

% run_ramp(parms,sarcE);
run_stretchshorten(parms,sarcE);
