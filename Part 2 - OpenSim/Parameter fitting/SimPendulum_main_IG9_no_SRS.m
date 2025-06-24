%% SimPenulum_main
% Input
% Change here subject and trial
info.subj   = 'TD5';           % Subject name
info.trial  = 1; % Trial number
info.option = 'Opt14_Activation_TanH_tresh_inertie20_IG9_ksrs200_wq1_end_kR0';   % Name to save results
info.wq     = 1;               % weight on q error
info.wqd    = 0.5;             % weight on qd error
info.kSRS   = 200; 
info.tau    = 0.080; 

% IG
IG.a_ext    = 0.1; 
IG.a_flex   = 0.05; 
IG.kFpe     = 0.2;
IG.B        = 0.03; 
IG.dt1      = 0.001; 
IG.kR       = 0.1; 

%% Import subject parameters and experimental data
% Path info -  Path to model and experimental data
% pathmain = pwd;
% [pathTemp,~,~] = fileparts(pathmain);
pathTemp = pwd;
% [pathRepo,~,~] = fileparts(pathTemp);
info.path      = [pathTemp '\Experimental data\' info.subj '\'];

% Import subject parameters
params_subject = ImportSubjectParameters(info);    % Input parameters (mtot,lc, l, age, m, RG, SE, Nmr, z)

% Import experimental data
bool_plot = 1;
dt_spline = 0.002; %0.0021;
data_exp  = ImportExperimentalData(info, bool_plot, params_subject, dt_spline);

% Define phases of pendulum (initial state, end of first swing)
bool_plot = 1;
[data_exp.x0, data_exp.N_1] = PendulumPhases(data_exp, bool_plot);

%% Import OpenSim model parameters
% OpenSim
addpath('MuscleModel');
muscles = {'rect_fem_','bifemlh_'};
[params_OS] = ReadOpenSimParams_OS4(info, params_subject, muscles);

%% Calculate LMT en Ma coefficients
bool_plot = 1;
[coeff_LMT_ma] = DefineLMTCoefficients(params_subject, info, muscles, bool_plot);

%% Create initial guess
bool_guess  = 1; % to create trial specific initial guess of lMtilda
[InitGuess] = CreateInitialGuess(bool_guess,params_OS, data_exp, muscles, coeff_LMT_ma);

%% Define states, controls, bounds and initial guess
addpath(genpath('C:\Users\timvd\Documents\casadi-windows-matlabR2016a-v3.5.5'))
import casadi.*;        % Import casadi libraries
opti = casadi.Opti();   % Initialise opti structure

N       = data_exp.Nspline;
N_1     = data_exp.N_1; 

% States
x            = opti.variable(1,N);
xd           = opti.variable(1,N); 
xdd          = opti.variable(1,N); 
lMtilda      = opti.variable(2,N);  

% (Slack) controls
vMtilda      = opti.variable(2,N);
lM_projected = opti.variable(2,N);
dt1          = opti.variable(1); 

% Parameters that will be optimized
a_ext         = opti.variable(1);
a_flex        = opti.variable(1); 
kFpe          = opti.variable(1);
B             = opti.variable(1); 
I_opt         = opti.variable(1); 

% Bounds
[Ub, Lb] = SelectBounds_Set1(data_exp);
opti.subject_to(Lb.x   <  x     < Ub.x);

opti.subject_to(Lb.a   < a_ext  < Ub.a);
opti.subject_to(Lb.a   < a_flex < Ub.a);
opti.subject_to(Lb.kFpe <  kFpe < Ub.kFpe);    % nominal value = 0.1
opti.subject_to(Lb.B    < B     < Ub.B); 
opti.subject_to(Lb.lM_projected  <  lM_projected(1,:));
opti.subject_to(Lb.lM_projected  <  lM_projected(2,:));
opti.subject_to(Lb.vMtilda  <  vMtilda  < Ub.vMtilda);
opti.subject_to(Lb.lMtilda  <  lMtilda  < Ub.lMtilda);
opti.subject_to(Lb.dt1 < dt1 < Ub.dt1); 
opti.subject_to(Lb.x_end < x(end) < Ub.x_end); % Constraint on resting angle (3 graden + en 3 graden -) 
opti.subject_to(params_OS.inert.I_OS - 0.2*params_OS.inert.I_OS < I_opt < params_OS.inert.I_OS + 0.2*params_OS.inert.I_OS); 

% Constraints on initial states
opti.subject_to(x(1)     == data_exp.x0(1));
opti.subject_to(xd(1)    == 0);

% Initial guess
opti.set_initial(x, data_exp.qspline);
opti.set_initial(xd,data_exp.qdspline);
opti.set_initial(a_ext, IG.a_ext );%0.01
opti.set_initial(a_flex, IG.a_flex  ); %0.01 guess(l)
opti.set_initial(kFpe, IG.kFpe);
opti.set_initial(B, IG.B); % 0.1 
opti.set_initial(lM_projected, InitGuess.lM_projected');
opti.set_initial(lMtilda, InitGuess.lMtilda');     
opti.set_initial(vMtilda, InitGuess.vM);
opti.set_initial(dt1,IG.dt1); 
opti.set_initial(I_opt, params_OS.inert.I_OS);

%% Define problem (muscle model)
% Calculate shift
kT     = 35;
shift  = getshift(kT);

% Calculate dlMdt
dlMdt_ext = vMtilda(1,:); 
dlMdt_flex = vMtilda(2,:); 

% Time phase 1
tF1 = N_1*dt1; 
dt2 = 0.002; 

opti.subject_to(tF1 > 0.2); 
opti.subject_to(xd(1:N_1) < 1e-4); 
opti.subject_to(1e-4 < xd(N_1+1)); % of kleiner als 1 er nog bij? 

% Skeletal dynamics fase 1 
[error_f1] = CalculateMusculoSkeletalDynamics_F1_no_SRS(x(1:N_1),xd(1:N_1),xdd(1:N_1), lMtilda(:,1:N_1), lM_projected(:,1:N_1),kFpe,vMtilda(:,1:N_1), ...
                a_ext, a_flex, data_exp, coeff_LMT_ma, params_OS, shift, B, info, I_opt); 

% Skeletal dynamics fase 2
[error_f2] = CalculateMusculoSkeletalDynamics_F1_no_SRS(x(N_1+1:end),xd(N_1+1:end),xdd(N_1+1:end),lMtilda(:,N_1+1:end), lM_projected(:,N_1+1:end), kFpe, vMtilda(:,N_1+1:end), ...
             a_ext, a_flex, data_exp, coeff_LMT_ma, params_OS, shift, B, info, I_opt); 

% Constraints - trapezoidal integration 
% dt = dt_spline; 
lMtilda_ext = lMtilda(1,:); 
lMtilda_flex= lMtilda(2,:); 
opti.subject_to((xd(1:N_1-1) + xd(2:N_1))*dt1/2 + x(1:N_1-1) == x(2:N_1));
opti.subject_to((xd(N_1:N-1) + xd(N_1+1:N))*dt2/2 + x(N_1:N-1) == x(N_1+1:N));
opti.subject_to((xdd(2:N_1)+xdd(1:N_1-1))*dt1/2 +xd(1:N_1-1) == xd(2:N_1));
opti.subject_to((xdd(N_1+1:N)+xdd(N_1:N-1))*dt2/2 +xd(N_1:N-1) == xd(1+N_1:N));

opti.subject_to((dlMdt_ext(2:N_1)+dlMdt_ext(1:N_1-1))*dt1/2 + lMtilda_ext(1:N_1-1) == lMtilda_ext(2:N_1));
opti.subject_to((dlMdt_ext(1+N_1:N)+dlMdt_ext(N_1:N-1))*dt2/2 + lMtilda_ext(N_1:N-1) == lMtilda_ext(1+N_1:N));
opti.subject_to((dlMdt_flex(2:N_1)+dlMdt_flex(1:N_1-1))*dt1/2 + lMtilda_flex(1:N_1-1) == lMtilda_flex(2:N_1));
opti.subject_to((dlMdt_flex(1+N_1:N)+dlMdt_flex(N_1:N-1))*dt2/2 + lMtilda_flex(N_1:N-1) == lMtilda_flex(1+N_1:N));

opti.subject_to(error_f1 == 0);
opti.subject_to(error_f2 == 0);

% Objective function
%J = DefineObjectiveFunction(x,xd,data_exp, info, vMtilda, a_ext, a_flex, kR); 
kR = 0;
J = DefineObjectiveFunction_time_Lars(x,xd,data_exp, info, vMtilda, a_ext, a_flex, kR, dt1, N_1, dt2, N);
opti.minimize(J); 
    
%% Solve problem
% options for IPOPT
options.ipopt.tol = 1*10^(-6);          
options.ipopt.linear_solver = 'mumps';
          
% Solve the OCP
opti.solver('ipopt',options);
sol = opti.solve();  
    
%% Results 
% Experimental data 
R.exp    = data_exp;
R.info   = info; 
R.subject= params_subject;
R.OS     = params_OS;
R.IG     = IG;

% States
R.x         = sol.value(x); 
R.xd        = sol.value(xd); 
R.xdd       = sol.value(xdd); 
R.lMtilda   = sol.value(lMtilda); 

% Controls
R.vMtilda     = sol.value(vMtilda);
R.lMprojected = sol.value(lM_projected);
R.dlMdt       = sol.value(vMtilda); 

% Parameters
R.a_ext   = sol.value(a_ext);
R.a_flex  = sol.value(a_flex);
R.a       = [R.a_ext R.a_flex]; 
R.kFpe    = sol.value(kFpe);
R.B       = sol.value(B); 
R.dt1     = sol.value(dt1); 
R.dt2     = dt2; 
R.tF1     = N_1*R.dt1; 
R.I_opt   = sol.value(I_opt);

% Objective function
R.wq      = info.wq;
R.wqd     = info.wqd;
R.J       = sol.value(J);

% Additional information
R.error_f1  = sol.value(error_f1);
R.error_f2  = sol.value(error_f2);
R.coeff     = coeff_LMT_ma; 
R.initGuess = InitGuess; 
R.shift     = shift; 

% Bounds
R.bounds.Ub = Ub; 
R.bounds.Lb = Lb; 

% Calculated Parameters
% [R.C] = RecalculateOutcomes2(R, info, N_1); 

% Write results 
% save([pathTemp,'/Results/',info.subj,'_T',num2str(info.trial),'_',info.option,'.mat'],'R')

%% Plot
close all
% Results
h = PlotResults(R, info);
% saveas(h,[pathTemp,'/Results/Figures/', info.subj, '_T',num2str(info.trial),info.option,'_1.Results.fig']);

% Params
g = PlotParams(R,info);
% saveas(g,[pathTemp,'/Results/Figures/', info.subj, '_T',num2str(info.trial),info.option,'_2.Params.fig']);

% Muscle geometry
f = PlotMuscleGeometry(R, info);
% saveas(f,[pathTemp,'/Results/Figures/', info.subj, '_T',num2str(info.trial),info.option,'_3.MuscleGeometry.fig']);

% Torques
p = PlotTorques(R); 
% saveas(p,[pathTemp,'/Results/Figures/', info.subj, '_T',num2str(info.trial),info.option,'_4.Torques.fig']);

% Reflexes
s = PlotReflexes_Force(R);
% saveas(s,[pathTemp,'/Results/Figures/', info.subj, '_T',num2str(info.trial),info.option,'_5.Reflexes.fig']);

%% Run forward simulation
% Params for Ode function
params.kFpe         = R.kFpe; % + 0.0037035; 
params.a_ext        = R.a_ext; %+ 0.0004;
params.a_flex       = R.a_flex; % - 0.010003;
params.data_exp     = R.exp;
params.coeff_LMT_ma = R.coeff;
params.params_OS    = R.OS;
params.shift        = R.shift;
params.B            = R.B; %- 0.007191;
params.info         = R.info;
params.I_opt        = R.I_opt;

time_opt1 = 0:R.dt1:20;

y0  = [R.x(1) R.xd(1) R.lMtilda(1,1) R.lMtilda(2,1)]';
yp0 = [R.xd(1) R.xdd(1) R.dlMdt(1,1) R.dlMdt(2,1)]';

% Ode Options
options = odeset('RelTol',1e-9,'AbsTol',1e-9);

% Simulate
[t1,x1] = ode15i(@(t,y,yp) OdeFun_F1_no_SRS(t,y,yp, params), [0 11], y0, yp0, options);

%%
figure(1)
plot(t1,x1(:,1)*180/pi,'b--')

%% Simulate biophysical model
y0  = [R.x(1) R.xd(1) R.lMtilda(1,1) R.lMtilda(2,1) 1e-4*ones(1,4)]';
yp0 = [R.xd(1) R.xdd(1) R.dlMdt(1,1) R.dlMdt(2,1) zeros(1,4)]';

[t2,x2] = ode15i(@(t,y,yp) OdeFun_F1_Bio(t,y,yp, params), [0 11], y0, yp0, options);


figure(1)
plot(t2,x2(:,1)*180/pi,'g--')



