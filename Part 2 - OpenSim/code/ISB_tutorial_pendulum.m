clear all; close all; clc
% addpath(genpath('C:\Users\timvd\Documents\BIOMUS\input'))
% addpath(genpath('C:\Users\timvd\Documents\BIOMUS\soft\forward_dynamics'))
M = 8;                 % aantal spieren

names = {'bifemlh_r'; 'bifemsh_r'; 'glut_max2_r'; 'rect_fem_r'; 'vas_int_r'; 'med_gas_r'; 'soleus_r'; 'tib_ant_r'};
vec_ksrs = [0 0 0 280 280 280];
vec_kT = [35 70 35 35 70 35];
k = 1;
% trial = 15;

[auxdata] = get_muscle_parms(M);

%% OpenSim
cd('C:\Users\u0167448\Documents\GitHub\ISB2025\Part 2 - OpenSim\input')

% Import the OpenSim modeling classes
import org.opensim.modeling.*

% Read in the osim model
osimModel = Model('leg39_path_actuators.osim'); %

% 3D model: only used to get muscle info
osimModel3D = Model('leg39_welded.osim'); %

params = zeros(5, length(names));
muscles = osimModel3D.getMuscles();

for i = 1:length(names)
    muscle = muscles.get(names{i});
    params(3,i) = muscle.getTendonSlackLength();
    params(2,i) = muscle.getOptimalFiberLength();
    params(1,i) = muscle.getMaxIsometricForce(); % one leg model
    params(4,i) = muscle.getPennationAngleAtOptimalFiberLength();
end

params(5,:) = 10*params(2,:);
auxdata.params = params;

% Initialize the model (this builds the system and initialize the state)
osimState = osimModel.initSystem();

% Update model state with current values
osimState.setTime(0);
% numVar = osimState.getNY();
% for i = 0:numVar-1
%     osimState.updY().set(i, s0(i+1,1));
% end

% Get the number of states, coordinates, muscles and controls from the model;
% in this case the number of controls equals the number of muscles
Nstates       = osimModel.getNumStateVariables();
Ncontrols     = osimModel.getNumControls();
Ncoord        = osimModel.getNumCoordinates();

% Get the names of the states from the model
stateNames = cell(Nstates,1);
for i = 1:Nstates
    stateNames(i,1) = cell(osimModel.getStateVariableNames().getitem(i-1));
    
    s0(i,1) = getY(osimState).get(i-1);
end

auxdata.NStates = Nstates;

% Get the names of the controls/muscles from the model
osimModel.computeStateVariableDerivatives(osimState); % To be able to access muscle fiber velocity.
Actuators = osimModel.getActuators();
% Muscles = osimModel.getMuscles();
controlNames = cell(Ncontrols,1);
controls_all = cell(Ncontrols,1);
lMT = zeros(1,M);
vMT = zeros(1,M);
for i = 1:Ncontrols
    currentActuator = Actuators.get(i-1);
    if strcmp(currentActuator.getConcreteClassName(), 'PathActuator')
        shouldBePathAct = osimModel.getActuators().get(currentActuator.getName());
        currentPathAct = PathActuator.safeDownCast(shouldBePathAct);
        lMT(i) = currentPathAct.getLength(osimState);
        vMT(i) = currentPathAct.getLengtheningSpeed(osimState);
    end
    controls_all(i,1) = cell(currentActuator.getName());
    q(i,1) = cell(currentActuator.getName());
end


%% simulate isometric contraction
lMtilda_isom = zeros(M,1);
fse_isom = zeros(M,1);

t_input = [0; 10];

modelnames = {'TendonForceOdeVecSRS', 'TendonForceOdeVecSRS_BP'};
x0 = {0, [0 0 0 0 0 0]};
x_isom = cell(M,2);
act = .05;

for m = 1:M
    figure(1)
    nexttile
    
    for j = 1:2
        A = [act act];
        LMT = [lMT(m); lMT(m)];
        VMT = [vMT(m); vMT(m)];
        
        % simulate
        [t,state] = ode15s(eval(['@', modelnames{j}]),[0,10],x0{j},[],t_input, A,LMT,VMT,params(:,m), auxdata.Fvparam, auxdata.Fpparam, auxdata.Faparam, lMT(m), 0, vec_kT(k));
        
        fse = state(:,1);
        fse_isom(m,j) = fse(end);
        x_isom{m,j} = state(end,:);
        
        plot(t, fse); hold on
        
        FMo = ones(size(fse,1),1)*params(1,m);
        lMo = ones(size(fse,1),1)*params(2,m);
        lTs = ones(size(fse,1),1)*params(3,m);
        alphao = ones(size(fse,1),1)*params(4,m);
        
        FT = fse .* FMo;
        lTtilda = fse/vec_kT(k) + 1;
        lM = sqrt((lMo.*sin(alphao)).^2+(lMT(m)-lTs.*lTtilda).^2);
        lMtilda = lM./lMo;
        lMtilda_isom(m,j) = lMtilda(end);
    end
    
  
end

return

%% simulate movement
close all
% s0 = [-1.04298e-09 -0.00528183 -0.464606 -0.00116913 0.0557369 -0.0255086 0 0 0 0 0 0];

auxdata.ksrs = vec_ksrs(k);
auxdata.kT = vec_kT(k);
modelnames = {'compute_state_derivatives_pendulum', 'compute_state_derivatives_BP_pendulum'};

input.act = ones(M,1) * act;

% s0 = getQ

for j = 1
    
    x0 = reshape([x_isom{:,j}], size(x_isom{1,j},2), M)';
    
    %             if j == 2
    %                 x0(:,end-1) = x0(:,end-1) * .9;
    %             end
    
    input.lMtilda_isom = lMtilda_isom(:,j);
    [tFW,sFW] = ode15s(eval(['@',modelnames{j}]),[0 5],[s0(:); x0(:)],[], input, osimModel, osimState, auxdata);
    
    fse = sFW(:,auxdata.NStates+1:auxdata.NStates+auxdata.NMuscles);
    
    
    figure(1)
    for m = 1:M
        subplot(3,3,m)
        plot(tFW, fse(:,m)); hold on
        box off
    end

    if j == 2
        Q0 = sFW(:,auxdata.NStates+auxdata.NMuscles+1:(auxdata.NStates+auxdata.NMuscles*2));
        Q2 = sFW(:,(auxdata.NStates+auxdata.NMuscles*2)+1:(auxdata.NStates+auxdata.NMuscles*3));
        Ld = sFW(:,(auxdata.NStates+auxdata.NMuscles*3)+1:(auxdata.NStates+auxdata.NMuscles*4));
        Non = sFW(:,(auxdata.NStates+auxdata.NMuscles*4)+1:(auxdata.NStates+auxdata.NMuscles*5));
        DRX = sFW(:,(auxdata.NStates+auxdata.NMuscles*5)+1:(auxdata.NStates+auxdata.NMuscles*6));
    
        figure(2)
        for m = 1:M
            subplot(3,3,m)
            plot(tFW, DRX(:,m)); hold on
            box off
        end
        
    end
    
    figure(3)
    plot(tFW, sFW(:,1:Nstates)); hold on

end

return


%% write things
s = [sFW(:,1) sFW(:,7) sFW(:,2) sFW(:,8) sFW(:,3) sFW(:,9) sFW(:,4) sFW(:,10) sFW(:,5) sFW(:,11) sFW(:,6) sFW(:,12)];

% Write model states to an OpenSim STO file
StatesData.name = 'Pendulum_FW_States_v3';
StatesData.nRows = length(tFW);
StatesData.nColumns = Nstates + 1;
StatesData.labels = [{'time'}; stateNames]';
StatesData.inDegrees = false;
StatesData.data = [tFW,s(:,1:Nstates)];
writeOpenSimStatesFile(StatesData)

%%
% Create data structure for the controls file
ControlData.name = 'Pendulum_FW_Controls';
ControlData.nRows = size(tFW);
ControlData.nColumns = Ncontrols+1; %All the controls + time
ControlData.inDegrees = false;
ControlData.labels= [{'time'}; controlNames]';
writeOpenSimControlFile(ControlData)
        

%%
return

function[auxdata] = get_muscle_parms(M)

% PARAMETERS
auxdata.NMuscles = M;

load('Fvparam.mat','Fvparam')
Fvparam(1) = 1.475*Fvparam(1);
Fvparam(2) = 0.25*Fvparam(2);
Fvparam(3) = Fvparam(3) + 0.75;
Fvparam(4) = Fvparam(4) - 0.027;
auxdata.Fvparam = Fvparam;

load('Faparam.mat','Faparam')
auxdata.Faparam = Faparam;

e0 = 0.6;
kpe = 4;
t50 = exp(kpe * (0.2 - 0.10e1) / e0);
pp1 = (t50 - 0.10e1);
t7 = exp(kpe);
pp2 = (t7 - 0.10e1);
Fpparam = [pp1;pp2];
auxdata.Fpparam = Fpparam;

end

function[input, s0, t0, tf] = load_data(trial, fs)

%% Initial state from inverse kinematics
load('PlatOnsets.dat','PlatOnsets')
eval(['load Trial' num2str(trial) '_InputDataTime'])
tOnset = PlatOnsets(trial);
time = InputDataTime;

indeces = find(time > tOnset);
sampleBegin = indeces(1) - 1 + round(fs * 0.02); % start simulation 100 ms before perturbation
t0 = InputDataTime(sampleBegin);
sampleEnd = indeces(1) - 1 + round(fs * 0.12); % end simulation 50 ms after platform onset
tf = InputDataTime(sampleEnd);

IKsolution = importdata(['SCP_gait10dof18musc_Trial' num2str(trial) '_ks.mot']);

% INPUT
pelvis_tilt_exp = IKsolution.data(sampleBegin,2);
pelvis_tilt_mod = IKsolution.data(sampleBegin,7) + IKsolution.data(sampleBegin,6) +IKsolution.data(sampleBegin,5);
delta_pelvis_tilt = pelvis_tilt_exp + pelvis_tilt_mod;

s0 = zeros(8,1);
s0(2:4) = [-IKsolution.data(sampleBegin,7)+delta_pelvis_tilt IKsolution.data(sampleBegin,6) -IKsolution.data(sampleBegin,5)]/180*pi;

%% Force data and muscle excitation
load(['Trial' num2str(trial)], 'data');

platpos = data.analog.platemotion(:,2)/100;
input.platacc = data.analog.platemotion(:,6);
input.time = [1:1:length(data.analog.platemotion(:,6))]/data.analog.samplerate;

% Muscle excitation
load('musdyn_SRS_dlmt_Trial15_ebase3_min10_50ms_ksrs0_web30_rev.mat', 'solution')
input.act =  max(solution.parameter(1:9), 0.01);

end

