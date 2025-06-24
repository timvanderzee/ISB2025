clear all; close all; clc
mainfolder = 'C:\Users\u0167448\Documents\GitHub\ISB2025';
addpath(genpath(cd))

% stype = 'standing_balance';
stype = 'pendulum_test';

[input, s0, t0, tf, models, muscle_names] = load_data(stype);

% generic muscle parameters
vec_ksrs = 0;
vec_kT = 35;
k = 1;
M = length(muscle_names);
auxdata = get_muscle_parms();
auxdata.NMuscles = M;

%% OpenSim
cd([mainfolder, '\Part 2 - OpenSim\Movement simulation\input\', stype])

% Import the OpenSim modeling classes
import org.opensim.modeling.*

% Read in the osim model
osimModel = Model(models{1}); %

% 3D model: only used to get muscle info
osimModel3D = Model(models{2}); %

params = zeros(5, length(muscle_names));
muscles = osimModel3D.getMuscles();

for i = 1:length(muscle_names)
    muscle = muscles.get(muscle_names{i});
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
if strcmp(stype,'standing_balance')
    numVar = osimState.getNY();
    for i = 0:numVar-1
        osimState.updY().set(i, s0(i+1,1));
    end
end

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
input.lMT = zeros(1,M);
input.vMT = zeros(1,M);

for i = 1:Ncontrols
    currentActuator = Actuators.get(i-1);
    if strcmp(currentActuator.getConcreteClassName(), 'PathActuator')
        shouldBePathAct = osimModel.getActuators().get(currentActuator.getName());
        currentPathAct = PathActuator.safeDownCast(shouldBePathAct);
        input.lMT(i) = currentPathAct.getLength(osimState);
        input.vMT(i) = currentPathAct.getLengtheningSpeed(osimState);
    end
    controls_all(i,1) = cell(currentActuator.getName());
    q(i,1) = cell(currentActuator.getName());
end


%% simulate movement
close all; clc
modelnames = {'Hill', 'Biophysical'};
ls = {'-','--'};

for j = 1:2
    [tFW,sFW] = sim_movement(modelnames{j}, t0, tf, s0, input, osimModel, osimState, auxdata);

    fse = sFW(:,auxdata.NStates+1:auxdata.NStates+auxdata.NMuscles);

    figure(1)
    for m = 1:M
        subplot(3,3,m)
        plot(tFW, fse(:,m)); hold on
        box off
        title(muscle_names{m})
    end
end

return


%% write things
s = [sFW(:,1) sFW(:,7) sFW(:,2) sFW(:,8) sFW(:,3) sFW(:,9) sFW(:,4) sFW(:,10) sFW(:,5) sFW(:,11) sFW(:,6) sFW(:,12)];

cd('C:\Users\u0167448\Documents\GitHub\ISB2025\Part 2 - OpenSim\Movement simulation\output')
% Write model states to an OpenSim STO file
StatesData.name = 'Pendulum_FW_States';
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
        

