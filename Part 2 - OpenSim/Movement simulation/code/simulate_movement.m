function [] = simulate_movement(stype, conditions, act, parms, mainfolder)

% get some settings
[input, s0, t0, tf, models, muscle_names] = load_data(stype, act);

% generic muscle parameters
M = length(muscle_names);
auxdata = get_muscle_parms();
auxdata.NMuscles = M;
auxdata.parms = parms;

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
ls = {'-','--'};
color = get(gca,'colororder');
auxdata.fse0 = zeros(1, M);
modelnames = {'Hill', 'Biophysical'};

% loop over models
for j = 1:length(modelnames)
    disp(modelnames{j})

    % loop over conditions
    for i = 1:length(conditions)
        disp(conditions{i})

        if strcmp(stype, 'pendulum_test')
            AMPs = zeros(M, 1);
            if strcmp(conditions{i}, 'pre_movement')
                AMPs(4) = .1;
                AMPs(5) = .2;
            end    
            auxdata.AMP = AMPs;
        end
        
        [tFW,sFW] = sim_movement(modelnames{j}, t0, tf, s0, input, osimModel, osimState, auxdata);
        fse = max(sFW(:,auxdata.NStates+1:auxdata.NStates+auxdata.NMuscles), 0);

        % used the steady-state for next simulations
        auxdata.fse0 = fse(1,:);
        
        % reorganize states
        s = nan(size(sFW,1), auxdata.NStates);
        
        for k = 1:(auxdata.NStates/2)
            kodd = (k-1)*2+1;
            keve = kodd + 1;
            
            s(:,kodd) = sFW(:,k);
            s(:,keve) = sFW(:,k+(auxdata.NStates/2));
        end
            
%         s = [sFW(:,1) sFW(:,7) sFW(:,2) sFW(:,8) sFW(:,3) sFW(:,9) sFW(:,4) sFW(:,10) sFW(:,5) sFW(:,11) sFW(:,6) sFW(:,12)];
        
        if strcmp(stype, 'pendulum_test')
            phi = s(:,1); % knee angle [?]
            plot(tFW, phi, ls{i}, 'color', color(j,:)); hold on
            box off
            title('Knee angle')

%             FSE(j,i) = min(phi);
        else
            for k = 1:3
            subplot(1,3,k)
            phi = sFW(:,k+1); % ankle angle [?]
            plot(tFW, phi, ls{i}, 'color', color(j,:)); hold on
            box off
%             title('Knee angle')
            end
            
        end
        % write to sto file
        if ~isfolder([mainfolder, '\Part 2 - OpenSim\Movement simulation\output\', stype])
            mkdir([mainfolder, '\Part 2 - OpenSim\Movement simulation\output\', stype])
        end
        cd([mainfolder, '\Part 2 - OpenSim\Movement simulation\output\', stype])

        if strcmp(stype, 'pendulum_test')
            if ~isfolder(['a=',num2str(act)])
                mkdir(['a=',num2str(act)])
            end
            cd(['a=',num2str(act)])
        end
        
        % Write model states to an OpenSim STO file
        StatesData.name = [stype, '_', modelnames{j}, '_', conditions{i}];
        StatesData.nRows = length(tFW);
        StatesData.nColumns = Nstates + 1;
        StatesData.labels = [{'time'}; stateNames]';
        StatesData.inDegrees = false;
        StatesData.data = [tFW,s(:,1:Nstates)];
        writeOpenSimStatesFile(StatesData)
        
    end
    disp(' ');
end
end
