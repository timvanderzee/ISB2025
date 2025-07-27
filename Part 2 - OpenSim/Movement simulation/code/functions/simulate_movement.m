function [FSE] = simulate_movement(stype, conditions, act, phi0, parms, mainfolder, data)

% define output
FSE = [];

% get some settings
[input, s0, t0, tf, models, muscle_names] = load_data(stype, act);

% generic muscle parameters
M = length(muscle_names);
auxdata = get_muscle_parms();
auxdata.NMuscles = M;
auxdata.parms = parms;

%% OpenSim
%cd([mainfolder, '\Part 2 - OpenSim\Movement simulation\input\', stype])
cd([mainfolder, [filesep 'Part 2 - OpenSim' filesep 'Movement simulation' filesep 'input' filesep] , stype])

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
else
    osimState.updY().set(0, phi0(2));
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
ls = {'--','-'};
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
            if strcmp(conditions{i}, 'premovement')
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
                  
        if strcmp(stype, 'pendulum_test')
            color = get(gca,'colororder');
            phi = (s(:,1) + phi0(1)) * 180/pi; % knee angle [?]
            FSE(j,i,2) = max(phi) - min(phi);
            FSE(j,i,1) = max(data(:,2,i)) - min(data(:,2,i));
            
            subplot(121)
            plot(data(:,1,i), data(:,2,i), ls{i}, 'color', [.5 .5 .5]); hold on
            title('Data')

            subplot(122)
            plot(tFW, phi, ls{i}, 'color', color(j,:), 'Displayname', [modelnames{j}, ' - ', conditions{i}]); hold on
            title('Model')

            for k = 1:2
                subplot(1,2,k)            
                box off
                xlabel('Time (s)')
                ylabel('Shank angle (deg)')
                xlim([0 5])
                ylim([-160 0])
            end
   
        else
            
            titles = {'Ankle', 'Knee', 'Hip'};
            for k = 1:3
                subplot(1,3,k)
                phi = sFW(:,k+1);
                plot(tFW, phi * 180/pi, '-', 'color', color(j,:),'Displayname', modelnames{j}); hold on
                box off
                title(titles{k})
                ylabel('Joint angle (deg)')
                xlabel('Time (s)')
            end
            
        end
        % write to sto file
        %if ~isfolder([mainfolder, '\Part 2 - OpenSim\Movement simulation\output\', stype])
        if ~isfolder([mainfolder, [filesep 'Part 2 - OpenSim' filesep 'Movement simulation' filesep 'output' filesep] , stype])
    
            %mkdir([mainfolder, '\Part 2 - OpenSim\Movement simulation\output\', stype])
            mkdir([mainfolder, [filesep 'Part 2 - OpenSim' filesep 'Movement simulation' filesep 'output' filesep] , stype])
        end
        %cd([mainfolder, '\Part 2 - OpenSim\Movement simulation\output\', stype])
        cd([mainfolder, [filesep 'Part 2 - OpenSim' filesep 'Movement simulation' filesep 'output' filesep] , stype])
        
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

legend
legend boxoff
end
