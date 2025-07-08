function [ds] = compute_state_derivatives(t, s, input, osimModel, osimState, auxdata, type)

act = input.act;

if isfield(input, 'platacc')
    acc = interp1(input.time, input.platacc, t)*9.81;
else
    acc = [];
end

body_states = s(1:auxdata.NStates);
MStates = length(s(auxdata.NStates+1:end)) / auxdata.NMuscles; % number of muscle states
muscle_states = reshape(s(auxdata.NStates+1:end), auxdata.NMuscles, MStates);
fse = max(muscle_states(:,1), 0);

% Import the OpenSim modeling classes
import org.opensim.modeling.*

% Update model state with current values  
osimState.setTime(t);
numVar = osimState.getNY();
for i = 0:numVar-1
   osimState.updY().set(i, body_states(i+1,1));
end

% Get the names of the controls/muscles from the model (same in this case)
osimModel.computeStateVariableDerivatives(osimState); % To be able to access muscle fiber velocity.

LMT = zeros(1,auxdata.NMuscles);
VMT = zeros(1,auxdata.NMuscles);
Actuators = osimModel.getActuators();
% Path actuators have to be before any other actuators.
for i = 1:auxdata.NMuscles
   currentActuator = Actuators.get(i-1);
   if strcmp(currentActuator.getConcreteClassName(), 'PathActuator')	
       shouldBePathAct = osimModel.getActuators().get(currentActuator.getName());
       currentPathAct = PathActuator.safeDownCast(shouldBePathAct);
       LMT(i) = currentPathAct.getLength(osimState);
       VMT(i) = currentPathAct.getLengtheningSpeed(osimState);
   end
end


dx = zeros(auxdata.NMuscles,MStates);
FT = zeros(auxdata.NMuscles,1);

for i = 1:auxdata.NMuscles
    [dx(i,:), FT(i)] = TendonForceOdeVecSRS(t, muscle_states(i,:),[0; 5],[act(i); act(i)],[LMT(i); LMT(i)], [VMT(i); VMT(i)], auxdata.params(:,i), auxdata.Fvparam, auxdata.Fpparam, auxdata.Faparam, auxdata.kT, type, auxdata.parms);   
end

M = 2000;
controls = [fse; M*acc/10000];

x_dot = computeOpenSimModelXdot(body_states,controls,t,osimModel,osimState);

ds = [x_dot; dx(:)];