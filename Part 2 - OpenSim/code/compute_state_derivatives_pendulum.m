function [ds] = compute_state_derivatives_pendulum(t, s, input, osimModel, osimState, auxdata)

act = input.act;
% acc = interp1(input.time, input.platacc, t)*9.81;


states = s(1:auxdata.NStates);
fse = s(auxdata.NStates+1:auxdata.NStates+auxdata.NMuscles);

% Import the OpenSim modeling classes
import org.opensim.modeling.*

% Update model state with current values  
osimState.setTime(t);
numVar = osimState.getNY();
for i = 0:numVar-1
   osimState.updY().set(i, states(i+1,1));
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


dfse = zeros(auxdata.NMuscles,1);
FT = zeros(auxdata.NMuscles,1);

for i = 1:auxdata.NMuscles
    [dfse(i), FT(i)] = TendonForceOdeVecSRS(t,fse(i),[0; 2],[act(i); act(i)],[LMT(i); LMT(i)], [VMT(i); VMT(i)], auxdata.params(:,i), auxdata.Fvparam, auxdata.Fpparam, auxdata.Faparam, input.lMtilda_isom(i), auxdata.ksrs, auxdata.kT);
end

% M = 2000;
% M = 0;
controls = fse;

x_dot = computeOpenSimModelXdot(states,controls,t,osimModel,osimState);

ds = [x_dot; dfse];