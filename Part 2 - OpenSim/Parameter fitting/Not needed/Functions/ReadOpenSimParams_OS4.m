function [params_OS] = ReadOpenSimParams_OS4(info, params, muscles)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
import org.opensim.modeling.*

% Model 
model_path = [info.path, info.subj,'_ScaledModel_ScaledForces.osim'];    
osimModel  = Model(model_path);
 
% Inertial params
    % Mass
    leg = params.leg; 
    bodies         = osimModel.getBodySet();
    if leg == 18
        tibia      = bodies.get('tibia_l');
    else
        tibia      = bodies.get('tibia_r');
    end

    params_inert.mass_OS = tibia.getMass(); 
    
    % COM
    % com            = ArrayDouble.createVec3(0);
    com            = tibia.getMassCenter;
    params_inert.lc_OS   = abs(com.get(1)); 
    
    %inertia        = Mat33(0);
    inertia         = tibia.get_inertia;
    params_inert.I_OS= inertia.get(0) + params_inert.mass_OS * params_inert.lc_OS ^2; %(apply steiner) 
    params_OS.inert = params_inert; 

% Muscle tendon properties
if leg == 18
    muscles_ind{1} = [char(muscles(1)),'l'];
    muscles_ind{2} = [char(muscles(2)),'l'];
else
    muscles_ind{1} = [char(muscles(1)),'r'];
    muscles_ind{2} = [char(muscles(2)),'r'];
end
   [params_Muscle,lOpt,L_TendonSlack,Fiso,PennationAngle]=ReadMuscleParameters(model_path,muscles_ind);
    params_MT    = params_Muscle; % 5 x 2 matrix: 1 - Fmo, 2 - Lmo, 3 - Lts, 4 - alphao, 5 - vmmax
    params_OS.MT = params_MT; 

% Force- velocity characteristics
load Fvparam.mat
Fvparam(1)       = 1.475*Fvparam(1); Fvparam(2) = 0.25*Fvparam(2);
Fvparam(3)       = Fvparam(3) + 0.75; Fvparam(4) = Fvparam(4) - 0.027;
params_Fvparam   = Fvparam;
params_OS.Fv     = params_Fvparam;

% Force length characteristics (active)
load Faparam.mat
params_Faparam   = Faparam;
params_OS.Fa     = params_Faparam;

% Force length characteristics (passive)
e0  = 0.6; kpe = 4; t50  = exp(kpe * (0.2 - 0.10e1) / e0);
pp1 = (t50 - 0.10e1); t7 = exp(kpe); pp2 = (t7 - 0.10e1);
params.Fpparam = [pp1;pp2];
params_OS.Fp   = params.Fpparam; 
end

