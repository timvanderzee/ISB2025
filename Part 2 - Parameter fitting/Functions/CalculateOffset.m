function [offset] = CalculateOffset(info, params_subject, t_span)
%Calculate offset voor bepalen LMT en MA
% x (BK) + offset = knee angle die nodig is voor lMT en ma te bepalen

% Inverse kinematics data 
IK         = importdata([info.path,'IK_Trial0',num2str(info.trial),'.mot']); 

% Selection ~ pendulum onset till end
start      = find(IK.data(:,1)== t_span(1)); 
stop       = find(IK.data(:,1)== t_span(2));
IK_sel     = IK.data(start:stop,:); 

% Find pelvis tilt angle 
col_pelT   = find(strcmp(IK.colheaders,'pelvis_tilt'));
pelv_tilt  = IK_sel(:,col_pelT); 

% Find hip flexion angle
if params_subject.leg == 18
    col_hipF   = find(strcmp(IK.colheaders,'hip_flexion_l'));
else
    col_hipF   = find(strcmp(IK.colheaders,'hip_flexion_r'));
end
hip_flex  = IK_sel(:,col_hipF); 

% Interpolate IK data to splined experimental data 
[pelv_tilt_spline, ~,~,~] = ExpDatAtDiscrTime(t_span,IK_sel(:,1),pelv_tilt,0.001); %0.002
[hip_flex_spline, ~,~,~]  = ExpDatAtDiscrTime(t_span,IK_sel(:,1),hip_flex,0.001);  %0.002

% Find offset 
start_pos      = ones(length(pelv_tilt_spline),1)*90;
offset         = start_pos - pelv_tilt_spline - hip_flex_spline;

end

