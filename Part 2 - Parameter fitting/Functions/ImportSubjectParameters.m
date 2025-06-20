function [params] = ImportSubjectParameters(info)

% subject
subjectname = info.subj; 

% Path to parameters 
% pathmain = pwd;
% [pathTemp,~,~] = fileparts(pathmain);
info           = load([info.path, subjectname, '.mat']);

% Define parameters
params.m_tot   = info.(char(subjectname)).m;               % Totall mass of patient
params.age     = info.(char(subjectname)).age;             % Age of subject
params.leg     = info.(char(subjectname)).z;               % 18 if left leg, 11 if right leg
params.Nmr     = info.(char(subjectname)).Nmr;             % First trial with MR
params.g       = 9.81;                              % Gravitational constant
%params.length_tibia = info.(char(subjectname)).lTibia; 

end
