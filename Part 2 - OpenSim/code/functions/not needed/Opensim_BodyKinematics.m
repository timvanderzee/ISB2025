function [out] = Opensim_COM_kinematics(filename_model,results_dir,event,kinematics_file)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% problemen

% De files zijn te groot om de analyse succesvol te laten lopen
% -> oplossing = inkorten files

% Er is nog altijd een probleem met memory leaks
%   => via command line laten lopen

%%



import org.opensim.modeling.*
function_path = which('settings_COM_kinematics.xml');
% get the analysis tool and change variables
path_generic_file=function_path;

tool=AnalyzeTool(path_generic_file);
osimModel=Model(filename_model);
tool.setModel(osimModel);
tool.setResultsDir(results_dir);
tool.setInitialTime(event(1));
tool.setFinalTime(event(2));
[path name ext]=fileparts(kinematics_file);
tool.setName(name);

% run the analysis
tool.setCoordinatesFileName(kinematics_file);
out_path_xml=fullfile(results_dir,['analysis_settings_' name '.xml']);
tool.setModelFilename(filename_model);
tool.print(out_path_xml);
%tool.run();
exe_path='"C:\OpenSim 3.2\bin\analyze.exe"';
Command = [exe_path ' -S ' out_path_xml];
system(Command)


out.body_kin.acc=importdata(fullfile(results_dir,[name '_BodyKinematics_acc_global.sto']));
out.body_kin.vel=importdata(fullfile(results_dir,[name '_BodyKinematics_vel_global.sto']));
out.body_kin.pos=importdata(fullfile(results_dir,[name '_BodyKinematics_pos_global.sto']));

clear tool


end

