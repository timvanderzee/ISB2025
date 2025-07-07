clear all; close all; clc

fullfilepath = which('main_part2.m');
filefolder = fullfilepath(1:end-12);

cd(filefolder)
cd ..
cd ..
cd ..

mainfolder = cd;
addpath(genpath(mainfolder));
import casadi.*;        % Import casadi libraries

%% load some parameters
cd('Part 2 - OpenSim\Movement simulation\input\common')
load('new_parms.mat','parms')

%% fit model parameters
opti = casadi.Opti();   % Initialise opti structure
vmax = 5; % maximal shortening velocity (L0/s)
RT = 0.1; % recovery time (s)
SRS_rel = 0.7; % relative SRS compared to no pre-movement (-) for the above recovery time
V_rel = 0.5; % relative velocity at which SRS changes are tested (relative to vmax) 
w = [10 100 1]; % weights for fitting (1) force, (2) history-dependence, (3) regularization

% parameters that are to be fitted
optparms = {'f', 'k11', 'k22', 'k21', 'JF'};

% running the fitting function
newparms = fit_model_parameters(opti, optparms, w, vmax, RT, SRS_rel, V_rel, parms);

%% simulate movement
close all

act = 1e-1; % activation
[FSE] = simulate_movement(act, newparms, mainfolder);

%% first swing excursion
figure(4)
bar(categorical({'Hill', 'Biophysical'}), -FSE)
legend('No pre-movement', 'With pre-movement', 'location', 'bestoutside')
legend boxoff
ylabel('First swing excursion (rad)')
box off