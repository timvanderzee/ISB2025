%% Assignment 0: Preparation
clear all; close all; clc

% Add github folder to path
filename = 'main_part2.m';
fullfilepath = which(filename);
filefolder = fullfilepath(1:end-length(filename));
cd(filefolder)
cd ..
cd ..
cd ..
mainfolder = cd;
addpath(genpath(mainfolder));

% Import casadi libraries
import casadi.*; 

%% Assignment 1.1: simulate standing balance
% load some default parameters to start with
cd([mainfolder, '\Part 2 - OpenSim\Movement simulation\input\common'])
load('parms.mat','parms')

% simulate standing balance
if ishandle(1), close(1); end; figure(1)
simulate_movement('standing_balance', {'reguluar'}, [], [], parms, mainfolder);
set(gcf,'units','normalized','position',[.2 .4 .6 .3])

%% Assignment 1.2: visualize effect of rate parameters
% visualize the effect of cross-bridge rate parameters on attachment and
% detachment functions
if ishandle(2), close(2); end; figure(2)
effect_of_rate_parms(parms)
set(gcf,'units','normalized','position',[.2 .1 .4 .7])

%% Assignment 2.1: fit model parameters
% load some default parameters to start with
cd([mainfolder, '\Part 2 - OpenSim\Movement simulation\input\common'])
load('new_parms.mat','parms')

% load short-range stiffness data (skinned rat soleus muscle fibers), 
% see Horslen et al. 2023: https://doi.org/10.1242/jeb.245456
cd([mainfolder, '\Part 2 - OpenSim\Movement simulation\input\common'])
load('thix_data.mat','SRSdata')

% initialise opti structure
opti = casadi.Opti(); 

% define weigth vector
w1 = 10;    % weight for fitting force-velocity
w2 = 100;   % weight for fitting short-range stiffness
w3 = 1; 	% weight for regularization
w = [w1 w2 w3];

% specify biophysical parameters to be fitted
optparms = {'f', 'k11', 'k22', 'k21'};

if ishandle(3), close(3); end; figure(3)
[newparms, out] = fit_model_parameters(opti, optparms, w, SRSdata, parms);
set(gcf,'units','normalized','position',[.2 .2 .4 .6])

%% Assignment 2.2: evaluate model force-velocity and short-range stiffness properties
if ishandle(4), close(4); end; figure(4)
eval_model_behavior(newparms, out, SRSdata)
set(gcf,'units','normalized','position',[.2 .4 .6 .3])

%% Assignment 2.3: compare the rate functions
if ishandle(5), close(5); end; figure(5)
compare_rate_funcs(newparms, parms)
set(gcf,'units','normalized','position',[.4 .4 .2 .2])

%% Assignment 3.1: simulate standing balance
if ishandle(6), close(6); end; figure(6)
simulate_movement('standing_balance', {'reguluar'}, [], [], newparms, mainfolder);
set(gcf,'units','normalized','position',[.2 .2 .6 .3])

%% Assignment 4.1: simulate pendulum test - typically developing (TD) child
% load TD data
cd([mainfolder, '\Part 2 - OpenSim\Movement simulation\input\pendulum_test'])
load('TD_data.mat','data')

% specify baseline activation of all muscles
act = 2e-2; 

% femur angle wrt horizontal (downward = negative). note: if changed from 
% default value (-.25), need to modify .osim file (see Assignment 7)
phi_femur = -.25; % [rad]

% initial knee angle (full extension = 0, flexion = negative). note: can be
% changed from default value, without needing to modify .osim file
phi_knee = 0; % [rad]

if ishandle(7), close(7); end; figure(7)
FSE_TD = simulate_movement('pendulum_test', {'reguluar', 'premovement'}, act, [phi_femur phi_knee], newparms, mainfolder, data);
set(gcf,'units','normalized','position',[.2 .4 .6 .4])

%% Assignment 4.2: simulate pendulum test - child with cerebral palsy (CP)
% load CP data
cd([mainfolder, '\Part 2 - OpenSim\Movement simulation\input\pendulum_test'])
load('CP_data.mat','data')

% specify baseline activation of all muscles
act = 2e-2;

% femur angle wrt horizontal (downward = negative). note: if changed from 
% default value (-.25), need to modify .osim file (see Assignment 7)
phi_femur = -.25; % [rad]

% initial knee angle (full extension = 0, flexion = negative). note: can be
% changed from default value, without needing to modify .osim file
phi_knee = 0; % [rad]

if ishandle(8), close(8); end; figure(8)
FSE_CP = simulate_movement('pendulum_test', {'reguluar', 'premovement'}, act, [phi_femur phi_knee], newparms, mainfolder, data);
set(gcf,'units','normalized','position',[.2 .2 .6 .4])

%% Assignment 6: visualize first swing excursion and its change with premovement
if ishandle(9), close(9); end; figure(9)
color = get(gca,'colororder');

% determine change in first swing excursion
dFSE_CP = FSE_CP(:,2,:) - FSE_CP(:,1,:);
dFSE_TD = FSE_TD(:,2,:) - FSE_TD(:,1,:);

% data
plot(FSE_CP(1,1,1), dFSE_CP(1,:,1), 'o', 'color', [.5 .5 .5], 'markerfacecolor', [.5 .5 .5]); hold on
plot(FSE_TD(1,1,1), dFSE_TD(1,:,1), 's', 'color', [.5 .5 .5], 'markerfacecolor', [.5 .5 .5])

% Hill-type model
plot(FSE_CP(1,1,2), dFSE_CP(1,:,2), 'o', 'color', color(1,:), 'markerfacecolor', color(1,:)); hold on
plot(FSE_TD(1,1,2), dFSE_TD(1,:,2), 's', 'color', color(1,:), 'markerfacecolor', color(1,:))

% Biophysical model
plot(FSE_CP(2,1,2), dFSE_CP(2,:,2),  'o', 'color', color(2,:), 'markerfacecolor', color(2,:)); hold on
plot(FSE_TD(2,1,2), dFSE_TD(2,:,2),  's', 'color', color(2,:), 'markerfacecolor', color(2,:))

% connecting lines
plot([FSE_CP(1,1,1) FSE_TD(1,1,1)], [dFSE_CP(1,:,1) dFSE_TD(1,:,1)], '--', 'color', [.5 .5 .5])
plot([FSE_CP(1,1,2) FSE_TD(1,1,2)], [dFSE_CP(1,:,2) dFSE_TD(1,:,2)], '--', 'color', color(1,:))
plot([FSE_CP(2,1,2) FSE_TD(2,1,2)], [dFSE_CP(2,:,2) dFSE_TD(2,:,2)], '--', 'color', color(2,:))

xlabel('First swing excursion (deg)')
ylabel('Change in first swing excursion (deg)')
set(gca,'Xdir','reverse')
box off
legend('Data - CP', 'Data - TD', 'Hill - CP', 'Hill - TD', 'Biophysical - CP', 'Biophysical - TD', 'location', 'best')
legend boxoff
set(gcf,'units','normalized','position',[.3 .3 .3 .4])
