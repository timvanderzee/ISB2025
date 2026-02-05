%% load parameters and set up simulations 
%
% Simulate force-length, force-velocity, force-pCa curves and save them %
% Code by Hansol X. Ryu (hansol.ryu@gatech.edu) for ISB 2025 tutorial % 

% load parameters % 
addpath(genpath([pwd,'/../..']))
warning('off')
load('parms.mat')
load('protocol.mat')
warning('on')

% You can un-comment this part to remove elastic components % 
% parms.forcible_detachment = 0; 
parms.kse = 1000; % series elastic element
% parms.kpe = parms.kpe*2; % parallel elastic element
parms.no_tendon = 0; % tendon
parms.kT = 0.0200;

odeopt = odeset('maxstep',1e-2);
half_s_len_norm = parms.s/2/parms.h;
nbins = 500;

parms.act = 1;
parms.cosa = 1;
parms.Noverlap = 1;


%% run 
hill_properties = struct();
ref_pCa = 6.089; %6.08; % << pCa during F-L, F-v simulation 

pCa = ref_pCa;
Ca = 10^(-pCa+6);

% initialize model
model = @fiber_dynamics;


for jj=1:10
    tic
parms.xi0 = linspace(-15,15,nbins);
parms.xi = parms.xi0;
parms.nbins = nbins;
parms.xss = zeros(1,parms.nbins + 4);

% parms.xss = zeros(1,7);
parms.n_func = @(xi,Q,eps)Q(1)...
    ./(sqrt(2*pi)*(sqrt(max(Q(3)/Q(1)-(Q(2)/Q(1))^2,eps))))*...
    exp(-((xi-(Q(2)/max(Q(1),eps))).^2)...
    /(2*(sqrt(max(Q(3)/Q(1)-(Q(2)/Q(1))^2,eps)))^2));
parms.xss(end-2) = 0.0909;

% set initial length according to the protocol

x0 = parms.xss;
t_sim = 0;
x_sim = x0;
pCa_sim = 0;

% solve ODE to get XB distributions over time
temp_idx = 1;

parms.vmtc = 0; % isometric 
turnoff_time = 0.2;
pCa_initial = 5.9;

Ts = [1 turnoff_time 4-turnoff_time 4];
pCas = [9 pCa_initial pCa 9];

for i = 1:length(Ts)
    T = Ts(i);
    Ca = 10^(-pCas(i)+6);
    [tnew,xnew] = ode15s(model, [0 T], x0, odeopt, parms, Ca);
    x0 = xnew(end,:);
    t_sim = [t_sim; tnew(2:end)+t_sim(end)];
    x_sim = [x_sim; xnew(2:end,:)];
    pCa_sim = [pCa_sim; Ca*ones(size(tnew(2:end)))];
end

% calculate force
F = nan(1,height(x_sim));
for i = 1:height(x_sim)
    [~,F(i)] = model(t_sim(i), x_sim(i,:)', parms, Ca);
end

toc
end
% plot time-length
% subplot(3,1,1)
% plot(t_sim, x_sim(:,end)/half_s_len_norm)
% hold on
% ylabel('l ce')

% plot time-pCa
subplot(2,1,1)
plot(t_sim, pCa_sim)
hold on
ylabel('\Delta Ca')

% plot time-force
subplot(2,1,2)
plot(t_sim, F' + pCa_sim*0.13)
hold on
plot(t_sim, pCa_sim*0.13)
pause(0.1)
F0(k) = F(temp_idx);
ylabel('F')