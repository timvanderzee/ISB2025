function [x0, N_1] = PendulumPhases(data_exp, bool_plot)
%Define initial state of pendulum and end state of phase 1
%   1. Define initial state
%   2. Define end of phase 1

q_exp  = data_exp.qspline; 
qd_exp = data_exp.qdspline;
N      = data_exp.Nspline;

% Initial state of pendulum
x0 = [q_exp(1) qd_exp(1)]; 

% Define end of first phase (for SRS)
for i = 1:N
    if qd_exp(i)*qd_exp(i+1) < 0 && i > 20
        break;
    end
end

N_1 = i; 

% Plot phase
if bool_plot == 1
    figure()
    plot(q_exp,'k','LineWidth',1.5)
    hold on
    line([i i],[-2.5 0],'Color','r','LineWidth',1.5)
    hold on
    title('Experimental data and end of first swing');
    hold on
    ylabel('Knee angle (rad)')
    legend('Exp data','End of first swing')
    xlabel('Time (frames)')

end

