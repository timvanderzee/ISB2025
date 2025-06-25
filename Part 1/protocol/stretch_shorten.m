function[t,x] = stretch_shorten(model, Ts, us, x0, parms, Ca)

odeopt = odeset('maxstep',1e-2);

if nargin < 5
    Ca = parms.Ca;
end

t = 0;
x = x0;
xd0 = zeros(size(x0));

% loop over phases
for i = 1:length(us)
    parms.vmtc = us(i);
    T = Ts(i);

    if T > 0
        
   
        [tnew,xnew] = ode15s(model, [0 T], x0, odeopt, parms, Ca);

        
        x0 = xnew(end,:);

        t = [t; tnew(2:end)+t(end)]; 
        x = [x; xnew(2:end,:)];
    end
end
end