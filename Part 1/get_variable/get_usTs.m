function[us, Ts] = get_usTs(v, AMPs, tiso, ISI, parms, type)

if nargin < 6
    type = 'stretch-shorten-stretch';
end

if strcmp(type,'stretch-shorten-stretch') % assume AMPs = n x 2

    % compute times
    T = AMPs./v; % time of stretch (s)
    Ts = [tiso T(1) T(2) ISI T(3) tiso];

    % compute velocities
    vs = [0 v(1) -v(2) 0 v(3) 0];
    us = vs * 0.5 * parms.s / parms.h; % fraction of h
    
else % assume AMPs = n x 6
    % compute times
    T = AMPs./v; % time of stretch (s)
    Ts = [tiso(1) T(1) T(2) ISI(1) T(3) tiso(2) T(4) tiso(3) T(5) ISI(2) T(6) tiso(4)];

    % compute velocities
    vs = [0 v(1) -v(2) 0 v(3) 0 -v(4) 0 -v(5) 0 v(6) 0];
    us = vs * 0.5 * parms.s / parms.h; % fraction of h
end

end