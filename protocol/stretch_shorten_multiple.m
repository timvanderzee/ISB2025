function[Fmodel, Ts, y, Y, Lmodel] = stretch_shorten_multiple(model, data, parms)

    Amps = data.AMPs;
    X0s = repmat(parms.xss,length(data.pCas),1);

    v = data.v;
    tiso = data.tiso;
    ISI = data.ISI;
    pCais = data.pCas;

    % assume the same across trials
    [us, Ts] = get_usTs(v(1,:), Amps(1,:), tiso(1,:), ISI(1,:), parms);
   
    for i = 1:length(pCais)
        Ca = 10^(-pCais(i)+6);
        [t,x] = stretch_shorten(model, Ts, us, X0s(i,:), parms, Ca);

        F = nan(1,length(x));
        for j = 1:length(x)
            if ~contains(func2str(model), 'implicit') % explicit
                [~,F(j),n(j,:),xi(j,:)] = model(t(j), x(j,:)', parms, Ca);
            else % implicit
                [~,F(j)] = model(t(j), x(j,:)', zeros(size(x(j,:)')), parms, Ca); % note: force does not depend on xdot
            end
        end

        lmtc = x(:,end);

        % get the force and SRS
        y(i) = get_force(t, Ts, F, lmtc, parms);    
    end
    
    if ~isfield(data, 'texp')
        data.texp(:,i) = (0:.001:(sum(Ts(1,:))))' - sum(Ts(1,1:4));
        Y.texp = data.texp;
    end
    
    for i = 1:length(pCais)
        % interpolate model to data
        Fmodel(:,i) = interp1(y(i).ti, y(i).Fi, data.texp(:,1));
        Lmodel(:,i) = interp1(y(i).ti, y(i).lmtc, data.texp(:,1));
        Y.SRS(i) = y(i).SRS;
        Y.Fpre(i) = y(i).Fpre;
%         Fmodel0(i) = y(i).F0;
    end
end