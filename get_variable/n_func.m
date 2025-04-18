function[n] = n_func(Q, xi)

    Q0 = Q(1);
    Q1 = Q(2);
    Q2 = Q(3);

    eps = 1e-6;

    Q0c = max(Q0, eps); % correct because can't be zero

    p = Q1/Q0c; % Eq. 52
    q = sqrt(max(Q2/Q0c - (Q1/Q0c)^2, eps));  % Eq. 52

    % n
    n = Q0 ./ (sqrt(2*pi)*q) * exp(-((xi-p).^2) / (2*q^2));  % Eq. 52, modified
end