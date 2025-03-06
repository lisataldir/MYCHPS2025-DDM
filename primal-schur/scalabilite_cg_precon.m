% Données du problème
L = 100;
S = 10;
E = 2*1e5;
Fd = 10;

% Plage de valeurs de N à tester
N_values = 2:20;  
iterations = zeros(size(N_values));  % Stocke le nombre d'itérations

for idx = 1:length(N_values)
    N = N_values(idx);
    H = L/N;
    n = 5;  % On garde un nombre d'éléments fixe pour simplifier
    h = H/(n-1);

    % ---------------------------------------------------------------
    % Construction des sous-domaines et matrices locales
    interface = (n-1)*(1:N-1) + 1;
    k0 = E*S/h;
    k = k0 * (2*eye(n) - diag(ones(n-1, 1), 1) - diag(ones(n-1, 1), -1));
    k(1,1) = k0;
    k(n,n) = k0;

    substruct = cell(N,1);
    kbb = cell(N,1);
    kib = cell(N,1);
    kii = cell(N,1);
    Sp = cell(N,1);
    Sd = cell(N,1);
    Rb = cell(N,1);
    A = cell(N,1);
    Ad = cell(N,1);
    bp = cell(N,1);
    num_local_interface = cell(N,1);

    for s = 1:N
        substruct{s} = (s-1)*(n-1) + (1:n);
        if s == 1
            kii{s} = k(2:n-1, 2:n-1);
            kib{s} = k(2:n-1, n);
            kbb{s} = k(n, n);
            A{s} = zeros(N-1,1);
            A{s}(1) = 1;
            Ad{s} = zeros(N-1,1);
            Ad{s}(1) = 1;
            bp{s} = zeros(1,1);
            num_local_interface{s} = 1;
        elseif s == N
            kii{s} = k(2:n, 2:n);
            kib{s} = k(2:n, 1);
            kbb{s} = k(1, 1);
            A{s} = zeros(N-1,1);
            A{s}(N-1) = 1;
            Ad{s} = zeros(N-1,1);
            Ad{s}(N-1) = -1;
            bp{s} = Fd*ones(1,1);
            num_local_interface{s} = 2*(N-1);
        else
            kii{s} = k(2:n-1, 2:n-1);
            kib{s} = k(2:n-1, [1, n]);
            kbb{s} = k([1, n], [1, n]);
            A{s} = zeros(N-1,2);
            A{s}(s-1,1) = 1;
            A{s}(s,2) = 1;
            Ad{s} = zeros(N-1,2);
            Ad{s}(s-1,1) = -1;
            Ad{s}(s,2) = 1;
            bp{s} = zeros(2,1);
            num_local_interface{s} = [max(num_local_interface{s-1})+1; max(num_local_interface{s-1})+2];
        end
        Sp{s} = kbb{s} - kib{s}' * (kii{s} \ kib{s});
        Sd{s} = pinv(Sp{s});
        Rb{s} = null(Sp{s});
    end
    concat_Sp = blkdiag(Sp{:});
    concat_Sd = blkdiag(Sd{:});
    concat_A = horzcat(A{:});
    concat_Ad = horzcat(Ad{:});
    concat_bp = zeros(size(concat_A,2),1);
    concat_bp(end) = Fd;

    Sp_assembled = concat_A * concat_Sp * concat_A';

    % ---------------------------------------------------------------
    % Construction du préconditionneur de Neumann
    Sp_tilde = inv(concat_Ad*concat_Ad')*(concat_Ad*concat_Sd*concat_Ad')*inv(concat_Ad*concat_Ad')';

    % ---------------------------------------------------------------
    % Initialisation du gradient conjugué préconditionné

    u = zeros(N*(n-1) + 1, 1);
    f = zeros(N*(n-1) + 1, 1);
    f(N*(n-1) + 1) = Fd;
    concat_ub = concat_A' * u(interface);

    rb = cell(N,1);
    for s=1:N
        idx_b = intersect(interface, substruct{s});
        idx_i = setdiff(substruct{s}, [1, interface]); 
        u(idx_i) = kii{s} \ (f(idx_i) - kib{s} * u(idx_b));
        rb{s} = bp{s} - Sp{s} * u(idx_b);
    end

    concat_rb = vertcat(rb{:});
    rb_assembled = concat_A * concat_rb;

    % Initialisation
    r = rb_assembled;
    db = Sp_tilde*r;

    iter=0;
    tol = 1e-3;
    iter_max = 1000;  % Sécurité contre boucle infinie
    
    while(norm(r)>tol && iter < iter_max)
        r_old = r;
        alpha = r'*Sp_tilde*r/(db'*Sp_assembled*db);
        u(interface) = u(interface) + alpha*db;
        r = r - alpha*Sp_assembled*db;
        beta = r'*Sp_tilde*r/(r_old'*Sp_tilde*r_old);
        db = Sp_tilde*r + beta*db;
        iter = iter + 1;
    end

    % Stockage du nombre d'itérations
    iterations(idx) = iter;
end

% ---------------------------------------------------------------
% Tracé du graphe nombre d'itérations en fonction de N

figure;
plot(N_values, iterations, 'LineWidth', 3);
grid on;
xlabel('Nombre de sous-domaines N');
ylabel('Nombre d''itérations');
title('Convergence du gradient conjugué préconditionné en fonction de N');