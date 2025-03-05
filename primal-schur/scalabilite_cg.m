% Scalabilité de la méthode du gradient conjugué pour N=5, N=8 et N=10
% Pour plus de commentaires sur la méthode, voir le fichier
% conjugate_gradient.m

% Données du problème
L = 100;
S = 10;
E = 2*1e5;
Fd = 10;

% Nombres d'élements
n_values = [10, 15, 20];
tab_N = ones(1,30);
tab_iter = ones(3,30);

figure;
hold on;
grid on;
title('Nombre d''itérations en fonction de N');
xlabel('N');
ylabel('Nombre d''itérations');
xlim([1, 30]);

for N = 2:30
    
    for idx = 1:length(n_values)

        % ---------------------------------------------------------------
        % Étape 1 - Remplissage des sous-domaines
        % ---------------------------------------------------------------
        n = n_values(idx);
        H = L / N;
        h = H / (n - 1);
        k0 = E * S / h;

        interface = (n-1)*(1:N-1) + 1;
        k = k0 * (2*eye(n) - diag(ones(n-1, 1), 1) - diag(ones(n-1, 1), -1));
        k(1,1) = k0;
        k(n,n) = k0;
        
        substruct = cell(N,1);
        kbb = cell(N,1);
        kib = cell(N,1);
        kii = cell(N,1);
        Sp_local = cell(N,1);
        Rb_local = cell(N,1);
        A_local = cell(N,1);
        bp_local = cell(N,1);
        num_local_interface = cell(N,1);
        
        for s = 1:N
            substruct{s} = (s-1)*(n-1) + (1:n);
            if s == 1
                kii{s} = k(2:n-1, 2:n-1);
                kib{s} = k(2:n-1, n);
                kbb{s} = k(n, n);
                A_local{s} = zeros(N-1,1);
                A_local{s}(1) = 1;
                bp_local{s} = zeros(1,1);
                num_local_interface{s} = 1;
            elseif s == N
                kii{s} = k(2:n, 2:n);
                kib{s} = k(2:n, 1);
                kbb{s} = k(1, 1);
                A_local{s} = zeros(N-1,1);
                A_local{s}(N-1) = 1;
                bp_local{s} = Fd*ones(1,1);
                num_local_interface{s} = 2*(N-1);
            else
                kii{s} = k(2:n-1, 2:n-1);
                kib{s} = k(2:n-1, [1, n]);
                kbb{s} = k([1, n], [1, n]);
                A_local{s} = zeros(N-1,2);
                A_local{s}(s-1,1) = 1;
                A_local{s}(s,2) = 1;
                bp_local{s} = zeros(2,1);
                num_local_interface{s} = [max(num_local_interface{s-1})+1; max(num_local_interface{s-1})+2];
            end
            Sp_local{s} = kbb{s} - kib{s}' * (kii{s} \ kib{s});
            Rb_local{s} = null(Sp_local{s});
        end
        assembled_Sp = blkdiag(Sp_local{:});
        assembled_A = horzcat(A_local{:});
        assembled_bp = zeros(size(assembled_A,2),1);
        assembled_bp(end) = Fd;
        
        Sp = assembled_A * assembled_Sp * assembled_A';

        % ---------------------------------------------------------------
        % Étape 2 - Algorithme du gradient conjugué
        % ---------------------------------------------------------------

        u = zeros(N*(n-1) + 1, 1);
        f = zeros(N*(n-1) + 1, 1);
        f(N*(n-1) + 1) = Fd;
        assembled_ub = assembled_A'*u(interface);
        rb_local = cell(N,1);
        for s=1:N  
            idx_b = intersect(interface, substruct{s});
            idx_i = setdiff(substruct{s}, [1, interface]); 
            u(idx_i) = kii{s} \ (f(idx_i) - kib{s} * u(idx_b));
            rb_local{s} = bp_local{s} - Sp_local{s}*u(idx_b);
        end
  
        assembled_rb = vertcat(rb_local{:});
        rb = assembled_A*assembled_rb;
        db = rb;
        
        iter=0;
        while(norm(rb)>1e-3)
            assembled_db = assembled_A'*db;
            for s=1:N
                idx_i = setdiff(substruct{s}, [1, interface]); 
                di = kii{s} \ (f(idx_i) - kib{s} * assembled_db(num_local_interface{s}));
            end
            alpha = rb'*rb/(db'*Sp*db);
            u(interface) = u(interface) + alpha*db;
            rb = rb - alpha*Sp*db;
            beta = rb'*Sp*db/(db'*Sp*db);
            db = rb + beta*db;
            iter = iter + 1;
        end

        % ---------------------------------------------------------------
        % Étape 3 - Stockage du nombre d'itérations
        % ---------------------------------------------------------------
        tab_iter(idx, N) = iter;

    end % fin de la boucle sur n
    tab_N(N) = N;
end % fin de la boucle sur N

% ---------------------------------------------------------------
% Étape 4 - Graphique
% ---------------------------------------------------------------
for idx = 1:length(n_values)
    plot(tab_N, tab_iter(idx, :), 'LineWidth', 3, 'DisplayName', sprintf('n=%d', n_values(idx)));
end

legend;