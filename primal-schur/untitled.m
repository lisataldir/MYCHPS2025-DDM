% Conditionnement de la matrice de Schur pour N=5, N=8 et N=10

% Données du problème
L = 100;
S = 10;
E = 2*1e5;

% Choix du nombre d'éléments
n_values = [10, 15, 20];
tab_N = ones(1,30);
tab_cond = ones(3,30);


figure;
hold on;
grid on;
title('Conditionnement de la matrice de Schur Sp');
xlabel('N');
ylabel('Conditionnement de Sp');
xlim([1, 20]);

for N = 2:30
    
    for idx = 1:length(n_values)
        
        n = n_values(idx);
        H = L / N;
        h = H / (n - 1);
        k0 = E * S / h;

        
        % Initialisation
        kbb = cell(N,1);
        kib = cell(N,1);
        kii = cell(N,1);
        Sp = cell(N,1);
        A = cell(N,1);
        
        for s = 1:N
            % Définition des matrices de rigidité locales
            k = k0 * (2 * eye(n) - diag(ones(n-1, 1), 1) - diag(ones(n-1, 1), -1));
            k(1,1) = k0;
            k(n,n) = k0;
            
            % Découpage de la matrice de rigidité locale
            if s == 1
                kii{s} = k(2:n-1, 2:n-1);
                kib{s} = k(2:n-1, n);
                kbb{s} = k(n, n);
            elseif s == N
                kii{s} = k(2:n, 2:n);
                kib{s} = k(2:n, 1);
                kbb{s} = k(1, 1);
            else
                kii{s} = k(2:n-1, 2:n-1);
                kib{s} = k(2:n-1, [1, n]);
                kbb{s} = k([1, n], [1, n]);
            end
            
            % Primal Schur complement
            Sp{s} = kbb{s} - kib{s}' * (kii{s} \ kib{s});
            
            % Primal assembly operator
            if s == 1
                A{s} = zeros(N-1,1);
                A{s}(1) = 1;
            elseif s > 1 && s < N
                A{s} = zeros(N-1,2);
                A{s}(s-1,1) = 1;
                A{s}(s,2) = 1;
            elseif s == N
                A{s} = zeros(N-1,1);
                A{s}(N-1) = 1;
            end
        end
        
        % Assemblage global
        assembled_Sp = blkdiag(Sp{:});
        assembled_A = horzcat(A{:});
        
        % Résolution du problème
        Sp_global = assembled_A * assembled_Sp * assembled_A';
        
        tab_cond(idx, N) = cond(Sp_global);

    end % fin de la boucle sur n
    tab_N(N) = N;
end % fin de la boucle sur N
for idx = 1:length(n_values)
    plot(tab_N, tab_iter(idx, :), 'LineWidth', 3, 'DisplayName', sprintf('n=%d', n_values(idx)));
end
legend;