% Conditionnement de la matrice de Schur pour N=5, N=8 et N=10

% Données du problème
L = 100;
S = 10;
E = 2*1e5;

% Choix des nombres de domaines
N_values = [5, 8, 10];

figure;
hold on;
grid on;
title('Conditionnement de la matrice de Schur Sp');
xlabel('n/N');
ylabel('Conditionnement de Sp');
ylim([0, 200]);
xlim([1, 30]);

for idx = 1:length(N_values)
    N = N_values(idx);
    tab_n = zeros(1,30);
    tab_cond = zeros(1,30);
    
    for i = 1:30
        % Constantes
        n = i * N;  
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
        
        tab_n(i) = i;
        tab_cond(i) = cond(Sp_global);
    end
  
    plot(tab_n, tab_cond, 'LineWidth', 4, 'DisplayName', sprintf('N=%d', N));
end

legend;
