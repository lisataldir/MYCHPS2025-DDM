% Conditionnement de la matrice de rigidité

% Données du problème
L = 100;
S = 10;
E = 2*1e5;

N = 1; % Problème non-structuré

tab_n = ones(1,50);
tab_cond = ones(1,50);

for n=2:50
    % Constantes
    h = N / (n - 1);
    k0 = E * S / h;
    
    % Définition matrice de rigidité globale
    k = k0 * (2*eye(n) - diag(ones(n-1, 1), 1) - diag(ones(n-1, 1), -1));
    k(1,1) = k0;
    k(n,n) = k0;

    % Calcul et stockage du conditionnement de K
    % On retire la 1ère ligne et la 1ère colonne pour ne pas avoir une
    % matrice singulière. On peut le faire car u1 = 0.
    tab_n(n) = n;
    tab_cond(n) = cond(k(2:n,2:n));
end
      
% Plot
figure;
subplot(2,1,1);
plot(tab_n, tab_cond, 'LineWidth', 4);
xlabel('n');
ylabel('Conditionnement de K');
ylim([0, 4500]);
xlim([1, 50]);
title('Conditionnement de la matrice de rigidité K');
grid on;
