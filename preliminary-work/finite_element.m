% Données du problème
L = 100;       
S = 10;
E = 2*1e5;
Fd = 10;

% Choix du nombre d'élements
n = input("Nombre d'élements n : ");
h = L/(n-1);

% Matrice de rigidité (on enlève la première ligne et la première colonne)
k = E*S/h * 2*eye(n-1,n-1) - E*S/h*diag(ones(n-2, 1), 1) - E*S/h*diag(ones(n-2, 1), -1);
k(n-1,n-1) = E*S/h;       

% Matrice des forces
F = zeros(n-1,1);
F(n-1) = Fd;

% Résolution du problème
u = k \ F;

% Affichage solution
disp("k = ")
disp(k)
disp("F = ")
disp(F)
disp("u = ")
disp(u)