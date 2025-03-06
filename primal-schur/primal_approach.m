% Résolution par la méthode primale

% Données du problème
L = 100;
S = 10;
E = 2*1e5;
Fd = 10;

% Choix du nombre de domaines
N = input("Nombre de domaines N : ");
H = L/N;

% Choix du nombre d'éléments par domaine
n = input("Nombre d'éléments n (> 1): ");
h = H/(n-1);

% Initialisation du vecteur des déplacements et des forces
u = zeros(N*(n-1) + 1, 1);
f = zeros(N*(n-1) + 1, 1);
f(N*(n-1) + 1) = Fd;

% Constante de rigidité
k0 = E*S/h;

% Définition des interfaces 
interface = (n-1)*(1:N-1) + 1;

% Définition des matrices de rigidités locales
k = k0 * (2*eye(n) - diag(ones(n-1, 1), 1) - diag(ones(n-1, 1), -1));
k(1,1) = k0;
k(n,n) = k0;

% Initialisation des sous-structures
substruct = cell(N,1);
kbb = cell(N,1);
kib = cell(N,1);
kii = cell(N,1);
Sp = cell(N,1);
Rb = cell(N,1);
A = cell(N,1);


for s = 1:N
    % Définition des sous-structures
    substruct{s} = (s-1)*(n-1) + (1:n);
    
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
    
    % Rigid body modes
    Rb{s} = null(Sp{s});
    
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

% Concatenation
concat_Sp = blkdiag(Sp{:});
concat_A = horzcat(A{:});
concat_bp = zeros(size(concat_A,2),1);
concat_bp(end) = Fd;

% Assemblage
Sp_assembled = concat_A * concat_Sp * concat_A';
bp_assembled = concat_A * concat_bp;

% Résolution du problème

% Calcul du déplacement aux interfaces
u(interface) = Sp_assembled \ bp_assembled;

% Calcul des déplacements internes
for s = 1:N
    idx_b = intersect(interface, substruct{s});
    idx_i = setdiff(substruct{s}, [1, interface]); % On retire aussi le point 1 puisqu'il vaut 0
    u(idx_i) = kii{s} \ (f(idx_i) - kib{s} * u(idx_b));
end

% Affichage des solutions
disp('Matrice primale de Schur :');
disp(Sp_assembled);
disp('bp :');
disp(bp_assembled);
disp('Déplacements :');
disp(u);
