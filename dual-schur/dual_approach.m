% Résolution du problème par la méthode duale

% Données du problème
L = 100;
S = 10;
E = 2*1e5;
Fd = 10;

% Choix du nombre de domaines
N = input("Nombre de domaines N (> 1): ");
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

% Initialisation
substruct = cell(N,1);
kbb = cell(N,1);
kib = cell(N,1);
kii = cell(N,1);
Sp = cell(N,1);
Sd = cell(N,1);
Rb = cell(N,1);
A = cell(N,1);
num_local_interface = cell(N,1);

for s = 1:N
    % Définition des sous-structures
    substruct{s} = (s-1)*(n-1) + (1:n);
    
    % Découpage de la matrice de rigidité locale
    if s == 1
        kii{s} = k(2:n-1, 2:n-1);
        kib{s} = k(2:n-1, n);
        kbb{s} = k(n, n);
        num_local_interface{s} = 1;
    elseif s == N
        kii{s} = k(2:n, 2:n);
        kib{s} = k(2:n, 1);
        kbb{s} = k(1, 1);
        num_local_interface{s} = 2*(N-1);
    else
        kii{s} = k(2:n-1, 2:n-1);
        kib{s} = k(2:n-1, [1, n]);
        kbb{s} = k([1, n], [1, n]);
        num_local_interface{s} = [max(num_local_interface{s-1})+1; max(num_local_interface{s-1})+2];
    end
    
    % Primal Schur complement
    Sp{s} = kbb{s} - kib{s}' * (kii{s} \ kib{s});
    if s==N
        Sp{s} = 0; % entraine une erreur de précision si absent
    end
    
    % Dual Schur complement
    Sd{s} = pinv(Sp{s});
    
    % Rigid body modes
    Rb{s} = null(Sp{s}, 1e-6);
    
    % Primal assembly operator
    if s == 1
        A{s} = zeros(N-1,1);
        A{s}(1) = 1;
    elseif s > 1 && s < N
        A{s} = zeros(N-1,2);
        A{s}(s-1,1) = -1;
        A{s}(s,2) = 1;
    elseif s == N
        A{s} = zeros(N-1,1);
        A{s}(N-1) = -1;
    end
end

% Concaténation des opérateurs
concat_Sd = blkdiag(Sd{:});
concat_A = horzcat(A{:});
concat_bd = zeros(size(concat_A,2),1);
concat_bp = zeros(size(concat_A,2),1);
concat_bp(end) = Fd;
concat_Rb = blkdiag(Rb{:});

% Assemblage des opérateurs
Sd_assembled = concat_A * concat_Sd * concat_A';
bd_assembled = concat_A * concat_bd;
G = concat_A * concat_Rb;
e = concat_Rb' * concat_bp;

% Résolution du problème
lambda_b = G' \ (-e);
alpha_b = G \ (-Sd_assembled*lambda_b - bd_assembled);
concat_lambda_b = concat_A' * lambda_b;

% Calcul du déplacement aux interfaces
for s = 1:N
    idx_b = intersect(interface, substruct{s});
    b = num_local_interface{s};
    if(s==1)
        u(idx_b) = Sd{s} * (concat_bp(b) + concat_lambda_b(b));
    else
        u(idx_b) = Sd{s} * (concat_bp(b) + concat_lambda_b(b)) + Rb{s} * alpha_b(s-1);
    end
    
end

% Calcul des déplacements internes
for s = 1:N
    idx_b = intersect(interface, substruct{s});
    idx_i = setdiff(substruct{s}, [1, interface]);
    u(idx_i) = kii{s} \ (f(idx_i) - kib{s} * u(idx_b));
end

% Affichage des solutions
disp('Matrice de Schur :');
disp(Sd_assembled);
disp('bd :');
disp(bd_assembled);
disp('Déplacements :');
disp(u);

