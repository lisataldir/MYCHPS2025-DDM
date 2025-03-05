% Données du problème
L = 100;
S = 10;
E = 2*1e5;
Fd = 10;

% Choix du nombre de domaines et du nombre d'élements par domaine
N = input("Nombre de domaines N : ");
H = L/N;
n = input("Nombre d'éléments n (> 1): ");
h = H/(n-1);


% ---------------------------------------------------------------


% Calcul préliminaire nécessaire pour le gradient conjugué (voir
% primal_approach.m pour avoir un code commenté)
interface = (n-1)*(1:N-1) + 1;

k0 = E*S/h;
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
% Fin du calcul préliminaire


% ---------------------------------------------------------------


% Début de l'algorithme du gradient conjugué 

% Initialisation 
u = zeros(N*(n-1) + 1, 1);
f = zeros(N*(n-1) + 1, 1);
f(N*(n-1) + 1) = Fd;
assembled_ub = assembled_A'*u(interface);
rb_local = cell(N,1);
for s=1:N
    % Résolution du problème de Dirichlet
    idx_b = intersect(interface, substruct{s});
    idx_i = setdiff(substruct{s}, [1, interface]); 
    u(idx_i) = kii{s} \ (f(idx_i) - kib{s} * u(idx_b));
    % Calcul du résidu rb local
    rb_local{s} = bp_local{s} - Sp_local{s}*u(idx_b);
end
% Calcul du résidu global
assembled_rb = vertcat(rb_local{:});
rb = assembled_A*assembled_rb;
db = rb;

iter=0;
while(norm(rb)>1e-3)
    assembled_db = assembled_A'*db;
    for s=1:N
        % Résolution du problème de Dirichlet
        idx_i = setdiff(substruct{s}, [1, interface]); 
        di = kii{s} \ (f(idx_i) - kib{s} * assembled_db(num_local_interface{s}));
    end
    % Mise-à-jour d'alpha, ub, rb, beta et db
    alpha = rb'*rb/(db'*Sp*db);
    u(interface) = u(interface) + alpha*db;
    rb = rb - alpha*Sp*db;
    beta = rb'*Sp*db/(db'*Sp*db);
    db = rb + beta*db;
    % Itération
    iter = iter + 1;
end

% Affichage des solutions
disp('Déplacements aux frontieres:');
disp(u(interface));

disp('Nombre d''itérations : ');
disp(iter);