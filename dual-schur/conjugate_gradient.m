% Algorithme du gradient conjugué pour la méthode duale

% Données du problème
L = 100;
S = 10;
E = 2*1e5;
Fd = 10;

% Choix du nombre de domaines et du nombre d'élements par domaine
N = input("Nombre de domaines N (> 1): ");
H = L/N;
n = input("Nombre d'éléments n (> 1): ");
h = H/(n-1);


% ---------------------------------------------------------------

% Calcul préliminaire nécessaire pour le gradient conjugué (voir
% dual_approach.m pour avoir un code commenté)

u = zeros(N*(n-1) + 1, 1);
f = zeros(N*(n-1) + 1, 1);
f(N*(n-1) + 1) = Fd;
k0 = E*S/h;
interface = (n-1)*(1:N-1) + 1;

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
num_local_interface = cell(N,1);

for s = 1:N
    substruct{s} = (s-1)*(n-1) + (1:n);
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
    Sp{s} = kbb{s} - kib{s}' * (kii{s} \ kib{s});
    if s==N
        Sp{s} = 0; 
    end

    Sd{s} = pinv(Sp{s});
    Rb{s} = null(Sp{s});

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

concat_Sd = blkdiag(Sd{:});
concat_A = horzcat(A{:});
concat_bd = zeros(size(concat_A,2),1);
concat_bp = zeros(size(concat_A,2),1);
concat_bp(end) = Fd;
concat_Rb = blkdiag(Rb{:});
Sd_assembled = concat_A * concat_Sd * concat_A';
bd_assembled = concat_A * concat_bd;
G = concat_A * concat_Rb;
e = concat_Rb' * concat_bp;

% Fin du calcul préliminaire

% ---------------------------------------------------------------


% Début de l'algorithme du gradient conjugué 


% Initialisation
P = eye(length(G)) - G*inv(G'*G)*G';
lambda = - G*inv(G'*G)*e;
r = P'*(-bd_assembled - Sd_assembled*lambda);
z = P*inv(Sd_assembled)*r;
d = z;

% Boucle
d_values = [];
iter = 1;
while (norm(r) > 1e-5)
    p = P'*Sd_assembled*d;
    alpha = r'*d/(d'*p);
    lambda = lambda + alpha*d;
    r = r - alpha*p;
    z = P*inv(Sd_assembled)*r;
    beta = zeros(iter+1);
    d_values = [d_values, d]; 
    for j=0:iter
        beta(j) = -z/d_values(j);
    end
    d = z + sum(beta)*d;
    % Itération
    iter = iter + 1;
end

% Post-processing
concat_alpha_b = inv(G'*G) * G' * (-bd_assembled - Sd_assembled*lambda);
concat_ub = concat_Sd*(concat_bp + concat_A'*lambda) + concat_Rb*concat_alpha_b;

disp(concat_ub);
disp(iter);