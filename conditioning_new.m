%% data
clc
clear all

E=200e9;        % Young Modulules Pa
S=10e-6;        % surface area 10e-6m²=10 mm²;
ES=E*S;
L=0.100;        % 100 mm length bar
F_d=10;
%ne       % element number per substructure
%Ns       % Number of substructure

%% Sous-structures

conditionnementS = zeros(20,30);
conditionnementK = zeros(20,30);
taille_h = zeros(20,30);
valeurs_h_sur_H = zeros(20,30); 

for Ns=1:10 % variable:nb de sous structures 
    for Ne=1:20 % variable:nb d'éléments par sous structures 
        n = Ne*Ns; %nb d'elements
        h = L/n; %taille des élements 
        H = L/Ns; %taille sous-structures
        coef = E*S/h; %pour les multiplicateurs de Lagrange
        taille_h(Ns,Ne)=h;
        valeurs_h_sur_H(Ns,Ne)=h/H;
        
        % rigidité globale 
        Ke = (E*S/h)*[1 -1; -1 1]; % contribution elementaire
        K = zeros(n+1);
        for i=1:n
            K(i,i)=K(i,i)+Ke(1,1);
            K(i+1,i)=K(i+1,i)+Ke(2,1);
            K(i,i+1)=K(i,i+1)+Ke(1,2);
            K(i+1,i+1)=K(i+1,i+1)+Ke(2,2);
        end
        %conditions aux limites deplacement
        K(1,:)=0;
        K(1,1)=(E*S/h);
        
        conditionnementK(Ns,Ne)=cond(K);

        
        % on decrit les interfaces 

        interface = zeros(Ns-1,2); % va etre la matrice contenant les couples des noeuds d'interface
        for i = 1:Ns-1
        interface(i,1) = (Ne+1)*i;
        interface(i,2) = (Ne+1)*i+1;   
        end

        % On definit les operateurs 
        % Les matrices  des sous structures sont toutes identiques
        K_sss = zeros(Ne+1);
        for i=1:Ne
            K_sss(i,i)=K_sss(i,i)+Ke(1,1);
            K_sss(i+1,i)=K_sss(i+1,i)+Ke(2,1);
            K_sss(i,i+1)=K_sss(i,i+1)+Ke(1,2);
            K_sss(i+1,i+1)=K_sss(i+1,i+1)+Ke(2,2);
        end

        % on reorganise  les lignes et colonnes afin d'avoir la notation par bloc
        % interne/externe

        % a est une variable de stockage
        % On change les lignes
        a = K_sss(1,:);
        K_sss(1:Ne-1,:) = K_sss(2:Ne,:);
        K_sss(Ne,:) = a;
        % Puis les colonnes
        a = K_sss(:,1);
        K_sss(:,1:Ne-1) = K_sss(:,2:Ne);
        K_sss(:,Ne) = a;

        %Calcul du complement de Schur primal de sous-structures.

        K_bb = K_sss(Ne:Ne+1,Ne:Ne+1);
        K_ii = K_sss(1:Ne-1,1:Ne-1);
        K_ib = K_sss(1:Ne-1,Ne:Ne+1);
        K_bi = K_sss(Ne:Ne+1,1:Ne-1);
        S_p_s = K_bb-K_bi*inv(K_ii)*K_ib;

        % Calcul des opérateurs assemblés
        S_p=zeros(Ns+1);
        for i = 1:Ns
            S_p(i,i)=S_p(i,i)+S_p_s(1,1);
            S_p(i+1,i)=S_p(i+1,i)+S_p_s(2,1);
            S_p(i,i+1)=S_p(i,i+1)+S_p_s(1,2);
            S_p(i+1,i+1)=S_p(i+1,i+1)+S_p_s(2,2);
        end
        %conditions aux limites deplacement
        S_p(1,:)=0;
        S_p(1,1)=(E*S/h);

        conditionnementS(Ns,Ne) = cond(S_p);
    end
end


%conditionnement? 
figure
X = reshape(taille_h,[numel(taille_h),1]);
Y = reshape(valeurs_h_sur_H,[numel(valeurs_h_sur_H),1]);
Z = reshape(conditionnementK,[numel(conditionnementK),1]);
scatter3(X,Y,Z,60,Z,'filled');
xlabel('h')
ylabel('h/H')
zlabel('Conditionnement de K en norme 2') 
%
figure
X = reshape(taille_h,[numel(taille_h),1]);
Y = reshape(valeurs_h_sur_H,[numel(valeurs_h_sur_H),1]);
Z = reshape(conditionnementS,[numel(conditionnementS),1]);
scatter3(X,Y,Z,60,Z,'filled');
xlabel('h')
ylabel('h/H')
zlabel(' Sp conditioning ') 

% xlabel('')

b_p = zeros(Ns+1,1);
b_p(Ns+1)= F_d;
% % résolution inconnues de bord
% u_p = linsolve(S_p, b_p);
