clear; clc;

%% Given integer (rational) matrix

% A numerical example

F = diag([2, 2, -1, -1]);


% or some random examples..
% F = randi(20, 10)-10;

% A1 = diag([4,4,4,4,4,4,4,4]) + diag([1,1,0,1,1,0,0],1);
% A2 = diag([3,3,3,3,3,3,3,3]) + diag([1,1,0,1,1,0,0],1);
% F = blkdiag(A1, A2);

%% Step 0 : factorize the char poly & min poly


syms x
minF = minpoly(F, x);                     % the minimal polynomial
[fac, eta_lst] = factorWithMult(minF, x); % factors and multiplicities of min poly
l = size(eta_lst, 2);                     % number of irreducible factors
n = size(F,1);                            % dimension of entire vector space

p_lst = {}; % list of irreducible factors
for i = 1:l
    p_lst{i} = sym2poly(fac(i));
end

v = cell(l, 1);                 % cell-of-cells: v{i}{j} = basis vector
delta_lst = cell(l, 1);         % cell-of-cells: delta{i}{j} = degree index
kappa_lst = zeros(1, l);        % number of blocks per primary divisor

%% Step 1 : Find F-cyclic bases
for i = 1:l
    pi = p_lst{i};               % the divisor
    di = length(pi) - 1;         % degree of pi
    Bi = [];                     % basis of ker(pi(F)^\eta_i);
    k = 0;                       % number of generators
    v{i} = {};                   % generators
    delta_lst{i} = {};           % length of F-cyclic basis
    
    for j = eta_lst(i):-1:1
        pF = polyvalm(pi, F);
        null_j = null(pF^j);        % basis for the null space ker(p(F)^j)
        null_jm1 = null(pF^(j-1));  % basis for the null space ker(p(F)^j-1)
        
        rj = (size(null_j,2) - size(null_jm1,2)) / di;
        
        % Find rj linearly independent vectors from ker(p(F)^j) \ ker(p(F)^j-1)
        W = null_j;

        trial = 1;
        while k < rj && trial <= size(W,2)
            w = W(:,trial);
            testSet = [Bi, pF^(j-1) * w];
            if rank(testSet) > rank(Bi) % if w is accepted
                for m = 0:(j*di-1)
                    Bi = [Bi, F^m * w];
                end
                k = k + 1;
                v{i}{k} = w;
                delta_lst{i}{k} = j;
            end
            trial = trial + 1;
        end
    end
    kappa_lst(i) = k;
end

%% Step 2 : Construct rational canonical basis
kappa_max = max(kappa_lst);
V = [];
delta = [];
for j = 1:kappa_max
    Ij = find(cellfun(@(ci) length(ci) >= j, v));  % return index where v{i}{j} exists

    vj = zeros(n,1);
    deltaj = 0;
    for ii = 1:length(Ij)
        idx = Ij(ii);
        d_i = length(p_lst{idx}) - 1;

        vj = vj + v{idx}{j};                        % generator for F-cyclic basis
        deltaj = deltaj + d_i * delta_lst{idx}{j};  % length of F-cyclic basis
    end

    Vj = [vj];
    for m = 1:deltaj-1
        Vj = [F * Vj, vj];
    end
    V = [Vj, V];
    delta = [deltaj, delta];
end
T = inv(V);

%% Coordinate transform

R = round(T * F * V); % desired rational canonical form


% Just copy and paste the result into Go code!
TrGo(R, kappa_max, delta) 

%% useful function
function [factors, multi] = factorWithMult(f, x)
    factors = []; 
    multi = [];
    fact = factor(f,  x, 'FactorMode', 'rational'); % factorize over rational numbers
    while ~isempty(fact)
        factors = [factors, fact(1)];
        indi = find((fact - fact(1)) == 0);
        multi = [multi, size(indi, 2)];
        fact = fact((find((fact - fact(1)) ~= 0)));
    end
end

function TrGo(F, kappa, delta)
    fprintf('F = [][]float64{\n');
    for i = 1:size(F,1)
        fprintf('\t{');
        fprintf('%.1i, ', F(i,1:end-1));
        fprintf('%.1i},\n', F(i,end));
    end
    fprintf('}\n\n');
    
    fprintf('kappa = %.1i', kappa);
    fprintf('\n\n');
    
    fprintf('r = []int{\n');
    for i = 0:kappa-1
        fprintf('%.1i, \n', sum(delta(1:i)));
    end
    fprintf('}\n');
end