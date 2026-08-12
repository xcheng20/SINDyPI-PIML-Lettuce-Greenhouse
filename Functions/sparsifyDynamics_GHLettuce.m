function [Xi_lambda,ODEs] = sparsifyDynamicsDW(Theta,dXdt,LHS_Sym,lambda,N,Sym_Struct,disp,NormalizeLib)
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

% compute Sparse regression: sequential least squares
%%
% Modified By: K
% Last Updated
%% Normalize the library data
if NormalizeLib==1
    % Change thisinto parfor to increase the speed if you have a large
    % library
   % Norm itu bukan bikin 0 tapi ngecut numbernya aja!

     parfor norm_k=1:size(Theta,2)
        normLib(norm_k) = norm(Theta(:,norm_k));
        Theta(:,norm_k) = Theta(:,norm_k)/normLib(norm_k);  %untuk RHS_Data ini pasti sama semua sih jadinya
     end
else
    for norm_k=1:size(Theta,2)
        normLib(norm_k) = norm(Theta(:,norm_k));
        Theta(:,norm_k) = Theta(:,norm_k)/normLib(norm_k);  %untuk RHS_Data ini pasti sama semua sih jadinya
     end
end

%% Peform sparse regression
Xi = Theta\dXdt;  % initial guess: Least-squares 
[n,m]=size(dXdt);

% lambda is our sparsification knob.  %How to determine lambda value? ->
% parameter coeff
for k=1:N
    smallinds = (abs(Xi)<lambda);   % find small coefficients
    Xi(smallinds)=0;                % and threshold
    for ind = 1:m                   % n is state dimension
        biginds = ~smallinds(:,ind);
        % Regress dynamics onto remaining terms to find sparse Xi
        Xi(biginds,ind) = Theta(:,biginds)\dXdt(:,ind);
    end
end

Xi_lambda = Xi;
%% Now output the SINDy Identified ODEs

% Now retrive the parameters
if NormalizeLib==1
    for norm_k=1:length(Xi)
        Xi(norm_k,:) = Xi(norm_k,:)/normLib(norm_k);   %How about the Xi does not need to be normalized?
    end
end

for i=1:m
     ODEs(i,1)=vpa(cell2sym(Sym_Struct)*Xi(:,i));  
end

%Numeric Array
% Num_Xi = nnz(Xi);

%% Choose whether you want to display the final discovered equation

if disp==1
     fprintf('The SINDy-PI discovered ODE is:\n')
     digits(6)
     for i=1:m
          fprintf(strcat('\t',char(LHS_Sym),'=',char(simplify(ODEs(i,1))),'\n'));  %this simplify does not make a different
     end
end
