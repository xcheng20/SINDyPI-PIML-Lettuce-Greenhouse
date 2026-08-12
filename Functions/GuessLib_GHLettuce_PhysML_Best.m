% This Guess library (Left Hand Side) is coded for lettuce dry weight example
% Last Updated: 10/04/2023
% Coded By: Fakhira

function [Data,Sym_Struct]=GuessLib_GHTomato_PhysML_Best(X,dX,u,d,iter,Highest_dXPoly_Order_Guess)

%% First get the size of the X matrix, determin the data length and the number of variables we have.
[Data_Length,Variable_Number_x]=size(X);
[~,Variable_Number_dX]=size(dX);
[~,Variable_Number_u]=size(u);
[~,Variable_Number_d]=size(d);
%Also create the symbolic variable
Symbol_x=sym('x',[Variable_Number_x,1]);
Symbol_dX=sym('dx',[Variable_Number_dX,1]);
Symbol_u=sym('u',[Variable_Number_u,1]);
Symbol_d=sym('d',[Variable_Number_d,1]);

%Now according the Highest Polynomial Order entered, we will calculate the data matrix.
Data=[];
Index=1; 

%POLY 2
% State 1
if Highest_dXPoly_Order_Guess == 1
    %Physical Knowledge for X1
    %d1^2dx1
    Basis_d     = d(:,1);
    Basis_d_Sym = Symbol_d(1);
        Data(:,Index)=dX(:,iter).*Basis_d(:,1).^2;
        Sym_Struct{1,Index}=Symbol_dX(iter,1)*Basis_d_Sym(1,1).^2;
        Index=Index+1;
end 

% State 2
if Highest_dXPoly_Order_Guess == 2
    %dx2x3
    Data(:,Index)=dX(:,iter).*X(:,3); %THE BEST
    Sym_Struct{1,Index}=Symbol_dX(iter,1)*Symbol_x(3,1);  
end

% State 3
if Highest_dXPoly_Order_Guess == 3
    %Without dX
    Data(:,Index)=dX(:,3);
    Sym_Struct{1,Index}=Symbol_dX(3,1);
    Index=Index+1;
end

%State 4
if Highest_dXPoly_Order_Guess == 4
      %dx4d4
      %Physical Knowledge for X4
      Basis_d     = d(:,4);
      Basis_d_Sym = Symbol_d(4);
        Data(:,Index)=dX(:,iter).*Basis_d(:,1);
        Sym_Struct{1,Index}=Symbol_dX(iter,1)*Basis_d_Sym(1,1);
        Index=Index+1;
end 
  