% This Guess library (Left Hand Side) is coded for lettuce dry weight example
% Last Updated: 10/04/2023
% Coded By: Fakhira

%The Guesslibrary is built based on the real model structure, which is:
% 1. Only individual terms of state and derivatives
% 2. Without U (control input)
% 3. Polynomial order 2 consists of only state variables
% Reasoning: The inclusion of XD, and U makes the algorithm fails to generate the
% model >> somehow stuck >> out of memory

function [Data,Sym_Struct]=GuessLib_GHLettuce_PhysML_Test(X,dX,u,d,iter,Highest_dXPoly_Order_Guess)

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

%POLY 1
if iter ==1
%Best GUESS of LHS from the previous test
    if Highest_dXPoly_Order_Guess >= 0
    %Physical Knowledge for X1
    Basis_x     = [X(:,1) X(:,2) X(:,3)]; 
    Basis_x_Sym = [Symbol_x(1) Symbol_x(2) Symbol_x(3)];
    Basis_d     = d(:,1);
    Basis_d_Sym = Symbol_d(1);
        Data(:,Index)=dX(:,iter).*Basis_d(:,1).^2;
        Sym_Struct{1,Index}=Symbol_dX(iter,1)*Basis_d_Sym(1,1).^2;
        Index=Index+1;
    end
    if Highest_dXPoly_Order_Guess >= 1
    %Physical Knowledge for X1
    Basis_x     = [X(:,1) X(:,2) X(:,3)]; 
    Basis_x_Sym = [Symbol_x(1) Symbol_x(2) Symbol_x(3)];
    Basis_d     = d(:,1);
    Basis_d_Sym = Symbol_d(1);
        %Poly 1 X
        for i=1:size(Basis_x,2)
            Data(:,Index)=dX(:,iter).*Basis_x(:,i);
            Sym_Struct{1,Index}=Symbol_dX(iter,1)*Basis_x_Sym(1,i);
            Index=Index+1;
        end
        %Poly 1 D
        for i=1:size(Basis_d,2)
            Data(:,Index)=dX(:,iter).*Basis_d(:,i);
            Sym_Struct{1,Index}=Symbol_dX(iter,1)*Basis_d_Sym(1,i);
            Index=Index+1;
        end
    end
    if Highest_dXPoly_Order_Guess >= 2
        %Poly 2 X
        for i=1:size(Basis_x,2)
            for j=i:size(Basis_x,2)
                Data(:,Index)=dX(:,iter).*Basis_x(:,i).*Basis_x(:,j);
                Sym_Struct{1,Index}=Symbol_dX(iter,1)*Basis_x_Sym(1,i)*Basis_x_Sym(1,j);
                Index=Index+1;
            end
        end
        Basis_x = [X(:,2) X(:,3)];
        Basis_x_Sym = [Symbol_x(2) Symbol_x(3)];
        for i=1:size(Basis_x,2)
            Data(:,Index)=dX(:,iter).*X(:,1).^2.*Basis_x(:,i);
            Sym_Struct{1,Index}=Symbol_dX(iter,1)*Symbol_x(1,1).^2*Basis_x_Sym(1,i);
            Index=Index+1;
        end
        Basis_x = [X(:,1) X(:,3)];
        Basis_x_Sym = [Symbol_x(1) Symbol_x(3)];
        for i=1:size(Basis_x,2)
            Data(:,Index)=dX(:,iter).*X(:,2).^2.*Basis_x(:,i);
            Sym_Struct{1,Index}=Symbol_dX(iter,1)*Symbol_x(2,1).^2*Basis_x_Sym(1,i);
            Index=Index+1;
        end
        Basis_x = [X(:,1) X(:,2)];
        Basis_x_Sym = [Symbol_x(1) Symbol_x(2)];
        for i=1:size(Basis_x,2)
            Data(:,Index)=dX(:,iter).*X(:,3).^2.*Basis_x(:,i);
            Sym_Struct{1,Index}=Symbol_dX(iter,1)*Symbol_x(3,1).^2*Basis_x_Sym(1,i);
            Index=Index+1;
        end
        %Poly 2 D
        for i=1:size(Basis_d,2)
           for j=i:size(Basis_d,2)
               Data(:,Index)=dX(:,iter).*Basis_d(:,i).*Basis_d(:,j);
               Sym_Struct{1,Index}=Symbol_dX(iter,1)*Basis_d_Sym(1,i)*Basis_d_Sym(1,j);
               Index=Index+1;
           end
        end 
    end
end 

if iter ==2 
     if Highest_dXPoly_Order_Guess >= 0
        Data(:,Index)=dX(:,iter).*X(:,3);
        Sym_Struct{1,Index}=Symbol_dX(iter,1)*Symbol_x(3,1);
        Index=Index+1;
     end 
     if Highest_dXPoly_Order_Guess >= 1
        %Physical Knowledge for X2
        Basis_x     = [X(:,1) X(:,2) X(:,3)]; 
        Basis_x_Sym = [Symbol_x(1) Symbol_x(2) Symbol_x(3)];
        Basis_d     = [d(:,1) d(:,2)];
        Basis_d_Sym =[Symbol_d(1) Symbol_d(2)];
        Basis_u     = [u(:,1) u(:,2)];
        Basis_u_Sym =[Symbol_u(1) Symbol_u(2)];
        %Poly 1 X
        for i=1:size(Basis_x,2)
            Data(:,Index)=dX(:,iter).*Basis_x(:,i);
            Sym_Struct{1,Index}=Symbol_dX(iter,1)*Basis_x_Sym(1,i);
            Index=Index+1;
        end
        %Poly 1 D
        for i=1:size(Basis_d,2)
            Data(:,Index)=dX(:,iter).*Basis_d(:,i);
            Sym_Struct{1,Index}=Symbol_dX(iter,1)*Basis_d_Sym(1,i);
            Index=Index+1;
        end
        %Poly 1 U
        for i=1:size(Basis_u,2)
            Data(:,Index)=dX(:,iter).*Basis_u(:,i);
            Sym_Struct{1,Index}=Symbol_dX(iter,1)*Basis_u_Sym(1,i);
            Index=Index+1;
        end
     end
    if Highest_dXPoly_Order_Guess >= 2
      %Poly 2 X
        for i=1:size(Basis_x,2)
            for j=i:size(Basis_x,2)
                Data(:,Index)=dX(:,iter).*Basis_x(:,i).*Basis_x(:,j);
                Sym_Struct{1,Index}=Symbol_dX(iter,1)*Basis_x_Sym(1,i)*Basis_x_Sym(1,j);
                Index=Index+1;
            end
        end
        Basis_x = [X(:,2) X(:,3)];
        Basis_x_Sym = [Symbol_x(2) Symbol_x(3)];
        for i=1:size(Basis_x,2)
            Data(:,Index)=dX(:,iter).*X(:,1).^2.*Basis_x(:,i);
            Sym_Struct{1,Index}=Symbol_dX(iter,1)*Symbol_x(1,1).^2*Basis_x_Sym(1,i);
            Index=Index+1;
        end
        Basis_x = [X(:,1) X(:,3)];
        Basis_x_Sym = [Symbol_x(1) Symbol_x(3)];
        for i=1:size(Basis_x,2)
            Data(:,Index)=dX(:,iter).*X(:,2).^2.*Basis_x(:,i);
            Sym_Struct{1,Index}=Symbol_dX(iter,1)*Symbol_x(2,1).^2*Basis_x_Sym(1,i);
            Index=Index+1;
        end
        Basis_x = [X(:,1) X(:,2)];
        Basis_x_Sym = [Symbol_x(1) Symbol_x(2)];
        for i=1:size(Basis_x,2)
            Data(:,Index)=dX(:,iter).*X(:,3).^2.*Basis_x(:,i);
            Sym_Struct{1,Index}=Symbol_dX(iter,1)*Symbol_x(3,1).^2*Basis_x_Sym(1,i);
            Index=Index+1;
        end
        %Poly 2 D
        for i=1:size(Basis_d,2)
           for j=i:size(Basis_d,2)
               Data(:,Index)=dX(:,iter).*Basis_d(:,i).*Basis_d(:,j);
               Sym_Struct{1,Index}=Symbol_dX(iter,1)*Basis_d_Sym(1,i)*Basis_d_Sym(1,j);
               Index=Index+1;
           end
        end 
       %Poly 2 U
       for i=1:size(Basis_u,2)
        for j=i:size(Basis_u,2)
           Data(:,Index)=dX(:,iter).*Basis_u(:,i).*Basis_u(:,j);
           Sym_Struct{1,Index}=Symbol_dX(iter,1)*Basis_u_Sym(1,i)*Basis_u_Sym(1,j);
           Index=Index+1;
         end
       end
    end
end

if iter ==3
  if Highest_dXPoly_Order_Guess >= 0
    %Without dX
    Data(:,Index)=dX(:,3);
    Sym_Struct{1,Index}=Symbol_dX(3,1);
    Index=Index+1;
  end
    if Highest_dXPoly_Order_Guess >= 1
    %Physical Knowledge for X3
    Basis_x     = X(:,3);
    Basis_x_Sym = Symbol_x(3);
    Basis_d     = [d(:,1) d(:,3)];
    Basis_d_Sym =[Symbol_d(1) Symbol_d(3)];
    Basis_u     = [u(:,2) u(:,3)];
    Basis_u_Sym =[Symbol_u(2) Symbol_u(3)];
    %Without dX
    Data(:,Index)=dX(:,3);
    Sym_Struct{1,Index}=Symbol_dX(3,1);
    Index=Index+1;
    %Poly 1 X
    for i=1:size(Basis_x,2)
        Data(:,Index)=dX(:,iter).*Basis_x(:,i);
        Sym_Struct{1,Index}=Symbol_dX(iter,1)*Basis_x_Sym(1,i);
        Index=Index+1;
    end
    %Poly 1 D
    for i=1:size(Basis_d,2)
        Data(:,Index)=dX(:,iter).*Basis_d(:,i);
        Sym_Struct{1,Index}=Symbol_dX(iter,1)*Basis_d_Sym(1,i);
        Index=Index+1;
    end
    %Poly 1 U
    for i=1:size(Basis_u,2)
        Data(:,Index)=dX(:,iter).*Basis_u(:,i);
        Sym_Struct{1,Index}=Symbol_dX(iter,1)*Basis_u_Sym(1,i);
        Index=Index+1;
    end
   if Highest_dXPoly_Order_Guess >= 2
        %Poly 2 D
        for i=1:size(Basis_d,2)
           for j=i:size(Basis_d,2)
               Data(:,Index)=dX(:,iter).*Basis_d(:,i).*Basis_d(:,j);
               Sym_Struct{1,Index}=Symbol_dX(iter,1)*Basis_d_Sym(1,i)*Basis_d_Sym(1,j);
               Index=Index+1;
           end
        end 
       %Poly 2 U
       for i=1:size(Basis_u,2)
        for j=i:size(Basis_u,2)
           Data(:,Index)=dX(:,iter).*Basis_u(:,i).*Basis_u(:,j);
           Sym_Struct{1,Index}=Symbol_dX(iter,1)*Basis_u_Sym(1,i)*Basis_u_Sym(1,j);
           Index=Index+1;
         end
       end 
    end
    end
end

if iter == 4
    if Highest_dXPoly_Order_Guess >= 0
    Basis_d     = d(:,4);
    Basis_d_Sym = Symbol_d(4);
        Data(:,Index)=dX(:,iter).*Basis_d(:,1);
        Sym_Struct{1,Index}=Symbol_dX(iter,1)*Basis_d_Sym(1,1);
        Index=Index+1;
    end 
    if Highest_dXPoly_Order_Guess >= 1
        %Physical Knowledge for X4
        Basis_x     = [X(:,1) X(:,3) X(:,4)]; 
        Basis_x_Sym = [Symbol_x(1) Symbol_x(3) Symbol_x(4)];
        Basis_d     = d(:,4);
        Basis_d_Sym = Symbol_d(4);
        Basis_u     = u(:,2);
        Basis_u_Sym = Symbol_u(2);
        %Poly 1 X
        for i=1:size(Basis_x,2)
            Data(:,Index)=dX(:,iter).*Basis_x(:,i);
            Sym_Struct{1,Index}=Symbol_dX(iter,1)*Basis_x_Sym(1,i);
            Index=Index+1;
        end
        %Poly 1 D
        for i=1:size(Basis_d,2)
            Data(:,Index)=dX(:,iter).*Basis_d(:,i);
            Sym_Struct{1,Index}=Symbol_dX(iter,1)*Basis_d_Sym(1,i);
            Index=Index+1;
        end
        %Poly 1 U
        for i=1:size(Basis_u,2)
            Data(:,Index)=dX(:,iter).*Basis_u(:,i);
            Sym_Struct{1,Index}=Symbol_dX(iter,1)*Basis_u_Sym(1,i);
            Index=Index+1;
        end
    end 
    if Highest_dXPoly_Order_Guess >= 2
        %Poly 2 X
        for i=1:size(Basis_x,2)
            for j=i:size(Basis_x,2)
                Data(:,Index)=dX(:,iter).*Basis_x(:,i).*Basis_x(:,j);
                Sym_Struct{1,Index}=Symbol_dX(iter,1)*Basis_x_Sym(1,i)*Basis_x_Sym(1,j);
                Index=Index+1;
            end
        end
        Basis_x = [X(:,1) X(:,3)];
        Basis_x_Sym = [Symbol_x(1) Symbol_x(3)];
        for i=1:size(Basis_x,2)
            Data(:,Index)=dX(:,iter).*X(:,4).^2.*Basis_x(:,i);
            Sym_Struct{1,Index}=Symbol_dX(iter,1)*Symbol_x(4,1).^2*Basis_x_Sym(1,i);
            Index=Index+1;
        end
        Basis_x = [X(:,1) X(:,4)];
        Basis_x_Sym = [Symbol_x(1) Symbol_x(4)];
        for i=1:size(Basis_x,2)
            Data(:,Index)=dX(:,iter).*X(:,3).^2.*Basis_x(:,i);
            Sym_Struct{1,Index}=Symbol_dX(iter,1)*Symbol_x(3,1).^2*Basis_x_Sym(1,i);
            Index=Index+1;
        end
        Basis_x = [X(:,3) X(:,4)];
        Basis_x_Sym = [Symbol_x(3) Symbol_x(4)];
        for i=1:size(Basis_x,2)
            Data(:,Index)=dX(:,iter).*X(:,1).^2.*Basis_x(:,i);
            Sym_Struct{1,Index}=Symbol_dX(iter,1)*Symbol_x(1,1).^2*Basis_x_Sym(1,i);
            Index=Index+1;
        end
        %Poly 2 D
        for i=1:size(Basis_d,2)
           for j=i:size(Basis_d,2)
               Data(:,Index)=dX(:,iter).*Basis_d(:,i).*Basis_d(:,j);
               Sym_Struct{1,Index}=Symbol_dX(iter,1)*Basis_d_Sym(1,i)*Basis_d_Sym(1,j);
               Index=Index+1;
           end
        end 
       %Poly 2 U
       for i=1:size(Basis_u,2)
        for j=i:size(Basis_u,2)
           Data(:,Index)=dX(:,iter).*Basis_u(:,i).*Basis_u(:,j);
           Sym_Struct{1,Index}=Symbol_dX(iter,1)*Basis_u_Sym(1,i)*Basis_u_Sym(1,j);
           Index=Index+1;
         end
       end
    end 
    end
end