% This library is coded for greenhouse lettuce example
% Last Updated: 29/05/2023
% Coded By: Fakhira
function [Data,Sym_Struct]= SINDyLib_GHLettuce_PhysML_Best(X,dX,u,d,iter,Highest_Poly_Order, Highest_dXPoly_Order)
% 
% load('SINDyPI_Greenhouselettuce_PhysML_Val_TRUE.mat');
% X = Data;
% dX = dData;
% iter =2;
% Highest_Poly_Order = 2;
% Highest_dXPoly_Order=0;

% First get the size of the X matrix, determine the data length and the number of variables we have.
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
Index=1;  %for iteration to pass to another cell

%% Candidate Libray Order 1 %for Reset

%% Candidate Libray Order 1 %for Reset

%STATE 1
if iter ==1
if Highest_Poly_Order == 1
    %Physical Knowledge for X1
    Basis_x     = [X(:,1) X(:,2) X(:,3)]; 
    Basis_x_Sym = [Symbol_x(1) Symbol_x(2) Symbol_x(3)];
    Basis_d     = d(:,1);
    Basis_d_Sym = Symbol_d(1);
    %Poly 1 X
    for i=1:size(Basis_x,2)
        Data(:,Index)=Basis_x(:,i);
        Sym_Struct{1,Index}=Basis_x_Sym(1,i);
        Index=Index+1;
    end
    %Poly 1 D
    for i=1:size(Basis_d,2)
        Data(:,Index)=Basis_d(:,i);
        Sym_Struct{1,Index}=Basis_d_Sym(1,i);
        Index=Index+1;
    end
    %Poly 1 XD
    for i=1:size(Basis_x,2)
        for j=1:size(Basis_d,2)
            Data(:,Index)=Basis_x(:,i).*Basis_d(:,j);
            Sym_Struct{1,Index}=Basis_x_Sym(1,i)*Basis_d_Sym(1,j);
            Index=Index+1;
        end
    end
end
end 

if iter ==2
%STATE 2
if Highest_Poly_Order == 1
    %Physical Knowledge for X2
    Basis_x     = [X(:,1) X(:,2) X(:,3)]; 
    Basis_x_Sym = [Symbol_x(1) Symbol_x(2) Symbol_x(3)];
    Basis_d     = [d(:,1) d(:,2)];
    Basis_d_Sym =[Symbol_d(1) Symbol_d(2)];
    Basis_u     = [u(:,1) u(:,2)];
    Basis_u_Sym =[Symbol_u(1) Symbol_u(2)];
 %Without XD-XUD
    %Poly 1 X
    for i=1:size(Basis_x,2)
        Data(:,Index)=Basis_x(:,i);
        Sym_Struct{1,Index}=Basis_x_Sym(1,i);
        Index=Index+1;
    end
    %Poly 1 D
    for i=1:size(Basis_d,2)
        Data(:,Index)=Basis_d(:,i);
        Sym_Struct{1,Index}=Basis_d_Sym(1,i);
        Index=Index+1;
    end
    %Poly 1 U
    for i=1:size(Basis_u,2)
        Data(:,Index)=Basis_u(:,i);
        Sym_Struct{1,Index}=Basis_u_Sym(1,i);
        Index=Index+1;
    end
%***XD Deleted****
    %Poly 1 XU
    for i=1:size(Basis_x,2)
        for j=1:size(Basis_u,2)
            Data(:,Index)=Basis_x(:,i).*Basis_u(:,j);
            Sym_Struct{1,Index}=Basis_x_Sym(1,i)*Basis_u_Sym(1,j);
            Index=Index+1;
        end
    end
    %Poly 1 UD
    for i=1:size(Basis_u,2)
        for j=1:size(Basis_d,2)
            Data(:,Index)=Basis_u(:,i).*Basis_d(:,j);
            Sym_Struct{1,Index}=Basis_u_Sym(1,i)*Basis_d_Sym(1,j);
            Index=Index+1;
        end
    end
%***XUD Deleted****
end
end

if iter==3
%STATE 3
if Highest_Poly_Order ==1
    %Physical Knowledge for X3
    Basis_x     = X(:,3);
    Basis_x_Sym = Symbol_x(3);
    Basis_d     = [d(:,1) d(:,3)];
    Basis_d_Sym =[Symbol_d(1) Symbol_d(3)];
    Basis_u     = [u(:,2) u(:,3)];
    Basis_u_Sym =[Symbol_u(2) Symbol_u(3)];
    %Poly 1 X
    for i=1:size(Basis_x,2)
        Data(:,Index)=Basis_x(:,i);
        Sym_Struct{1,Index}=Basis_x_Sym(1,i);
        Index=Index+1;
    end
    %Poly 1 D
    for i=1:size(Basis_d,2)
        Data(:,Index)=Basis_d(:,i);
        Sym_Struct{1,Index}=Basis_d_Sym(1,i);
        Index=Index+1;
    end
    %Poly 1 U
        Data(:,Index)=Basis_u(:,2);
        Sym_Struct{1,Index}=Basis_u_Sym(1,2);
        Index=Index+1;
    %BEST lib
    Data(:,Index)=Basis_d(:,2).*Basis_u(:,1);
    Sym_Struct{1,Index}=Basis_d_Sym(1,2)*Basis_u_Sym(1,1);
    Index=Index+1;

    Data(:,Index)=Basis_x(:,1).*Basis_u(:,1);
    Sym_Struct{1,Index}=Basis_x_Sym(1,1)*Basis_u_Sym(1,1);
    Index=Index+1;


%     %Poly 1 XD
%     for i=1:size(Basis_x,2)
%         for j=1:size(Basis_d,2)
%             Data(:,Index)=Basis_x(:,i).*Basis_d(:,j);
%             Sym_Struct{1,Index}=Basis_x_Sym(1,i)*Basis_d_Sym(1,j);
%             Index=Index+1;
%         end
%     end 
%     %Poly 1 XU
%     for i=1:size(Basis_x,2)
%         for j=1:size(Basis_u,2)
%             Data(:,Index)=Basis_x(:,i).*Basis_u(:,j);
%             Sym_Struct{1,Index}=Basis_x_Sym(1,i)*Basis_u_Sym(1,j);
%             Index=Index+1;
%         end
%     end
%     %Poly 1 UD
%     for i=1:size(Basis_u,2)
%         for j=1:size(Basis_d,2)
%             Data(:,Index)=Basis_u(:,i).*Basis_d(:,j);
%             Sym_Struct{1,Index}=Basis_u_Sym(1,i)*Basis_d_Sym(1,j);
%             Index=Index+1;
%         end
%     end 
%     %Poly 1 XUD
%     for i=1:size(Basis_x,2)
%         for j=1:size(Basis_u,2)
%             for k=1:size(Basis_d,2)
%                 Data(:,Index)=Basis_x(:,i).*Basis_u(:,j).*Basis_d(:,k);
%                 Sym_Struct{1,Index}=Basis_x_Sym(1,i)*Basis_u_Sym(1,j)*Basis_d_Sym(1,k);
%                 Index=Index+1;
%             end
%         end
%     end
end
end

if iter ==4
%STATE 4
if Highest_Poly_Order ==1
    %Physical Knowledge for X4
    Basis_x     = [X(:,1) X(:,3) X(:,4)]; 
    Basis_x_Sym = [Symbol_x(1) Symbol_x(3) Symbol_x(4)];
    Basis_d     = d(:,4);
    Basis_d_Sym = Symbol_d(4);
    Basis_u     = u(:,2);
    Basis_u_Sym = Symbol_u(2);
%Without XD-UD
    %Poly 1 X
    for i=1:size(Basis_x,2)
        Data(:,Index)=Basis_x(:,i);
        Sym_Struct{1,Index}=Basis_x_Sym(1,i);
        Index=Index+1;
    end
    %Poly 1 D
    for i=1:size(Basis_d,2)
        Data(:,Index)=Basis_d(:,i);
        Sym_Struct{1,Index}=Basis_d_Sym(1,i);
        Index=Index+1;
    end
    %Poly 1 U
    for i=1:size(Basis_u,2)
        Data(:,Index)=Basis_u(:,i);
        Sym_Struct{1,Index}=Basis_u_Sym(1,i);
        Index=Index+1;
    end
%***XD Deleted****
    %Poly 1 XU
    for i=1:size(Basis_x,2)
        for j=1:size(Basis_u,2)
            Data(:,Index)=Basis_x(:,i).*Basis_u(:,j);
            Sym_Struct{1,Index}=Basis_x_Sym(1,i)*Basis_u_Sym(1,j);
            Index=Index+1;
        end
    end
%***UD Deleted****
    %Poly 1 XUD
    for i=1:size(Basis_x,2)
        for j=1:size(Basis_u,2)
            for k=1:size(Basis_d,2)
                Data(:,Index)=Basis_x(:,i).*Basis_u(:,j).*Basis_d(:,k);
                Sym_Struct{1,Index}=Basis_x_Sym(1,i)*Basis_u_Sym(1,j)*Basis_d_Sym(1,k);
                Index=Index+1;
            end
        end
    end
end
end

if iter ==1
    %Use Derivative Best Guess using previous pretest
    if Highest_dXPoly_Order == 1
        Data(:,Index)=dX(:,iter).*X(:,1);
        Sym_Struct{1,Index}=Symbol_dX(iter,1)*Symbol_x(1,1);
        Index=Index+1;
    
        Data(:,Index)=dX(:,iter).*X(:,1).^2;
        Sym_Struct{1,Index}=Symbol_dX(iter,1)*Symbol_x(1,1).^2;
        Index=Index+1;
    end 
end

if iter ==2
    if Highest_dXPoly_Order == 2
        Data(:,Index)=dX(:,iter).*X(:,3); %THE BEST
        Sym_Struct{1,Index}=Symbol_dX(iter,1)*Symbol_x(3,1);  
    end
end

%dX
if iter==3
    if Highest_dXPoly_Order == 3
        %Without dX
        Data(:,Index)=dX(:,3);
        Sym_Struct{1,Index}=Symbol_dX(3,1);
        Index=Index+1;
    end 
end

if iter ==4
    if Highest_dXPoly_Order == 4
        Data(:,Index)=dX(:,iter).*d(:,4);
        Sym_Struct{1,Index}=Symbol_dX(iter,1)*Symbol_d(4,1);
        Index=Index+1;
    end 
end