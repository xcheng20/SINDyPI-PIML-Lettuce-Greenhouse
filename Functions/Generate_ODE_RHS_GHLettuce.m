%% This function will generate the ODE function based on calculation result of the symbolic expression.
% This file has been modified to perform the parallel computing.
% Last Update: 2019/05/15
% Coded By: K

%% Test Data
% load('Test_Data\ODEs_DouPen.mat')
% 
% Sindy_ODEs=ODEs
% var_num_state=1;
% var_num_control=2;
% var_num_disturbance=1;

%% Function
function Generate_ODE_RHS_DW(Sindy_ODEs,var_num_state,var_num_control,var_num_disturbance)
x=sym('x',[var_num_state,1]);
u=sym('u',[var_num_control,1]);
d=sym('d',[var_num_disturbance,1]);
syms t
f= matlabFunction(Sindy_ODEs,'File','Sindy_ODE_RHS_GHLettuce','Optimize',true,'Vars',{t,x,u,d});










