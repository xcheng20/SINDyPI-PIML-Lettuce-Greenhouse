%% This function will calculate the l2 norm error between the discovered system and actual test data
% Last Updated: 2019/04/22
% Coded By: K

function [Score]=Get_ScoreGHLettuce(dData,Data,u,d,Control,Disturbance,tspan,state0,Shuffle)
% Get the ODE simulation result
Noise=0;
% Using ODEs to get the simulation data
%[dData_Es,Data_Es]=Get_Sim_Data(@(t,z)Sindy_ODE_RHS(t,z,u),state0,u,tspan,Noise,Control,Shuffle);

% Get the score of the result
%[n,m]=size(Data_Es);
%Score=sum(norm(dData_test(1:n,:)-dData_Es))+sum(norm(Data_test(1:n,:)-Data_Es));

%% Using one step prediction to get the simulation data
dData_Es=Sindy_ODE_RHS_GHLettuce(0,Data',u',d')';   %kenapa ini gak u_test?

% Get the score of the result
% Error = dData-dData_Es;
% Score=sum(norm(Error));
Score=sum(norm(dData-dData_Es)/norm(dData));





