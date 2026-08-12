%% This function will generate the simulation data of an ODE function.
%You could determine the noise level by input variable "Noise". If you do
%not want any noise, set noise to zero. Please indicate whether your ODE
%function have control input, if the answer is yes, please set the
%"Control" as 1.

% Last Update: 2019/04/21
% Coded By: K

function [dData,Data]=Get_Sim_Data_GHLettuce(ODE,x0,u,d,tspan,Control,Disturbance,Noise,Shuffle)
%% Get the size of the state and control
[N1,M1]=size(x0);
[N2,M2]=size(u);

%% Extract the Data from Struct to Double
% u later

%% Get simulation data by simulating the system using ODE45
% Determine the left hand side derivative
    x_list(1,:)=x0; %transpose 
    dx_list(1,:)=ODE(0,x_list(1,:)',u(1,:)',d(1,:)'); 
    for i=2:length(u) 
        [t_1,x_1] = ode45(@(t_1,x_1)ODE(t_1,x_1,u(i-1,:)',d(i-1,:)'),tspan(1,i-1:i),x0'); 
        %, u and disturbance becomes integration, tapi ode45 should be
        %column vector. x_1 should be row vector
        x_list(i,:)  =x_1(end,:);
        dx_list(i,:) =ODE(0,x_list(i,:)',u(i,:)',d(i,:)');
        x0=x_list(i,:)'; 
    end

%% Add some noise to the system
for i=1:N1
    Data(:,i)=x_list(:,i)+Noise*randn(size(x_list(:,i)));
end
%
for i=1:N1
    dData(:,i)=dx_list(:,i)+Noise*randn(size(dx_list(:,i)));
end

%% Shuffle the data
if Shuffle==1
    Sequence=randperm(size(Data,1));
    Data_Es=Data_Es(Sequence,:);
    dData=dData(Sequence,:);
end