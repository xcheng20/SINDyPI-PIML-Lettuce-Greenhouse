% SINDy-PI Test to Lettuce Greenhouse Model
%Date: 28/05/2023
%Code by: Fakhira

% Close all, clear all, clc
close all;clear all; clc;
[status,msg] = mkdir('Results');
addpath('Functions')
set(0,'defaulttextInterpreter','latex')

%% ===============================================PREPARING FOR SIMULATION========================================================
% Simulate the greenhouse and gather the simulation data
% Model parameters
ops.nx  = 4;    % #states
ops.nu  = 3;    % #controllable inputs
ops.nd  = 4;    % #disturbance inputs

% simulation parameters
c         = 86400;                          % number of seconds in a day (24*60*60)
nDays     = 40;                             % #days in simulation (max=320)
ops.h     = 15*60;                          % sample period in seconds: #minutes*60 = dt per hour
ops.L     = nDays*c;                        % final time simulation
ops.t     = 0:ops.h:ops.L;                  % initial time vector  %units in seconds
ops.N     = length(ops.t);                  % number of samples in initial time vector
ops.N0    = 1;     

% control signals
x0          = [0.0035; 1e-03 ;15; 0.008];      % initial state (biomass)
x0_test     = [0.0040; 3e-03; 12; 0.009];      
x0_val      = [0.0042; 2e-03; 13; 0.007];       

u_0_40Days      = load('bin/data/control signals/u_0_40Days.mat');
u               = u_0_40Days.u_0_40Days;           % make a cell
u               = u';                              % transpose the u
u               = u(2:3842,:);
u_40_80Days     = load('bin\data\control signals\u_40_80Days.mat');
u_test          = u_40_80Days.u_40_80Days;           % make a cell
u_test          = u_test';                           % transpose the u     
u_test          = u_test(2:3842,:);
u_80_120Days    = load('bin\data\control signals\u_80_120Days.mat');
u_val           = u_80_120Days.u_80_120Days;           % make a cell
u_val           = u_val';                           % transpose the u
u_val           = u_val(2:3842,:);
u_new           = [u ; u_test ; u_val]; 

%disturbances
d_0_40Days  = load("bin/data/disturbances/d_0_40Days.mat");
d           = d_0_40Days.d_0_40Days;                            
d           = d';                            
% convert CO2 from ppm.10^3
d(:,2)      = co2ppm2dens(d(:,3), d(:,2))*1e3; 
d(:,4)      = rh2vaporDens(d(:,3), d(:,4));
d           = d(1:3841,:);

d_40_80Days  = load("bin\data\disturbances\d_40_80Days.mat");
d_test       = d_40_80Days.d_40_80Days;                           
d_test       = d_test';                           
% convert CO2 from ppm.10^3
d_test(:,2)    = co2ppm2dens(d_test(:,3), d_test(:,2))*1e3;
d_test(:,4)    = rh2vaporDens(d_test(:,3), d_test(:,4));
d_test         = d_test(1:3841,:);

d_80_120Days  = load("bin\data\disturbances\d_80_120Days.mat");
d_val         = d_80_120Days.d_80_120Days;                           
d_val         = d_val';                             
% convert CO2 from ppm.10^3
d_val(:,2)    = co2ppm2dens(d_val(:,3), d_val(:,2))*1e3;
d_val(:,4)    = rh2vaporDens(d_val(:,3), d_val(:,4));
d_val         = d_val(1:3841,:);

d_new         = [d; d_test ; d_val]; 

%Redefine control signals and disturbances for 30 days
u       = u_new(1:ops.N,:);
u_test  = u_new(ops.N+1:ops.N+ops.N,:);
u_test2 = u_new(ops.N+ops.N+1:ops.N+ops.N+ops.N,:);

d       = d_new(1:ops.N,:);
d_test  = d_new(ops.N+1:ops.N+ops.N,:);
d_test2 = d_new(ops.N+ops.N+1:ops.N+ops.N+ops.N,:);

% Test Plot u and d
CPUTime   = zeros(1,ops.N);
% x(1)    = x.dw
% x(2)    = x.co2In
% x(3)    = x.tempIn
% x(4)    = x.vapIn

% d(1)    = d.rad;
% d(2)    = d.co2Out
% d(3)    = d.tempOut
% d(4)    = d.vapOut

% u(1)    = u.co2
% u(2)    = u.vent
% u(3)    = u.heat

%% RUN THE ODE FILES AND GATHER THE SIMULATION DATA

% Add noise level
% Noise = [0;1e-12;1e-11;1e-10;1e-9;1e-8;1e-7];
Noise=1e-7; %Choose noise level 

% Shuffle the data
Shuffle = 0;

for n=1:length(Noise)
%Perform Simulation for the Training
[dData,Data]=Get_Sim_Data_GHLettuce(@(t,x,u,d)GHLettuce_ODE(t,x,u,d),x0,u,d,ops.t,ops.nu,ops.nd,Noise(n,1),Shuffle);

%Perform Simulation for the Test Data
[dData_test,Data_test]=Get_Sim_Data_GHLettuce(@(t,x,u,d)GHLettuce_ODE(t,x,u,d),x0_test,u_test,d_test,ops.t,ops.nu,ops.nd,Noise(n,1),Shuffle);

%% Plot Data
% figure(1)
% PlotResults_GHLettuce_train
% figure(2)
% PlotResults_GHLettuce_test
% figure(3)
% PlotResults_GHLettuce_val

%% SINDy-PI Preparation
% Now perform sparse regression of non-linear dynamics
% Get the number of states we have
[dtat_length,n_state]=size(Data);

% Define the control input
n_control=3;

% Define the disturbance
n_disturbance=4;

% Choose whether you want to display actual ODE or not
disp_actual_ode=1;

% If the ODEs you want to display is the actual underlyting dynamics of the
% system, please set actual as 1
actual=1;

%% Print the actual ODE we try to discover %investigate this one
digits(4)
Print_ODEs = Print_ODEs_GHLettuce(@(t,x,u,d)GHLettuce_ODE(t,x,u,d),n_state,n_control,n_disturbance,disp_actual_ode,actual);

%% Create symbolic states
dx=sym('dx',[n_state,1]);

% Theta Library Preparation
lam=[1e-16];

N_iter=20;  
disp=1;
NormalizeLib=1; 
tic

%% =============================================== SPARSE REGRESSION ========================================================
%RESULT AFTER MANY ITERATIONS
%LHS Guess
dXPoly_Order_Guess=[1;2;3;4];
%RHS Library
Poly_Order=[1;1;1;1]; %design for second order
dXPoly_Order=[1;2;3;4]; %design for longer model structure

    %Sparse Regression from X1 to X4
    for iter=1:n_state
        fprintf('\n \n Calculating the %i expression...\n',iter)
        % Change the library base on what equation you want to work on
        Highest_dXPoly_Order_Guess= dXPoly_Order_Guess(iter);
        Highest_dXPoly_Order= dXPoly_Order(iter); 
        Highest_Poly_Order = Poly_Order(iter); 

        % According to the previous parameter generate the left hand side guess
        [LHS_Data,LHS_Sym]=GuessLib_GHLettuce_PhysML_Best(Data,dData,u,d,iter,Highest_dXPoly_Order_Guess);  
        %Store the result
        LHS_Sym_all{:,iter}=LHS_Sym; 

        %Generate RHS Library %Test or Best
        [SINDy_Data,SINDy_Struct]= SINDyLib_GHLettuce_PhysML_Test(Data,dData,u,d,iter,Highest_Poly_Order, Highest_dXPoly_Order);

        % Store the last library result
        SINDy_Struct_all{:,iter}=SINDy_Struct; 
        SINDy_Data_all{:,iter}=SINDy_Data; 
           
        % Run the for loop and try all the left hand guess
        for i=1:length(LHS_Sym)
             if iter==1 && i==1
                Xi=cell(n_state,length(LHS_Sym),length(lam));
                ODE=cell(n_state,length(LHS_Sym),length(lam));
                ODEs=cell(n_state,length(LHS_Sym),length(lam));
            end
            
            % Print the left hand side that we are testing
            fprintf('\t Testing the left hand side as %s:\n',cell2sym(LHS_Sym(1,i)));
            
            % Exclude the guess from SINDy library
            [RHS_Data,RHS_Struct]=ExcludeGuess(SINDy_Data,SINDy_Struct,LHS_Sym{i});
            
            for j=1:length(lam)
                fprintf('\n\t\t Testing the lambda as %d... \n',lam(j,1))
                % Select the sparse threshold
                lambda=lam(j);
    
                % Perform the sparse regression problem, Datasets for only
                % Model selection
                [Xi_lambda{iter,i,j},ODE{iter,i,j}]=sparsifyDynamics_GHLettuce(RHS_Data,LHS_Data(:,i),LHS_Sym{i},lambda,N_iter,RHS_Struct,disp,NormalizeLib);
    
                %% Perform sybolic calculation and solve for dX
                Eqn=LHS_Sym{i}==ODE{iter,i,j};
                digits(4)
                ODE_Guess=(vpasolve(Eqn,dx(iter))); %Simplify does not give any difference in this case
            
                % Print the discovered ODE
                fprintf(strcat('\t The corresponding ODE we found is: ',char(dx(iter,1)),'=',char(simplify(ODE_Guess)),'\n \n'));   
                % It is ok to use simplify since it simplify the equation but dont put it together with vpa  
                
                % Store the result
                ODEs{iter,i,j}=ODE_Guess; 
            end
        end
    end 
    %% Now generate the ODE function file and test the accuracy of the
    % identified system
    
    fprintf('\v Start calculating the best model that could represent the training data...\n \n')
    
    for iter=1:n_state
        % Print which expression are you working on
        fprintf('\t Calculating the best model for the %d expression...\n',iter)
    
     for i=1:length(LHS_Sym_all{1,iter})
            % Print the process
            fprintf('\t Calculating the score of previously found ODE on the test data, %d %% finished. \n',round((i/length(LHS_Sym))*100))
    
            for j=1:length(lam)
                %If the previous ODE is 0, set the score as NaN, else calculate
                %it.
                if isempty(ODEs{iter,i,j})   %must be integer or logical values
                    ODE_Not_Exist=1;
                    %fprintf('\n ODE at X, Y: %d, %d', i, j)
                    Score(iter,i,j)=NaN;
                else
                %Generate the ODE file
                Generate_ODE_RHS_GHLettuce(ODEs{iter,i,j},n_state,n_control,n_disturbance)
    
                % Calculate the accuracy of the file 
                Score(iter,i,j)=Get_ScoreGHLettuce(dData_test(:,iter),Data_test,u_test,d_test,ops.nu,ops.nd,ops.t,x0_test,Shuffle);
                end 
            end 
            %Get the best lambda, after all of the LHS_Sym
            [minVal1(iter,i),minIndex1(iter,i)]=min(Score(iter,i,:));  %minIndex1 = lambda location, minIndex2= LHS guess, val=value
     end
         
        % Get the best score of all LHS and use this ODE file
        [minVal2(iter,1),minIndex2(iter,1)]=min(minVal1(iter,1:length(LHS_Sym_all{1,iter})));   
        
        % Store the best ODE
        ODE_Best(iter,1)=ODEs{iter,minIndex2(iter,1),minIndex1(iter,minIndex2(iter,1))};  %{LHS_Guess, lambda}
        
        % Print the Result
        fprintf('\n\n\n\t The SINDy-PI discovered Best ODE for the %d expression is:\n',iter)
        fprintf('\t %s = %s \n\n\n',char(dx(iter)),char((ODE_Best(iter)))')
    end
   
    

    %% Now generate this best guess ODE and print its result
    disp_best=1;
    if disp_best==1
        fprintf('\n\n\n\v The SINDy-PI discovered Best ODE for the whole system is:\n')
        digits(4)
        for iter=1:n_state
            fprintf(strcat('\v ******\v\t',char(dx(iter)),'=',char(simplify(ODE_Best(iter,1))),'\n'));
        end
    end
    
    % Also Print the actual ODE for comparison
    digits(4)
    fprintf('\n\n\n')
    Print_ODEs_GHLettuce(@(t,x,u,d)GHLettuce_ODE(t,x,u,d),n_state,n_control,n_disturbance,disp_actual_ode,actual);

    
    % Get the simulation result
    % Generate the ODE file after training
    fprintf('\n\n\n\v Generating the Best Model for comparison for...')
    Generate_ODE_RHS_GHLettuce(ODE_Best(:,1),n_state,n_control,n_disturbance);
    
    % Gather the new simulation data from the best ODE model (same dataset)
    [dData_True_val,Data_True_val]=Get_Sim_Data_GHLettuce(@(t,x,u,d)GHLettuce_ODE(t,x,u,d),x0_val,u_val,d_val,ops.t,ops.nu,ops.nd,Noise(n,1),Shuffle);
    [dData_Es_val,Data_Es_val]=Get_Sim_Data_GHLettuce(@(t,x,u,d)Sindy_ODE_RHS_GHLettuce(t,x,u,d),x0_val,u_val,d_val,ops.t,ops.nu,ops.nd,Noise(n,1),Shuffle);
    
    % Prediction Error 
    dScore_1_val_Best=(norm(dData_True_val(:,1)-dData_Es_val(:,1)))/norm(dData_True_val(:,1));
    dScore_2_val_Best=(norm(dData_True_val(:,2)-dData_Es_val(:,2)))/norm(dData_True_val(:,2));
    dScore_3_val_Best=(norm(dData_True_val(:,3)-dData_Es_val(:,3)))/norm(dData_True_val(:,3));
    dScore_4_val_Best=(norm(dData_True_val(:,4)-dData_Es_val(:,4)))/norm(dData_True_val(:,4));
    
    % Store the result of Simulation
    dScore_all_Best=[dScore_1_val_Best dScore_2_val_Best dScore_3_val_Best dScore_4_val_Best];
    
    %% Store individual result
    %STATE 1
    Lambda_X1_Best= lam(minIndex1(1,minIndex2(1,1)));
    LHS_X1_Best= LHS_Sym_all{1,1}{1,minIndex2(1,1)};
    dx1_1 = {'dx1*x1','dx1*x1^2'};
    Sym_Struct_X1_Best= strjoin(dx1_1);
    dScore_X1_Best = dScore_1_val_Best;
    ODEs_X1_Best = ODE_Best(1,1);
    X1_BEST = [string(LHS_X1_Best) Sym_Struct_X1_Best Lambda_X1_Best string(dScore_X1_Best) string(ODEs_X1_Best)];
     ODEs_X1_Best(n,1)= ODE_Best(1,1);
    dScore_X1_Noise(1,n)=dScore_1_val_Best;
    Data_Es_val_1(:,n)=Data_Es_val(:,1);
    dData_Es_val_1(:,n)=dData_Es_val(:,1);

    %STATE 2
    Lambda_X2_Best= lam(minIndex1(2,minIndex2(2,1)));
    LHS_X2_Best= LHS_Sym_all{1,2}{1,minIndex2(2,1)};
    Sym_Struct_X2_Best= SINDy_Struct_all{1,2}(1,end);
    dScore_X2_Best = dScore_2_val_Best;
    ODEs_X2_Best = ODE_Best(2,1);
    X2_BEST = [string(LHS_X2_Best) Sym_Struct_X2_Best Lambda_X2_Best string(dScore_X2_Best) string(ODEs_X2_Best)];
    ODEs_X2_Best(n,1) = ODE_Best(2,1);
    dScore_X2_Noise(1,n)=dScore_2_val_Best;
    Data_Es_val_2(:,n)=Data_Es_val(:,2);
    dData_Es_val_2(:,n)=dData_Es_val(:,2);

    %STATE 3
    Lambda_X3_Best= lam(minIndex1(3,minIndex2(3,1)));
    LHS_X3_Best= LHS_Sym_all{1,3}{1,minIndex2(3,1)};
    Sym_Struct_X3_Best= SINDy_Struct_all{1,3}(1,end);
    dScore_X3_Best = dScore_3_val_Best;
    ODEs_X3_Best= ODE_Best(3,1);
    X3_BEST = [string(LHS_X3_Best) Sym_Struct_X3_Best Lambda_X3_Best string(dScore_X3_Best) string(ODEs_X3_Best)];
    ODEs_X3_Best(n,1)= ODE_Best(3,1);
    dScore_X3_Noise(1,n)=dScore_3_val_Best;
    Data_Es_val_3(:,n)=Data_Es_val(:,3);
    dData_Es_val_3(:,n)=dData_Es_val(:,3);

    % STATE 4
    Lambda_X4_Best= lam(minIndex1(4,minIndex2(4,1)));
    LHS_X4_Best= LHS_Sym_all{1,4}{1,minIndex2(4,1)};
    Sym_Struct_X4_Best= SINDy_Struct_all{1,4}(1,end);
    dScore_X4_Best = dScore_4_val_Best;
    ODEs_X4_Best= ODE_Best(4,1);
    X4_BEST = [Noise(n,1) string(LHS_X4_Best) Sym_Struct_X4_Best Lambda_X4_Best string(dScore_X4_Best) string(ODEs_X4_Best)];
    ODEs_X4_Best(n,1)= ODE_Best(4,1);
    dScore_X4_Noise(1,n)=dScore_4_val_Best;
    Data_Es_val_4(:,n)=Data_Es_val(:,4);
    dData_Es_val_4(:,n)=dData_Es_val(:,4);

    sprintf('\n Finished Simulating The Noise %d...',n) 

    % Show process
    fprintf('\n Simulation finished, plotting the result...... \n')
    close all
    % Create the new directory to save the plot
    % [fld_status, fld_msg, fld_msgID]=mkdir('Figures');
    
    figure(4)
    plot(ops.t/c,Data_True_val(:,1),'linewidth',0.5,'Color','black')
    hold on
    plot(ops.t/c,Data_Es_val(:,1),'linewidth',0.5,'linestyle','--','color','blue')
    title('Validation')
    ylabel('$x_1$ (g/m$^2$)','fontsize',12,'interpreter','latex')
    xlabel('Time (days)','fontsize',13,'interpreter','latex')
    legend('True Prediction','Estimated Prediction')
    
    figure(5)
    plot(ops.t/c,dData_True_val(:,1),'linewidth',0.5,'Color','black')
    hold on
    plot(ops.t/c,dData_Es_val(:,1),'linewidth',0.5,'linestyle','--','color','blue')
    title('Validation')
    ylabel('$\dot{x_1}$','fontsize',12,'interpreter','latex')
    xlabel('Time (days)','fontsize',13,'interpreter','latex')
    legend('True Prediction','Estimated Prediction')
    
    figure(6)
    plot(ops.t/c,Data_True_val(:,2),'linewidth',0.5,'Color','black')
    hold on
    plot(ops.t/c,Data_Es_val(:,2),'linewidth',0.5,'linestyle','--','color','blue')
    title('Validation')
    ylabel('$x_2$ (kg/m$^3$)','fontsize',12,'interpreter','latex')
    xlabel('Time (days)','fontsize',13,'interpreter','latex')
    legend('True Prediction','Estimated Prediction')
    
    figure(7)
    plot(ops.t/c,dData_True_val(:,2),'linewidth',0.5,'Color','black')
    hold on
    plot(ops.t/c,dData_Es_val(:,2),'linewidth',0.5,'linestyle','--','color','blue')
    title('Validation')
    ylabel('$\dot{x_2}$','fontsize',12,'interpreter','latex')
    xlabel('Time (days)','fontsize',13,'interpreter','latex')
    legend('True Prediction','Estimated Prediction')
    
    figure(8)
    plot(ops.t/c,Data_True_val(:,3),'linewidth',0.5,'Color','black')
    hold on
    plot(ops.t/c,Data_Es_val(:,3),'linewidth',0.5,'linestyle','--','color','blue')
    title('Validation')
    ylabel('$x_3$ ($^\circ$C)','fontsize',12,'interpreter','latex')
    xlabel('Time (days)','fontsize',13,'interpreter','latex')
    legend('True Prediction','Estimated Predictio')
    
    figure(9)
    plot(ops.t/c,dData_True_val(:,3),'linewidth',0.5,'Color','black')
    hold on
    plot(ops.t/c,dData_Es_val(:,3),'linewidth',0.5,'linestyle','--','color','blue')
    title('Validation')
    ylabel('$\dot{x_3}$','fontsize',12,'interpreter','latex')
    xlabel('Time (days)','fontsize',13,'interpreter','latex')
    legend('True Prediction','Estimated Prediction')
    
    figure(10)
    plot(ops.t/c,Data_True_val(:,4),'linewidth',0.5,'Color','black')
    hold on
    plot(ops.t/c,Data_Es_val(:,4),'linewidth',0.5,'linestyle','--','color','blue')
    title('Validation')
    ylabel('$x_4$ (kg/m$^3$)','fontsize',12,'interpreter','latex')
    xlabel('Time (days)','fontsize',13,'interpreter','latex')
    legend('True Prediction','Estimated Prediction')
    
    figure(11)
    plot(ops.t/c,dData_True_val(:,4),'linewidth',0.5,'Color','black')
    hold on
    plot(ops.t/c,dData_Es_val(:,4),'linewidth',0.5,'linestyle','--','color','blue')
    title('Validation')
    ylabel('$\dot{x_4}$','fontsize',12,'interpreter','latex')
    xlabel('Time (days)','fontsize',13,'interpreter','latex')
    legend('True Prediction','Estimated Prediction')

toc

end

%% Test
fprintf('\n Simulation finished, plotting the result...... \n')

dScore_X_Noise = [dScore_X1_Noise;dScore_X2_Noise;dScore_X3_Noise;dScore_X4_Noise];

% % Concantenate the BEST model
% BEST_MODEL=[X1_BEST;X2_BEST;X3_BEST;X4_BEST];