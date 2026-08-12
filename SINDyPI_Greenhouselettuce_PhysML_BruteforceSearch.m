% Physics-Informed ML of Greenhouse Lettuce Model Identification
% This code is the modification from SINDY-PI code
%Date: 01/04/2026
%Code by: Fakhira
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

% Close all, clear all, clc
close all;clear all; clc;
[status,msg] = mkdir('Results');
addpath('Functions')
set(0,'defaulttextInterpreter','latex')

%% ===============================================PREPARING FOR SIMULATION========================================================
%% Gather the greenhouse weather data for 3 datasets: training, testing, and
% validation
% Model parameters
ops.nx  = 4;    % #states
ops.nu  = 3;    % #controllable inputs
ops.nd  = 4;    % #disturbance inputs

% simulation parameters
c         = 86400;                          % number of seconds in a day (24*60*60)
nDays     = 40;                             % days in simulation (max=320)
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

%% Run the ODE files and gather the simulation data
% Add noise level to train/test and validation data
Noise = 0.00; %-9
Noise_val=0.00;

% Shuffle the data
Shuffle = 0;

%Perform Simulation for the Training
[dData,Data]=Get_Sim_Data_GHLettuce(@(t,x,u,d)GHLettuce_ODE(t,x,u,d),x0,u,d,ops.t,ops.nu,ops.nd,Noise,Shuffle);

%Perform Simulation for the Test Data
[dData_test,Data_test]=Get_Sim_Data_GHLettuce(@(t,x,u,d)GHLettuce_ODE(t,x,u,d),x0_test,u_test,d_test,ops.t,ops.nu,ops.nd,Noise,Shuffle);

%Perform Simulation for the Validation Data
[dData_val,Data_val]=Get_Sim_Data_GHLettuce(@(t,x,u,d)GHLettuce_ODE(t,x,u,d),x0_val,u_val,d_val,ops.t,ops.nu,ops.nd,Noise_val,Shuffle);

% Plot the Weather Data
figure(5)
PlotResults_GHLettuce_train
figure(6)
PlotResults_GHLettuce_test
figure(7)
PlotResults_GHLettuce_val

%% Sparse Regression Preparation
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

% Print the actual ODE we try to discover 
digits(4)
Print_ODEs = Print_ODEs_GHLettuce(@(t,x,u,d)GHLettuce_ODE(t,x,u,d),n_state,n_control,n_disturbance,disp_actual_ode,actual);

% Create symbolic states
dx=sym('dx',[n_state,1]);

% Theta Library Preparation
lam=[1e-20;1e-19;1e-18;1e-17;1e-16;1e-15;1e-14;1e-13;1e-12]; %set many lambda values
N_iter=20;  
disp=1;
NormalizeLib=1; 
tic

%% =============================================== SPARSE REGRESSION USING BRUTE-FORCE SEARCH OF CANDIDATE LIBRARIES ========================================================
% Define the LHS and RHS libraries for each states
% 1-4 represents the polynomial order. Below is the default library settings, you
% can play changing the library design and see how it affects the
% identification process

%LHS Guess
dXPoly_Order_Guess=[2;2;2;2]; 
%RHS Library
Poly_Order=[1;1;1;1]; 
dXPoly_Order=[1;2;3;4]; 

    %Sparse Regression of X1 to X4
    for iter=1:n_state
        fprintf('\n \n Calculating the %i expression...\n',iter)
        % Change the library base on what equation you want to work on
        Highest_dXPoly_Order_Guess= dXPoly_Order_Guess(iter);
        Highest_dXPoly_Order= dXPoly_Order(iter); 
        Highest_Poly_Order = Poly_Order(iter); 

        % According to the previous parameter generate the left hand side guess
        [LHS_Data,LHS_Sym]=GuessLib_GHLettuce_PhysML_Test(Data,dData,u,d,iter,Highest_dXPoly_Order_Guess);  
        %Store the result
        LHS_Sym_all{:,iter}=LHS_Sym; 

        %Generate RHS Library
        [SINDy_Data,SINDy_Struct]= SINDyLib_GHLettuce_PhysML_Best(Data,dData,u,d,iter,Highest_Poly_Order, Highest_dXPoly_Order);

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
                %If the previous ODE is 0, set the score as NaN, else
                %calculate it
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
   
    
%% =============================================== BEST IDENTIFIED MODEL GENERATION ========================================================
% Generate the best model from previous search   
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
    [dData_True_val,Data_True_val]=Get_Sim_Data_GHLettuce(@(t,x,u,d)GHLettuce_ODE(t,x,u,d),x0_val,u_val,d_val,ops.t,ops.nu,ops.nd,Noise_val,Shuffle);
    [dData_Es_val,Data_Es_val]=Get_Sim_Data_GHLettuce(@(t,x,u,d)Sindy_ODE_RHS_GHLettuce(t,x,u,d),x0_val,u_val,d_val,ops.t,ops.nu,ops.nd,Noise_val,Shuffle);
    
    % Prediction Error 
    dScore_1_val_Best=(norm(dData_True_val(:,1)-dData_Es_val(:,1)))/norm(dData_True_val(:,1));
    dScore_2_val_Best=(norm(dData_True_val(:,2)-dData_Es_val(:,2)))/norm(dData_True_val(:,2));
    dScore_3_val_Best=(norm(dData_True_val(:,3)-dData_Es_val(:,3)))/norm(dData_True_val(:,3));
    dScore_4_val_Best=(norm(dData_True_val(:,4)-dData_Es_val(:,4)))/norm(dData_True_val(:,4));
    
    % Store the result of Simulation
    dScore_all_Best=[dScore_1_val_Best dScore_2_val_Best dScore_3_val_Best dScore_4_val_Best];
    
    % Store individual result
    %STATE 1
    Lambda_X1_Best= lam(minIndex1(1,minIndex2(1,1)));
    LHS_X1_Best= LHS_Sym_all{1,1}{1,minIndex2(1,1)};
    dx1_1 = {'dx1*x1','dx1*x1^2'};
    Sym_Struct_X1_Best= strjoin(dx1_1);
    dScore_X1_Best = dScore_1_val_Best;
    ODEs_X1_Best = ODE_Best(1,1);
    X1_BEST = [string(LHS_X1_Best) Sym_Struct_X1_Best Lambda_X1_Best string(dScore_X1_Best) string(ODEs_X1_Best)];
    
    %STATE 2
    Lambda_X2_Best= lam(minIndex1(2,minIndex2(2,1)));
    LHS_X2_Best= LHS_Sym_all{1,2}{1,minIndex2(2,1)};
    Sym_Struct_X2_Best= SINDy_Struct_all{1,2}(1,end);
    dScore_X2_Best = dScore_2_val_Best;
    ODEs_X2_Best = ODE_Best(2,1);
    X2_BEST = [string(LHS_X2_Best) Sym_Struct_X2_Best Lambda_X2_Best string(dScore_X2_Best) string(ODEs_X2_Best)];
    
    %STATE 3
    Lambda_X3_Best= lam(minIndex1(3,minIndex2(3,1)));
    LHS_X3_Best= LHS_Sym_all{1,3}{1,minIndex2(3,1)};
    Sym_Struct_X3_Best= SINDy_Struct_all{1,3}(1,end);
    dScore_X3_Best = dScore_3_val_Best;
    ODEs_X3_Best = ODE_Best(3,1);
    X3_BEST = [string(LHS_X3_Best) Sym_Struct_X3_Best Lambda_X3_Best string(dScore_X3_Best) string(ODEs_X3_Best)];
    
    %STATE 4
    Lambda_X4_Best= lam(minIndex1(4,minIndex2(4,1)));
    LHS_X4_Best= LHS_Sym_all{1,4}{1,minIndex2(4,1)};
    Sym_Struct_X4_Best= SINDy_Struct_all{1,4}(1,end);
    dScore_X4_Best = dScore_4_val_Best;
    ODEs_X4_Best = ODE_Best(4,1);
    X4_BEST = [string(LHS_X4_Best) Sym_Struct_X4_Best Lambda_X4_Best string(dScore_X4_Best) string(ODEs_X4_Best)];
    
    % Concantenate the BEST model
    BEST_MODEL=[X1_BEST;X2_BEST;X3_BEST;X4_BEST];

    % Show process
    fprintf('\n Simulation finished! \n')
    close all
    % Create the new directory to save the plot
    % [fld_status, fld_msg, fld_msgID]=mkdir('Figures');
    
    %% Plot the result of simulation
    fprintf('\n Now, plotting the result...... \n')
    figure(8)
    plot(ops.t/c,Data_True_val(:,1),'linewidth',0.75,'Color','green')
    hold on
    plot(ops.t/c,Data_Es_val(:,1),'linewidth',0.75,'linestyle','--','color','blue')
    title('Validation')
    ylabel('$x_1$ (g/m$^2$)','fontsize',16,'interpreter','latex')
    xlabel('Time (days)','fontsize',14,'interpreter','latex')
    legend('Nominal Prediction','Estimated Prediction')
    
    figure(9)
    plot(ops.t/c,dData_True_val(:,1),'linewidth',0.6,'Color','green')
    hold on
    plot(ops.t/c,dData_Es_val(:,1),'linewidth',0.6,'linestyle','--','color','blue')
    title('Validation')
    ylabel('$\dot{x_1}$','fontsize',16,'interpreter','latex')
    xlabel('Time (days)','fontsize',14,'interpreter','latex')
    legend('Nominal Prediction','Estimated Prediction')
   
    figure(10)
    plot(ops.t/c,Data_True_val(:,2),'linewidth',0.5,'Color','black')
    hold on
    plot(ops.t/c,Data_Es_val(:,2),'linewidth',0.5,'linestyle','--','color','blue')
    title('Validation')
    ylabel('$x_2$ (kg/m$^3$)','fontsize',12,'interpreter','latex')
    xlabel('Time (days)','fontsize',13,'interpreter','latex')
    legend('True Prediction','Estimated Prediction')
    
    figure(11)
    plot(ops.t/c,dData_True_val(:,2),'linewidth',0.5,'Color','black')
    hold on
    plot(ops.t/c,dData_Es_val(:,2),'linewidth',0.5,'linestyle','--','color','blue')
    title('Validation')
    ylabel('$\dot{x_2}$','fontsize',12,'interpreter','latex')
    xlabel('Time (days)','fontsize',13,'interpreter','latex')
    legend('True Prediction','Estimated Prediction')
    
    figure(12)
    plot(ops.t/c,Data_True_val(:,3),'linewidth',0.5,'Color','black')
    hold on
    plot(ops.t/c,Data_Es_val(:,3),'linewidth',0.5,'linestyle','--','color','blue')
    title('Validation')
    ylabel('$x_3$ ($^\circ$C)','fontsize',12,'interpreter','latex')
    xlabel('Time (days)','fontsize',13,'interpreter','latex')
    legend('True Prediction','Estimated Predictio')
    
    figure(13)
    plot(ops.t/c,dData_True_val(:,3),'linewidth',0.5,'Color','black')
    hold on
    plot(ops.t/c,dData_Es_val(:,3),'linewidth',0.5,'linestyle','--','color','blue')
    title('Validation')
    ylabel('$\dot{x_3}$','fontsize',12,'interpreter','latex')
    xlabel('Time (days)','fontsize',13,'interpreter','latex')
    legend('True Prediction','Estimated Prediction')
    
    figure(14)
    plot(ops.t/c,Data_True_val(:,4),'linewidth',0.5,'Color','black')
    hold on
    plot(ops.t/c,Data_Es_val(:,4),'linewidth',0.5,'linestyle','--','color','blue')
    title('Validation')
    ylabel('$x_4$ (kg/m$^3$)','fontsize',12,'interpreter','latex')
    xlabel('Time (days)','fontsize',13,'interpreter','latex')
    legend('True Prediction','Estimated Prediction')
    
    figure(15)
    plot(ops.t/c,dData_True_val(:,4),'linewidth',0.5,'Color','black')
    hold on
    plot(ops.t/c,dData_Es_val(:,4),'linewidth',0.5,'linestyle','--','color','blue')
    title('Validation')
    ylabel('$\dot{x_4}$','fontsize',12,'interpreter','latex')
    xlabel('Time (days)','fontsize',13,'interpreter','latex')
    legend('True Prediction','Estimated Prediction')

    sprintf('\n Finished Simulating The Best Model...') 
toc

%% =============================================== PLOTTING LAMBDA SCORE FOR EACH LHS GUESS ========================================================
    for iter=1:n_state
        %iter=3; %activate to only evaluate x3 
    % Plot all lambdas value
    figure(iter)
        for i=1:length(LHS_Sym_all{1,iter})
            for j=1:length(lam)
                plot(lam(j,1),Score(iter,i,j),'color',[1 iter*0.2 iter*0.25],'linewidth',1.5)
                s1 = scatter(lam(j,1),Score(iter,i,j),70,'MarkerFaceColor',[1 iter*0.2 iter*0.25],'MarkerEdgeColor',[0 0 0],'linewidth',1.5);
                set(gca,'XScale','log')
                ylim([-2 60])
                hold on
                grid on
            end
        end
        %Only Plot the best lambda and add legend
        if iter ==1 || iter ==2 || iter ==4
            hold on
            grid on
            plot(lam(1,1),Score(1,1,1),'color',[1 iter*0.2 iter*0.25],'linewidth',1.5);
            s2 = scatter(lam(1,1),Score(1,1,1),70,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 0],'linewidth',1.5,'DisplayName','Correct LHS Guess');
            ylabel('Prediction Error','fontsize',20,'interpreter','latex');
            xlabel('Sparsity Hyperparameter $\lambda$','fontsize',20,'interpreter','latex');
        title(sprintf('Model Selection of $x_%d$',iter),'fontsize',20)
        legend([s1,s2],'Incorrect LHS Guess', 'Correct LHS Guess');
        elseif iter ==3
            plot(lam(5,1),Score(1,1,5),'color',[1 iter*0.2 iter*0.25],'linewidth',1.5);
            s2 = scatter(lam(5,1),Score(1,1,5),70,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 0],'linewidth',1.5,'DisplayName','Correct LHS Guess');
            ylabel('Prediction Error','fontsize',20,'interpreter','latex');
            xlabel('Sparsity Hyperparameter $\lambda$','fontsize',20,'interpreter','latex');
        title(sprintf('Model Selection of $x_%d$',iter),'fontsize',20)
        legend([s1,s2],'Incorrect LHS Guess', 'Correct LHS Guess');
        end 
    end

%% Bar graph
xlab = {'x_1','x_2','x_3','x_4'};
vals = [18 10;
        37 18;
        6  6;
        16 12];
figure
b = bar(1:size(vals,1), vals);
xticks(1:size(vals,1))
xticklabels(xlab)

Leg = {'Nominal Model','Identified Model'};
legend(Leg, 'Location', 'best')   % or 'northwest'

ylabel('Number of terms (monomials)', 'FontSize', 14, 'Interpreter', 'latex')
xlabel('ODE Model', 'FontSize', 14, 'Interpreter', 'latex')
set(gca, 'FontSize', 14)

xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1, ytips1, labels1, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom')

xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = string(b(2).YData);
text(xtips2, ytips2, labels2, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom')