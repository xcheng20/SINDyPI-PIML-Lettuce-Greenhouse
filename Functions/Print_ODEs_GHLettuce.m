%% This function will display the ODE function.

% Last Updated: 2019/04/21
% Coded By: K

function ODEs=Print_ODEs_GHLettuce(ODE,n_state,n_control,n_disturbance,disp,actual)

x_vars=sym('x',[n_state;1]);
u_vars=sym('u',[n_control;1]);
d_vars=sym('d',[n_disturbance;1]);
dx_vars=sym('dx',[n_state;1]);

%Replace the input arguments in ODE with symbolic vars then, assigned to ODEs
if isa(ODE,'function_handle')
    if n_control~=0
        ODEs=ODE(0,x_vars,u_vars,d_vars); %This is just assigning the parameter value and symbolic function to ODE
    else
        ODEs=ODE(0,x_vars,d_vars);
    end
end

if disp==1
    if actual==1
        fprintf('\v The actual ODE of the system is/are :\n')
    else
        fprintf('\v The discovered ODE of the system is/are :\n')
    end
     for i=1:n_state
          digits(4)
          fprintf(strcat('\t',char(dx_vars(i,1)),'=',char(vpa((ODEs(i,1)))),'\n')); %concatenate 
     end
end