function p = DefineParameters
% Model parameters
% parameter 					description 									[unit] 					nominal value
p.satH2O1 = 9348; 				% saturation water vapour parameter 			[J m^{-3}] 				9348
p.satH2O2 = 17.4; 				% saturation water vapour parameter 			[-] 					17.4
p.satH2O3 = 239; 				% saturation water vapour parameter 			[°C] 					239
p.satH2O4 = 10998;  			% saturation water vapour parameter 			[J m^{-3}] 				10998
p.R = 8314; 					% ideal gas constant 							[J K^{-1} kmol^{-1}] 	8314
p.T = 273.15; 					% conversion from C to K 						[K] 					273.15

p.leak = 0.75e-4; 				% ventilation leakage through the cover 		[m s^{-1}] 				0.75e-4
p.CO2cap = 4.1; 				% CO2 capacity of the greenhouse 				[m^3{air} m^{-2}{gh}]   4.1
p.H2Ocap = 4.1; 				% Vapor capacity of the greenhouse 				[m^3{air} m^{-2}{gh}]   4.1
p.aCap = 3e4; 					% effective heat capacity of the greenhouse air [J m^{-2}{gh} °C^{-1}]  3e4
p.ventCap = 1290; 				% heat capacity per volume of greenhouse air 	[J m^{-3}{gh} °C^{-1}]  1290
p.trans_g_o = 6.1; 				% overall heat transfer through the cover 		[W m^{-2}{gh} °C^{-1}]  6.1
p.rad_o_g = 0.2; 				% heat load coefficient due to solar radiation 	[-] 					0.2

p.alfaBeta = 0.544; 			% yield factor 									[-] 					0.544
p.Wc_a = 2.65e-7; 				% respiration rate 								[s^{-1}] 				2.65e-7
p.CO2c_a = 4.87e-7; 			% respiration coefficient 						[s^{-1}]  				4.87e-7
p.laiW = 53; 					% effective canopy surface 						[m^2{leaf} kg^{-1}{dw}] 53
p.photI0 = 3.55e-9; 			% light use efficiency 							[kg{CO2} J^{-1}]  		3.55e-9
p.photCO2_1=5.11e-6;  			% temperature influence on photosynthesis 		[m s^{-1} °C^{-2}] 		5.11e-6
p.photCO2_2=2.3e-4;				% temperature influence on photosynthesis 		[m s^{-1} °C^{-1}] 		2.3e-4
p.photCO2_3=6.29e-4; 			% temperature influence on photosynthesis 		[m s^{-1}] 				6.29e-4
p.photGamma = 5.2e-5; 			% carbon dioxide compensation point 			[kg{CO2} m^{-3}{air}] 	5.2e-5
p.evap_c_a = 3.6e-3; 			% coefficient of leaf-air vapor flow 			[m s^{-1}] 				3.6e-3
p.energyCost = 6.35e-9/2.20371; % price of energy                               [€ J^{-1}]              6.35e-9 [Dfl J^{-1}] (division by 2.20371 represents currency conversion)
p.co2Cost = 42e-2/2.20371;      % price of CO2                                  [€ kg^{-1}{CO2}]        42e-2 [Dfl kg^{-1}{CO2}] (division by 2.20371 represents currency conversion)
p.productPrice1 = 1.8/2.20371;  % parameter for price of product                [€ m^{-2}{gh}]          1.8 [Dfl kg^{-1}{gh}] (division by 2.20371 represents currency conversion)
p.productPrice2 = 16/2.20371;   % parameter for price of product                [€ kg^{-1}{gh} m^{-2}{gh}] 16 (division by 2.20371 represents currency conversion)

% Perturbation parameters
% parameter 					description 									[unit] 					nominal value
p.c.photo = 1;                  % perturbation of gross photosynthesis          [-]                     1
p.c.grossPhot = 1;              % perturbation of net photosynthesis (bad name!)[-]                     1
p.c.resp = 1;                   % perturbation of respiration rate              [-]                     1
p.c.transp = 1;                 % perturbation of transpiration rate            [-]                     1
p.c.rad = 1;                    % perturbation of outdoor global radiation      [-]                     1
p.c.rhAdd = 0;                  % perturbation of outdoor relative humidity     [%]                     0
p.c.tempAdd = 0;                % perturbation of outdoor temperature           [°C]                    0
p.c.co2 = 1;                    % perturbation of outdoor CO2 concentration     [-]                     1


p.lue = 7.5e-8;
p.heatLoss = 1;
p.heatEff = 0.1;
p.gasPrice = 4.55e-4;
p.lettucePrice = 136.4;
p.heatMin = 0;
p.heatMax = 100;
