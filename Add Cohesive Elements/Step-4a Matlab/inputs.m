% Eralp Demir
% 19/06/2025

% Inputs to the cohesive element addition

% Input filename that is used as the input
filename = 'abq_hex';

% Text used for grain-sets
% GRAINSETNAME= 'GRAIN-';%'grain';%'GRAIN-';
GRAINSETNAME = 'GRAIN-';

% Text used for material name (just before the number)
% MATERIALNAME= 'GRAIN-';%'grain-';%'MATERIAL-GRAIN';
% MATERIALNAME = 'grain-';
% MATERIALNAME = 'grain-';
MATERIALNAME = 'MATERIAL-GRAIN';

% Modified file name with cohesive elements
newfilename = 'Job_coh';

% Material-ID
matID = 2;

% Number of SDV
noDepvar = 12;

% Integration method
% 1: GAUSS-QUADRATURE, 2: NEWTON-COTES, 3: NEWTON-COTES (ABAQUS DEFAULT)
INTMTD = 1;


% Inputs for the bilinear cohesive model
% INPUTS 
% *: REQUIRED ONLY FOR MIXED MODE; "0" OTHERWISE.
% `: REQURED FOR 2D PROBLEMS ONLY 
% PROPS(1) / MODE / FRACTURE MODE / (-): 
%   = "1" FOR OPENING MODE ONLY
%   = "2" FOR shear MODE ONLY
%   = "3" FOR MIXED MODE 
% PROPS(2) / KI / PENALTY STIFFNESS FOR OPENING MODE / (N/MM)
% PROPS(3) / KII / PENALTY STIFFNESS FOR SHEARING MODE / (N/MM)
% PROPS(4) / SI / MAXIMUM NORMAL STRESS BEFORE FAILURE / (MPA)
% *PROPS(5) / SII / MAXIMUM SHEAR STRESS BEFORE FAILURE / (MPA)
% PROPS(6) / GCI / FRACTURE TOUGHNESS FOR OPENING MODE / (J)
% *PROPS(7) / GCII / FRACTURE TOUGHNESS FOR SHEARING MODE / (J)
% *PROPS(8) / ETA / BENZEGGAGH-KENANE (B-K) EXPONENT FOR MIXED MODE / (-)
% `PROPS(9) / HEIGHT / HEIGHT OF 2D ELEMENTS / (MM) 
PROPS = [3., 1.d7, 1.d7, 100., 100., 5., 5., 2., 1.];


