function Neper2Abaqus(filename,matID,noDepvar)





%% Input
% filename - name of the material parameter file 

% Set material ID:
% - enter "0" to use Dream3D number output
% - enter "1-10" material ID number if user defined!



% Number of state variables
% Number of outputs
% noDepvar


% The name of the excel file:
% inputfile_info.xlsx
    
% --------------------------
% Convert the Dream3D outputs to Abaqus input file
% Designed for input structure of DBF_code - A UMAT subroutine for crystal
% plasticity

% Feb. 2nd, 2022
% written by
% Eralp Demir
% eralp.demir@eng.ox.ac.uk
% --------------------------
    


% .msh file format

% $Nodes
% <number_of_nodes>
% <node_id> <node_x> <node_y> <node_z>
% ...
% $EndNodes

% $Elements
% <number_of_elements>
% <elt_id> <elt_type> <number_of_tags> <tag1> ... <elt_id_node1> ...
% ...
% $EndElements


% $NSets
% <number_of_nsets>
% <nset1_label>
% <nset_node_nb>
% <nset_node1>
% <nset_node2>
% ...
% <nset2_label>
% ...
% $EndNSets


% $PhysicalNames
% <number_of_physical_names>
% <physical_dimension> <physical_id> <physical_name>
% ...
% $EndPhysicalNames


% $ElsetOrientations
% <number_of_sets> <orientation_descriptor>
% <element_id> ori_des1> ...
% ...
% $EndElsetOrientations



tic  

format shortg



% Read  ".msh" file
fid = fopen([filename '.msh'],'r+');

% Read line by line
tline = fgetl(fid);
tlines = cell(0,1);
while ischar(tline)
    tlines{end+1,1} = tline;
    tline = fgetl(fid);
end

fclose(fid);


% Convert cell to string arrays
slines=string(tlines);


% % % Read .INP file
% % gid = fopen([filename '.inp'],'r+');
% % 
% % % Read line by line
% % ttline = fgetl(gid);
% % ttlines = cell(0,1);
% % while ischar(ttline)
% %     ttlines{end+1,1} = ttline;
% %     ttline = fgetl(gid);
% % end
% % 
% % fclose(gid);
% % 
% % 
% % % Convert cell to string arrays
% % sslines=string(ttlines);
% % 
% % 
% % %%%%%%%%%%%% PROCESSING .INP FILE %%%%%%%%%%%%%%%%%%%%%%%
% % 
% % % Find the number of elemnent from .INP file
% % for i=1:size(sslines,1)
% %     % Find the first line that start with "*Elset"
% %     idx = strfind(sslines(i),'*Elset');
% %     if idx == 1
% %         nd=i;
% %         break
% %     end
% % end
% % 
% % % Go back two lines and convert to number
% % aa = str2num(sslines(nd-2));
% % 
% % % The first element is the number of elements
% % tot_els = aa(1);
% % 
% % 
% % % Find the line number right before "*PART" from .INP file
% % for i=1:size(sslines,1)
% %     % Find the first line that start with "*END PART"
% %     idx = strfind(sslines(i),'*End Part');
% %     if idx == 1
% %         nd=i;
% %         break
% %     end
% % end
% % 
% % nd = nd - 1;
% % 
% % % The variable 'nd' will be used later when creating the "filename_.INP" file


%%%%%%%%%%%% PROCESSING .MSH FILE %%%%%%%%%%%%%%%%%%%%%%%

% Find the header line
st =find(slines=='$Nodes');

% Find the cell that starts with the total node number
tot_nodes = str2num(slines(st+1));

% Read node coordinates
crds = zeros(tot_nodes,4);
for i = st+2 : 1 : st+2+tot_nodes-1

    dummy = str2num(slines(i));
    
    % Node index
    ii = dummy(1);
    
    % Coordinates
    crds(ii,1:4) = dummy(1:4);
end





% Find the header line
st =find(slines=='$EndElements');

% Assuming the same element type throughout the mesh!!!
% Read the first line
dummy = str2num(slines(st-1));

% Number of tags
ntag = dummy(3);


neltyp = dummy(2);

% Element type (Neper documentation)
% <elt_type> is an integer specifying the type of elements: 15 for a 0D element,
% 1 for a 1st-order 1D element (2 nodes), 8 for a 2nd-order 1D element (3 nodes),
% 2 for a 1st-order triangular element (3 nodes), 
% 3 for a 1st-order quadrangular element (4 nodes), 
% 9 for a 2nd-order triangular element (6 nodes), 
% 16 for a 2nd-order quadrangular element (8 nodes), 
% 10 for a 2nd-order quadrangular element (9 nodes), 
% 4 for a 1st-order tetrahedral element (4 nodes),
% 5 for a 1st-order hexahedral element (8 nodes), 
% 11 for a 2nd-order tetrahedral element (10 nodes), 
% 17 for a 2nd-order hexahedral element (20 nodes), 
% 6 for a 1st-order prismatic element (6 nodes), 
% 18 for a 2nd-order prismatic element (15 nodes).

% If 2D problem
% PS: plane stress
% PE: plane strain
ch='PS'; % plane stress by default

switch neltyp
    
    case 2
        
        eltyp = ['C',ch,'3'];
        nnpel = 3;
        numpt = 1;
        
    case 3
        
        eltyp = ['C',ch,'4'];
        nnpel = 4;
        numpt = 4;
        
    case 9
        
        eltyp = ['C',ch,'6'];
        nnpel = 6;
        numpt = 3;
        
    case 16
        
        eltyp = ['C',ch,'8'];
        nnpel = 8;
        numpt = 4;
        
    case 4
        
        eltyp = 'C3D4';
        nnpel = 4;
        numpt = 1;
        
    case 6
        
        eltyp = 'C3D6';
        nnpel = 6;
        numpt = 2;
        
    case 5 
        
        eltyp = 'C3D8';
        nnpel = 8;
        numpt = 8;
        
    case 11
        
        eltyp = 'C3D10';
        nnpel = 10;
        numpt = 4;
        
    case 18 
        
        eltyp = 'C3D15';
        nnpel = 15;
        numpt = 9;
        
    case 17
        
        eltyp = 'C3D20';
        nnpel = 20;
        numpt = 27;
        
end






% Total number of elements
% Read from bottom line until the  <eltyp changes>
i=0; etyp = neltyp;
while etyp == neltyp 
    
    i = i + 1;
    dummy = str2num(slines(st-i));
    
    if size(dummy,2)>1
        etyp = dummy(2);
    % Reached the top of the line
    else
        break
    end
end
tot_els  = i-1;


% Read connectivity
conn = zeros(tot_els,nnpel+1);
% Grain ids
grains = zeros(tot_els,1);
% Phase ids
phases = zeros(tot_els,1);

% Material-ID
materials = ones(tot_els,1)*matID;




 ii=tot_els;
for i = st-1 : -1 : st-1-tot_els+1

    dummy = str2num(slines(i));
    

    
    % Connectivity
    conn(ii,1:nnpel+1) = [dummy(1), dummy(3+ntag+1 : 1 : 3+ntag+nnpel)];
    
    
    grains(ii) = dummy(4);
    
    % Last tag is assumed to be the phase id
    % It is normally zero so set to some number
    phases(ii) = dummy(3+ntag) ;
    
    % Element index
    ii = ii -1;
    
end






% Read set names
% Set name
% Find the header line
% Start from the bottom
st =find(slines=='$EndPhysicalNames');
cc = char(slines(st-1));
[aa,~,~,nextindex] = sscanf(cc,'%d %d');
ee=char(aa(2));
dd = cc(nextindex:end);
setname = dd(1:end-size(ee,2));












% Euler angles for each geain
% Find the header line
st =find(slines=='$ElsetOrientations');
tot_grains = sscanf(slines(st+1),'%d');
eulers = zeros(tot_grains,3);
for i=st+2:1:tot_grains+st+1
    
    dummy=str2num(slines(i));
    
    eulers(dummy(1),1:3) = dummy(2:4);
    

end


% Euler angles for each element
% Assign Euler angles to each element
euler = zeros(tot_els,3);
for i=1:1:tot_els

    % Grain ID
    grnid = grains(i);
    
    % Euler angle
    euler(i,1:3) = eulers(grnid,1:3);
    
end














[grain_order, grain_record]=unique(grains);

grain_order=grain_order';

grain_record=grain_record';



phase_order=phases(grain_record);

material_order = materials(grain_record);

euler_angle1=euler(grain_record,1);
euler_angle2=euler(grain_record,2);
euler_angle3=euler(grain_record,3);
%     



%%% Generate the rotation matrix        
% for ii=1:length(grain_order)
%     %%
%     %need to create the rotation matrix for each euler angle for each
%     %grain.  This matrix is then used to rotate a global orientation to
%     %the what is is n the local orientation.
%     zrot=[cosd(euler_angle1(ii)), sind(euler_angle1(ii)), 0; -sind(euler_angle1(ii)), cosd(euler_angle1(ii)),0; 0,0,1];
%     xrot=[1,0,0;0,cosd(euler_angle2(ii)),sind(euler_angle2(ii));0,-sind(euler_angle2(ii)),cosd(euler_angle2(ii))];
%     zrot2=[cosd(euler_angle3(ii)),sind(euler_angle3(ii)),0;-sind(euler_angle3(ii)),cosd(euler_angle3(ii)),0;0,0,1];
% 
%     %total rotation matrix - crystal to sample transformation
%     total_rot=transpose(zrot2*xrot*zrot);
% 
% 
% 
% end




% Write the overall element and node sets to input file
% open inp file and write keywords 
inpFile = fopen([filename '.inp'],'wt');
fprintf(inpFile,'** Generated by Neper and modified by: Neper2Abaqus.m\n');
fprintf(inpFile,'**PARTS\n**\n');
fprintf(inpFile,'*Part, name=NEPER\n');

% write nodes
fprintf(inpFile,'*NODE\n');
fprintf(inpFile,'%d,\t%e,\t%e, \t%e\n',crds');

% write elements
fprintf(inpFile,['*Element, type=',eltyp,'\n']);
str=[];
for i=1:1:nnpel
    str = [str, '%8d,'];
end
str = [str, '%8d\n'];
fprintf(inpFile,str,conn');


%% Write the elements sets for each grain to the input file
% create element sets containing grains
for ii = 1:numel(unique(grains))
    %%
    fprintf(inpFile,'\n*Elset, elset=GRAIN-%d\n',grain_order(ii));
    fprintf(inpFile,'%d, %d, %d, %d, %d, %d, %d, %d, %d\n',conn(grains==grain_order(ii))');
    numels=0;

    for tt=1:length(conn(grains==grain_order(ii)))
        %%
        numels=numels+1;
    end
   numels_total(grain_order(ii))=numels;
end

%% Write element set for each phase to input file
uniPhases = unique(phases);
for ii = 1:numel(unique(phases))
    fprintf(inpFile,'\n*Elset, elset=Phase-%d\n',ii);
    fprintf(inpFile,'%d, %d, %d, %d, %d, %d, %d, %d, %d\n',conn(phases==uniPhases(ii))');
end

% %% Calculate grain spherical equivalent diameter
% % calulate diamater in microns
% % additionally, the dimaters for each ground are written to a separate
% % text file to be used to developed a grain size histogram
% diameterID=fopen('diameter.txt','w');
% for ii=1:numel(unique(grains))
%     %%
%     diameter(grain_order(ii))=((((6.0/pi)*(numels_total(grain_order(ii))))^(1/3)));
%     fprintf(diameterID, '%d\n', diameter(grain_order(ii)));
% end
% fclose(diameterID);

%% write sections to each grain
for ii=1:length(grain_order)
    %%
    fprintf(inpFile,'\n**Section: Section_Grain-%d\n*Solid Section, elset=GRAIN-%d, material=MATERIAL-GRAIN%d\n,\n',grain_order(ii),grain_order(ii),grain_order(ii));
end
%% Continue writing the input file with assembly information
% write a closing keyword
fprintf(inpFile,'*End Part');

%writing assembly
fprintf(inpFile,'\n**\n**ASSEMBLY\n**');
fprintf(inpFile,'\n*Assembly, name=Assembly\n**');
fprintf(inpFile,'\n*Instance, name=NEPER-1, part=NEPER\n');

% write nodes
fprintf(inpFile,'*NODE\n');
fprintf(inpFile,'%d,\t%e,\t%e, \t%e\n',crds');

% write elements
fprintf(inpFile,['*Element, type=', eltyp, '\n']);
str=[];
for i=1:1:nnpel
    str = [str, '%8d,'];
end
str = [str, '%8d\n'];
fprintf(inpFile,str,conn');

fprintf(inpFile,'\n*End Instance\n**');


%% Closing the assembly component of the input file

fprintf(inpFile,'\n*End Assembly');


fprintf(inpFile, '\n**MATERIALS\n**');

%import material parameters to be used in the development of materials
%for each grain.
xlRange='A1:A6';
[A16]=readmatrix('PROPS.xlsx','Sheet','Material_parameters','DataRange',xlRange);


%% Finalising the input file

% Flag for reading the inputs from the file or material library
% "0": material library in usermaterial.f will be used
% "1": use the material parameters in excel file



% Are the material properties given in the excel file?
% Read from excel file (if read_all_props==true)
if A16(6)==0
    
    
    A = strings(6,1);

    
    % Flag for reading the inputs from the file or material library
    % "0": material library in usermaterial.f will be used
    % "1": use the material parameters in excel file
    A(6)=0;
    
    % Number variables in DEPVAR
    noPROPS = 6;
    
    % Do not define anything further
    
else
    
    A = strings(250,1);


        
    % Flag for reading the inputs from the file or material library
    % "0": material library in usermaterial.f will be used
    % "1": use the material parameters in excel file
    A(6)=1;   
    
    
    % Number variables in DEPVAR
    % Has a fixed size - including additional space for extra variables
    noPROPS = 250;
    
    
    
end








for ii=1:length(grain_order)

    fprintf(inpFile, '\n*Material, name=MATERIAL-GRAIN%d',grain_order(ii));
    fprintf(inpFile, ['\n*Depvar\n', num2str(noDepvar), ',']);
    fprintf(inpFile, ['\n*User Material, constants=',num2str(noPROPS),'\n']);

    % Euler angles
    A(1:3) = [euler_angle1(ii), euler_angle2(ii), euler_angle3(ii)];
    % Grain - ID
    A(4) = grain_order(ii);
    
    
    % Phase - ID
    % IF DEFINED BY THE USER (>0)
    if matID>0
        
        A(5) = material_order(ii);
        
    % Use Dream3D output 
    else
        
        A(5) = phase_order(ii);
        
    end

%     % Adding the centroid information in x,y,z coordinates
%     A(9:11)=centroid(grain_order(ii),:);
%     % center element

%     %adding the calculated equivalent spherical diameter for each grain
%     A(13)=diameter(grain_order(ii));

    % Read the properties from the PROPS if desired
    if A16(6) ==1
        
        % Loop through all different phases
        for iph = 1:numel(unique(phases))
            
            % Column character (read next column for each phase)
            letter = char(iph+ 64);

            xlRange = [letter,'1:', letter, '250']; % A1-A250
            %[B]=xlsread('inputfile_info.xlsx','Material_parameters',xlRange);
            [B]=readmatrix('PROPS.xlsx','Sheet','Material_parameters','DataRange',xlRange);

            

            A(7:250) = B(7:250);
            
            
            
        end
        
    end





    % Printing this information to file
    fprintf(inpFile, '%s, %s, %s, %s, %s, %s, %s, %s\n',A);

end

fprintf(inpFile,'\n**');
fprintf(inpFile, '\n**\n** STEP: Loading\n**\n*Step, name=Loading, nlgeom=YES, inc=10000\n*Static\n0.01, 10., 1e-05, 1.');
fprintf(inpFile, '\n**\n** OUTPUT REQUESTS\n**');
fprintf(inpFile, '\n*Restart, write, frequency=0\n**');
fprintf(inpFile, '\n** FIELD OUTPUT: F-Output-1\n**\n*Output, field, variable=PRESELECT\n**');
if noDepvar>0
    fprintf(inpFile, '\n** FIELD OUTPUT: F-Output-2\n**\n*Element Output, directions=YES\nSDV,\n**');
end
fprintf(inpFile, '\n** HISTORY OUTPUT: H-Output-1\n**\n*Output, history, variable=PRESELECT\n**');
fprintf(inpFile, '\n*End Step');

% close the file
fclose(inpFile);


toc

return

end




