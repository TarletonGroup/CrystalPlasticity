% Eralp Demir
% 12/08/2024
% 14/04/2025 updated
% 22/04/2025 updated for 3D
%
function addcohesive
%
% ***************************************************
% INPUTS
%%
inputs
%
% Available element types
% 2D
% linear elements: CPS3, CPE3, CPS4, CPE4, CPS4R, or CPE4R, 
% quadratic elements: CPS6, CPE6, CPS6R, CPE6R, CPS8, CPE8, CPS8R, or CPE8R.
% 3D
% linear triangular prism: C3D4,
% linear hexahedral: C3D8, C3D8R
% quadratic triangular prism: C3D10
% quadratic hexahedral: C3D20, C3D20R.
%
%
% NOTE: Cohesive element type will be set according to the selected solid 
% element type.
%
%
% % Material-ID
% matID=2;
% %
% % Number of SDV
% noDepvar=1;
% %
% % Filename
% filename='Job-1';
% %
% %
% % New filename
% newfilename='testcase_coh';
% %
% %
% % Text used for grain-sets
% GRAINSETNAME= 'GRAIN-';%'grain';%'GRAIN-';
% %
% % Text used for material name (just before the number)
% MATERIALNAME= 'GRAIN-';%'grain-';%'MATERIAL-GRAIN';
% %
% % Integration method
% % 1: GAUSS-QUADRATURE, 2: NEWTON-COTES, 3: NEWTON-COTES (ABAQUS DEFAULT)
% INTMTD=1;
%
% ***************************************************
%
%
%
% Calculations start HERE!
%%
%
tic
%


% Read the INP file
%%
fid = fopen([filename '.inp'],'r+');
%
% Read line by line
tline = fgetl(fid);
tlines = cell(0,1);
while ischar(tline)
    tlines{end+1,1} = tline;
    tline = fgetl(fid);
end
%
fclose(fid);
%
%
% Convert cell to string arrays
slines=string(tlines);
%
%
% Find the header line - node coordinates
st =find(slines=='*Node');
if isempty(st)
    st =find(slines=='*NODE');
end
if isempty(st)
    disp('Could not locate the node list!')
    return
end
%
%
%
%
% Connectivity
crds=[];i=st(1);
cmp=false;
while ~cmp
    i=i+1;
    crds=[crds; str2num(slines(i))];
    cmp=strncmpi(slines(i),'*Element, type=',15);
end
%
%
elchar=char(slines(i));
elchar=elchar(16:end);
%
switch elchar
%
%
    case {'CPS3', 'CPE3'}
        %
        eltyp=3*2; % 6
        nnpel = 3;
        %
        nsurf=3;
        nnps=2;
        surfs=[ 1, 2;
                2, 3;
                3, 1];
        % Number of integration points of the cohesive element 
        numpt=2; % INDEPENDENT OF INTMTD
        %
    case {'CPS4', 'CPE4', 'CPS4R', 'CPE4R'}
        %
        eltyp=4*2; % 8
        nnpel = 4;
        %
        nsurf=4;
        nnps=2;
        %
        surfs=[1, 2;
               2, 3;
               3, 4;
               4, 1];
        % Number of integration points of the cohesive element 
        numpt=2; % INDEPENDENT OF INTMTD
        %
    case {'CPS6', 'CPE6','CPS6R', 'CPE6R'}
        %
        eltyp=6*2; % 12
        nnpel = 6;
        nsurf=3;
        nnps=3;
        surfs=[ 1, 2, 4;
                2, 3, 5;
                3, 1, 6];
        % Number of integration points of the cohesive element 
        if INTMTD==1
            numpt=3;
        elseif INTMTD==2
            numpt=3;
        elseif INTMTD==3
            numpt=2;
        end
        %
    case {'CPS8', 'CPE8', 'CPS8R', 'CPE8R'}
        %
        eltyp=8*2; % 16
        nnpel = 8;
        nsurf=4;
        nnps=3;
        surfs=[ 1, 2, 5;
                2, 3, 6;
                3, 4, 7;
                4, 1, 8];
        % Number of integration points of the cohesive element 
        if INTMTD==1
            numpt=3;
        elseif INTMTD==2
            numpt=3;
        elseif INTMTD==3
            numpt=2;
        end
        %
    case {'C3D4'}
        %
        eltyp=4*3; % 12 - same as CPS6/..
        nnpel = 4;
        %
        % Surfacea of the solid element
        nsurf=4;
        nnps=3;
        %
        surfs=[1,2,3;
               1,3,4;
               3,2,4;
               1,4,2;];
       % Number of integration points of the cohesive element 
       numpt=1; % INDEPENDENT OF INTMTD
       %
    case {'C3D8','C3D8R'}
        %
        eltyp=8*3; % 24
        nnpel = 8;
        %
        nsurf=6;
        nnps=4;
        %
        surfs=[4,3,2,1;
               5,6,7,8;
               1,5,8,4;
               2,3,7,6;
               1,2,6,5;
               3,4,8,7];
        % Number of integration points of the cohesive element
        numpt=4; % INDEPENDENT OF INTMTD
        %
     case {'C3D10'}
        %
        eltyp=10*3; % 30
        nnpel = 10;
        %
        % Surfacea of the solid element
        nsurf=4;
        nnps=6;
        %
        surfs=[1,2,3,5,6,7;
               1,3,4,7,10,8;
               3,2,4,6,9,10;
               1,4,2,8,9,5;];
       % Number of integration points of the cohesive element 
       if INTMTD==1
           numpt=4;
       elseif INTMTD==2
           numpt=6;  % corrected on 07/08/25
       elseif INTMTD==3
           numpt=3;
       end
       %
     case {'C3D20','C3D20R'}
        %
        eltyp=20*3; % 60
        nnpel = 20;
        %
        nsurf=6;
        nnps=8;
        %
        surfs=[4,3,2,1,11,10,9,12;
               5,6,7,8,13,14,15,16;
               1,5,8,4,17,16,20,12;
               2,3,7,6,10,19,14,18;
               1,2,6,5,9,18,13,17;
               3,4,8,7,11,20,15,19];
        % Number of integration points of the cohesive element
        if INTMTD==1
            numpt=9;
        elseif INTMTD==2
            numpt=9;
        elseif INTMTD==3
            numpt=4;
        end
end
%
%
%
%
%
%
% Number of nodes
tot_nodes=size(crds,1);
%
% Dimension
numdim=size(crds,2)-1;
%
% Coordinates
conn=[];
cmp=false;cmp1=false;cmp2=false;
% GoTo connecitivty - the first line
i=i+1;
%
% Decide on wheter the
% Assuming that the number of elements are greater than "1":
% A single line connectivity
if length(str2num(slines(i)))==length(str2num(slines(i+1)))
    singleline=true;
% C3D20 element the next line is used for connectivity    
else
    singleline=false;
end
%
% Read the connectivity
while ~cmp
    %
    % A single line connectivity
    if singleline
        conn=[conn; str2num(slines(i))];
        i=i+1; % goto next line
    % C3D20 element the next line is used for connectivity
    else
        aa=[str2num(slines(i)), str2num(slines(i+1))];
        conn=[conn; aa];
        i=i+2;
    end
    %
    cmp1=strncmpi(slines(i),'*Elset',6);
    cmp2=strncmpi(slines(i),'*Nset',5);
    cmp=cmp1 || cmp2;
    %
end
%
%
% Number of solid elements
tot_els = size(conn,1);
%
% 
GRAINSETCHAR=['*Elset, elset=',GRAINSETNAME];
NGSET=length(GRAINSETCHAR);
% Find all element sets
a=strncmpi(slines,GRAINSETCHAR,NGSET);
%
ind=find(a==1);
ngrain=length(ind);
grains=zeros(1,tot_els);
%
% Assuming 
for i=1:ngrain
    % character array
    ch=char(slines(ind(i)));
    % Elements are entered using generate function
    if endsWith(slines(ind(i)),', generate')
        gID=str2num(ch(NGSET+1:end-10));
        %
        k=ind(i)+1;
        %
        val=str2num(slines(k));
        %
        grains(val(1):val(3):val(2))=gID;%i;
        %
        % Elements are entered one-by-one
    else
       gID=str2num(ch(NGSET+1:end));
       j=ind(i)+1;
       while ~startsWith(slines(j),'*')
           
           grains(str2num(slines(j)))=gID;%i;
           j=j+1;
       end
       %
    end
    %
end
%
%
% Phase
phases=ones(1,tot_els);
%
%
% Material-ID
materials = ones(tot_els,1)*matID;
%
%
% Find the line with material id
a=strncmpi(slines,'*User Material,',15);
% Read the grain-id
% Euler angles per grain
ind=find(a==1);
eulers=zeros(length(ind),3);
%
MATERIALNAME=['name=', MATERIALNAME];
NMATNAME=length(MATERIALNAME);
a=strfind(slines,MATERIALNAME);
ind=[];
for i=1:length(a)
    if ~isempty(a{i})
        ind=[ind;i];
    end
end
% Find the grain numbers
% Not necessarily in order
sno=NMATNAME+11+1;
for i=1:length(ind)
    aa=slines(ind(i));
    bb=char(aa);
    gID=str2num(bb(sno:end));
    % Read the property line
    b=str2num(slines(ind(i)+4));
    eulers(gID,1:3)=b(1:3);
end
%
% Euler angles for each element
% Assign Euler angles to each element
euler = zeros(tot_els,3);
for i=1:1:tot_els
    %
    % Grain ID
    grnid = grains(i);
    %
    % Euler angle
    euler(i,1:3) = eulers(grnid,1:3);
    %
end
%
%
%
%
[grain_order, grain_record]=unique(grains);
%
grain_order=grain_order';
%
grain_record=grain_record';
%
%
%
phase_order=phases(grain_record);
%
material_order = materials(grain_record);
%
euler_angle1=euler(grain_record,1);
euler_angle2=euler(grain_record,2);
euler_angle3=euler(grain_record,3);
%
%
% Add cohesive elements
%
%
%% Cohesive modifications
%
% For each element search for its neighbors
%
% Cohesive elemetn connectivity ==> element no. to element no.
cohel_surf=[];
cohel_comp=[];
%
conn_=conn;
crds_=crds;
tot_nodes_=tot_nodes;
%
%
conn_coh=[];
%
%
% Loop through all of the nodes
for i=1:1:size(crds,1)
    %
    % Node number
    inod=crds(i,1);
    %
    % find the elements that connect to that node
    [ieles, inodes]=find(conn_(:,2:end)==inod);
    %
    % grain IDs
    gIDs=grains(ieles);
    %
    ugIDs=unique(gIDs);
    %
    % if the node is at a GB
    if length(ugIDs)>1
        %
        % Excluding the first grain
        for j=2:length(ugIDs)
            %
            % Unique grain-ID
            gID=ugIDs(j);
            %
            % Find the index in the array having the same grain-ID
            igid=find(gIDs==gID);
            %
            % Find the element numbers
            iele=ieles(igid);
            %
            % Node numbers
            inode=inodes(igid);
            %
            % Increment node number
            tot_nodes_=tot_nodes_+1;
            %
            % Change the solid connectivity
            for k=1:1:length(iele)
                %
                conn_(iele(k),inode(k)+1)=tot_nodes_;
                %
                % Add node coordinates
                crds_(tot_nodes_,1:numdim+1)=[tot_nodes_, crds(inod,2:end)];
                %
            end
            %
        end
    end
%
end
%
%
disp('Modification of the solid elements is completed!')
%
% loop through all elements
% Identify cohesive interfaces for each element
coh_int=zeros(tot_els,nsurf);
for iele0=1:tot_els
    %
    % home element connectivity
    nodes0(1,:) = conn(iele0,2:nnpel+1);
    %
    % home element center point coordinate
    xyz0=mean(crds(nodes0,2:end));
    %
    % home element neighbors
    neighs0=[];
    for i=1:nnpel
        [ele, nod]=find(nodes0(i)==conn(:,2:end));
        neighs0=[neighs0;ele];
    end
    neighs0=unique(neighs0);
    %
    % loop through home element surfaces
    for i0=1:nsurf
        %
        % home element surface nodes
        snodes0 = nodes0(surfs(i0,1:nnps));
        %
        % coordinates
        sxyz0=crds_(snodes0,2:end);
        %
        % search for the reaminig elements having the same surface connectivity
        for iele=1:length(neighs0)
            %
            % element number
            iele1=neighs0(iele);
            %
            % If the element is different than the home element
            if ~(iele1==iele0)
                %
                % If it is a grain boundary
                if ~(grains(iele0)==grains(iele1))                      
                    %
                    % neighbor element nodes
                    nodes1 = conn(iele1,2:nnpel+1);
                    %    
                    xyz1=mean(crds(nodes1,2:end));
                    %
                    % loop through home element surfaces
                    for i1=1:nsurf
                        %
                        % neighbor element surface nodes
                        snodes1 = nodes1(surfs(i1,1:nnps));
                        %
                        % mean center coordinates of the surface
                        sxyz1=crds_(snodes1,2:end);
                        %
                        % If the same surface is shared with another element
                        if norm(sort(sxyz0)-sort(sxyz1))<1e-19
                            %
                            cohel_surf=[cohel_surf; ...
                                iele0, i0, iele1, i1];
                            %
                            cohel_comp=[cohel_comp; ...
                                (iele0-1)*nsurf+i0, (iele1-1)*nsurf+i1];
                            %
                        end
                        %
                    end
                    %
                end
                %
            end
            %
        end
        %
     end
    %
    %
    disp([ num2str(iele0/tot_els*100),' % of cohesive surfaces is identified!'])
    %
end
%
% 
% sort the cohesive surfaces
cohel_comp=sort(cohel_comp,2);
%
[cohel_comp,ind]=unique(cohel_comp,'rows');
%
cohel_surf = cohel_surf(ind,:);
%
% Total number of cohesive elements
cohel = size(cohel_surf,1);
%
% Connectivity using original node numbers
for i=1:cohel
    iele0=cohel_surf(i,1);
    snodes0=conn_(iele0,surfs(cohel_surf(i,2),:)+1);
    nodes0=conn_(iele0,2:end);
    xyz0=mean(crds_(nodes0,2:end));
    sxyz0=crds_(snodes0,2:end);
    %
    iele1=cohel_surf(i,3);
    snodes1=conn_(iele1,surfs(cohel_surf(i,4),:)+1);
    nodes1=conn_(iele1,2:end);
    xyz1=mean(crds_(nodes1,2:end));
    sxyz1=crds_(snodes1,2:end);
    %
    % Loop though the nodes to find the matching coordinates
    snodes1_=zeros(1,size(snodes1,2));
    for inod0=1:size(snodes0,2)
        for inod1=1:size(snodes1,2)
            if norm(sxyz1(inod1,:)-sxyz0(inod0,:))<1e-19
                snodes1_(inod0)=snodes1(inod1);
            end
        end
    end
    snodes1=snodes1_;


%     % check the ordering of the surface nodes
%     if 
%         % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%         % DOES IT DIFFER FOR ELEMENT TYPES - WATCH FOR THIS!
%         % 2D linear quadrilateral - cohesive element
%         if (eltyp==6) || (eltyp==8)
%             snodes1=[snodes1(2), snodes1(1)];
%         % 2D quadratic quadrilateral - cohesive element
%         elseif (eltyp==12) || (eltyp==16)
%             snodes1=[snodes1(2), snodes1(1), snodes1(3)];
%         % 3D 8-node linear brick - cohesive element
%         elseif (eltyp==24)
% 
% 
%         end
%         %
%     end  
    %
    conn_coh=[conn_coh; ...
        i+tot_els,snodes0,snodes1];
    %
    % check the direction
    % normal direction
    ndir=(xyz1-xyz0)/norm(xyz1-xyz0);
    %
    % 2D linear quadrilateral - cohesive element
    if (eltyp==6) || (eltyp==8)
        %
        % coordinates of the unchanged nodes
        xyza=crds_(conn_coh(i,4),2:end);
        xyzb=crds_(conn_coh(i,5),2:end);
        %
        sdir=(xyzb-xyza)/norm(xyzb-xyza);
        %
        ndir_=[-sdir(2), sdir(1)];
        %
        sgn=dot(ndir,ndir_);
        %
        % Switch the order of connectivity
        if sgn<0
            %
            conn_coh(i,2:end)=[conn_coh(i,3), conn_coh(i,2), ...
                conn_coh(i,5), conn_coh(i,4)];
            %
        end
        %
    % 2D quadratic quadrilateral - cohesive element
    elseif ((eltyp==12) || (eltyp==16)) && numdim==2
        %
        % coordinates of the unchanged nodes
        xyza=crds_(conn_coh(i,5),2:end);
        xyzb=crds_(conn_coh(i,6),2:end);
        %
        sdir=(xyzb-xyza)/norm(xyzb-xyza);
        %
        ndir_=[-sdir(2), sdir(1)];
        %
        sgn=dot(ndir,ndir_);
        %
        % Switch the order of connectivity
        if sgn<0
            % 2 1 3 - 5 4 6
            conn_coh(i,2:end)=[conn_coh(i,3), conn_coh(i,2), conn_coh(i,4), ...
                conn_coh(i,6), conn_coh(i,5), conn_coh(i,7)];
            %
        end
        %
    % 3D linear triangular prism - cohesive element
    elseif eltyp==12 && numdim==3 % C3D4
         %
        % coordinates of the unchanged nodes
        xyza=crds_(conn_coh(i,5),2:end);
        xyzb=crds_(conn_coh(i,6),2:end);
        xyzc=crds_(conn_coh(i,7),2:end);
        %
        sdir=(xyzb-xyza)/norm(xyzb-xyza);
        tdir=(xyzc-xyza)/norm(xyzc-xyza);
        %
        ndir_=cross(sdir,tdir);
        %
        sgn=dot(ndir,ndir_);
        %
        % Switch the order of connectivity
        if sgn<0
            % 1 3 2 - 5 7 6
            conn_coh(i,2:end)=[conn_coh(i,2), conn_coh(i,4), conn_coh(i,3), ...
                conn_coh(i,5), conn_coh(i,7), conn_coh(i,6)];
            %
        end
        %
    % 3D quadratic triangular prism - cohesive element
    elseif eltyp==30 % C3D10
         %
        % coordinates of the unchanged nodes
        xyza=crds_(conn_coh(i,8),2:end);
        xyzb=crds_(conn_coh(i,9),2:end);
        xyzc=crds_(conn_coh(i,10),2:end);
        %
        sdir=(xyzb-xyza)/norm(xyzb-xyza);
        tdir=(xyzc-xyza)/norm(xyzc-xyza);
        %
        ndir_=cross(sdir,tdir);
        %
        sgn=dot(ndir,ndir_);
        %
        % Switch the order of connectivity
        if sgn<0
            % 1 3 2 6 5 4 - 7 9 8 12 11 10
            conn_coh(i,2:end)=[conn_coh(i,2), conn_coh(i,4), conn_coh(i,3), conn_coh(i,7), conn_coh(i,6), conn_coh(i,5),...
                conn_coh(i,8), conn_coh(i,10), conn_coh(i,9), conn_coh(i,13), conn_coh(i,12), conn_coh(i,11)];
            %
        end
        %
    % 3D linear brick - cohesive element
    elseif eltyp==24
         %
        % coordinates of the unchanged nodes
        xyza=crds_(conn_coh(i,6),2:end);
        xyzb=crds_(conn_coh(i,7),2:end);
        xyzc=crds_(conn_coh(i,9),2:end);
        %
        sdir=(xyzb-xyza)/norm(xyzb-xyza);
        tdir=(xyzc-xyza)/norm(xyzc-xyza);
        %
        ndir_=cross(sdir,tdir);
        %
        sgn=dot(ndir,ndir_);
        %
        % Switch the order of connectivity
        if sgn<0
            % 1 4 3 2 - 5 8 7 6
            % Note that the 1st column is the element number so,
            % "1" is added to each column number!
            conn_coh(i,2:end)=[conn_coh(i,2), conn_coh(i,5), conn_coh(i,4), conn_coh(i,3), ...
                conn_coh(i,6), conn_coh(i,9), conn_coh(i,8), conn_coh(i,7)];
            %
        end
        %
    % 3D quadratic - cohesive element
    elseif eltyp==60
         %
        % coordinates of the unchanged nodes
        xyza=crds_(conn_coh(i,6),2:end);
        xyzb=crds_(conn_coh(i,7),2:end);
        xyzc=crds_(conn_coh(i,9),2:end);
        %
        sdir=(xyzb-xyza)/norm(xyzb-xyza);
        tdir=(xyzc-xyza)/norm(xyzc-xyza);
        %
        ndir_=cross(sdir,tdir);
        %
        sgn=dot(ndir,ndir_);
        %
        % Switch the order of connectivity
        if sgn<0
            % 1 2 3 4 / 5 6 7 8   -   9 10 11 12 / 13 14 15 16
            % 1 4 3 2 / 8 7 6 5 - 9 12 11 10 / 16 15 14 13
            % Note that the 1st column is the element number so,
            % "1" is added to each column number!
            conn_coh(i,2:end)=[conn_coh(i,2), conn_coh(i,5), conn_coh(i,4), conn_coh(i,3), ...
                conn_coh(i,9), conn_coh(i,8), conn_coh(i,7), conn_coh(i,6), ...
                conn_coh(i,10), conn_coh(i,13), conn_coh(i,12), conn_coh(i,11), ...
                conn_coh(i,17), conn_coh(i,16), conn_coh(i,15), conn_coh(i,14)];
            %
        end
        %
    end
    %
end
%
%
%
%% Generate ABAQUS INP file
%
% NOTE: MUST BE EXTENDED FOR 3D CASE
% Find the boundary nodes
Xmin=min(crds_(:,2));
X0 = find(crds_(:,2)==Xmin);
% The indices of coordinates has changed
X0 = crds_(X0,1);
%
Ymin=min(crds_(:,3));
Y0 = find(crds_(:,3)==Ymin);
% The indices of coordinates has changed
Y0 = crds_(Y0,1);
%
Xmax=max(crds_(:,2));
X1 = find(crds_(:,2)==Xmax);
% The indices of coordinates has changed
X1 = crds_(X1,1);
%
% Incase of 3D fix also z-direction
if numdim==3
    Zmin=min(crds_(:,4));
    Z0 = find(crds_(:,4)==Zmin);
    % The indices of coordinates has changed
    Z0 = crds_(Z0,1);
end
%
%
% Write INP file
% Write the overall element and node sets to input file
% open inp file and write keywords 
inpFile = fopen([newfilename '.inp'],'wt');
fprintf(inpFile,'** Generated by Abaqus and modified by: CreateAbqINPFile.m\n');
fprintf(inpFile,'**PARTS\n**\n');
fprintf(inpFile,'*Part, name=Part-1\n');
%
% write nodes
fprintf(inpFile,'*NODE\n');
if numdim==2
    fprintf(inpFile,'%d,\t%e,\t%e\n',crds_');
elseif numdim==3
    fprintf(inpFile,'%d,\t%e,\t%e,\t%e\n',crds_');
end
%
% UEL definitions
fprintf(inpFile,['*User element, nodes=', num2str(nnps*2), ...
    ', type=U1, properties=9, coordinates=', num2str(numdim), ...
    ', variables=', num2str(numpt), '\n']);
if numdim==2
    fprintf(inpFile,' 1, 2\n');
elseif numdim==3
    fprintf(inpFile,' 1, 2, 3\n');
end
fprintf(inpFile,'**\n');
%
fprintf(inpFile,'*Element, type=U1, elset=UEL\n');
%
% CORR!!!
str=[];
for i=1:1:2*nnps
    str = [str, '%8d,'];
end
%
fprintf(inpFile,[str, '%8d\n'],conn_coh');
fprintf(inpFile,'*Uel Property, elset=UEL\n');
fprintf(inpFile, ...
    '\t%e, \t%e, \t%e, \t%e, \t%e, \t%e, \t%e, \t%e, \t%e\n', PROPS);
fprintf(inpFile,'*Elset, elset=UEL\n');
fprintf(inpFile,'%8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d\n',...
    (tot_els+1:1:tot_els+cohel)');
fprintf(inpFile,'\n**\n');
%
%
% write elements
fprintf(inpFile,['*Element, type=',elchar,'\n']);
str=[];
for i=1:1:nnpel
    str = [str, '%8d,'];
end
str = [str, '%8d\n'];
fprintf(inpFile,str,conn_');
%
%
% Write the elements sets for each grain to the input file
% create element sets containing grains
for ii = 1:numel(unique(grains))
    %
    fprintf(inpFile,'\n*Elset, elset=GRAIN-%d\n',grain_order(ii));
    fprintf(inpFile,'%d, %d, %d, %d, %d, %d, %d, %d, %d\n',conn(grains==grain_order(ii))');
    numels=0;
    %
    for tt=1:length(conn(grains==grain_order(ii)))
        %
        numels=numels+1;
    end
   numels_total(grain_order(ii))=numels;
end
%
% Write element set for each phase to input file
uniPhases = unique(phases);
for ii = 1:numel(unique(phases))
    fprintf(inpFile,'\n*Elset, elset=Phase-%d\n',ii);
    fprintf(inpFile,'%d, %d, %d, %d, %d, %d, %d, %d, %d\n',conn(phases==uniPhases(ii))');
end
%
% % Calculate grain spherical equivalent diameter
% % calulate diamater in microns
% % additionally, the dimaters for each ground are written to a separate
% % text file to be used to developed a grain size histogram
% diameterID=fopen('diameter.txt','w');
% for ii=1:numel(unique(grains))
%     
%     diameter(grain_order(ii))=((((6.0/pi)*(numels_total(grain_order(ii))))^(1/3)));
%     fprintf(diameterID, '%d\n', diameter(grain_order(ii)));
% end
% fclose(diameterID);
%
% write sections to each grain
for ii=1:length(grain_order)
    %
    fprintf(inpFile,'\n**Section: Section_Grain-%d\n*Solid Section, elset=GRAIN-%d, material=MATERIAL-GRAIN%d\n,\n',grain_order(ii),grain_order(ii),grain_order(ii));
end
% Continue writing the input file with assembly information
% write a closing keyword
fprintf(inpFile,'*End Part');
%
% writing assembly
fprintf(inpFile,'\n**\n**ASSEMBLY\n**');
fprintf(inpFile,'\n*Assembly, name=Assembly\n**');
fprintf(inpFile,'\n*Instance, name=Part-1-1, part=Part-1\n');
%
%
fprintf(inpFile,'\n*End Instance\n**');
%
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% NEEDS TO BE MODIFIED FOR THE 3D CASE
% Add the boundary conditions
fprintf(inpFile,'**\n');
fprintf(inpFile,'*Nset, nset=X0, internal, instance=Part-1-1\n');
fprintf(inpFile,'%d, %d, %d, %d, %d, %d, %d, %d, %d\n',X0');
%
fprintf(inpFile,'\n**\n');
fprintf(inpFile,'*Nset, nset=Y0, internal, instance=Part-1-1\n');
fprintf(inpFile,'%d, %d, %d, %d, %d, %d, %d, %d, %d\n',Y0');
%
if numdim==3
    fprintf(inpFile,'\n**\n');
    fprintf(inpFile,'*Nset, nset=Z0, internal, instance=Part-1-1\n');
    fprintf(inpFile,'%d, %d, %d, %d, %d, %d, %d, %d, %d\n',Z0');    
end
%
fprintf(inpFile,'\n**\n');
fprintf(inpFile,'*Nset, nset=X1, internal, instance=Part-1-1\n');
fprintf(inpFile,'%d, %d, %d, %d, %d, %d, %d, %d, %d\n',X1');
%
%
% Closing the assembly component of the input file
%
fprintf(inpFile,'\n*End Assembly');
%
%
fprintf(inpFile, '\n**MATERIALS\n**');
%
% import material parameters to be used in the development of materials
% for each grain.
xlRange='A1:A6';
[A16]=readmatrix('PROPS.xlsx','Sheet','Material_parameters','DataRange',xlRange);
%
%
% Finalising the input file
%
% Flag for reading the inputs from the file or material library
% "0": material library in usermaterial.f will be used
% "1": use the material parameters in excel file
%
% Are the material properties given in the excel file?
% Read from excel file (if read_all_props==true)
if A16(6)==0
    %   
    A = strings(6,1);
    %   
    % Flag for reading the inputs from the file or material library
    % "0": material library in usermaterial.f will be used
    % "1": use the material parameters in excel file
    A(6)=0;
    %
    % Number variables in DEPVAR
    noPROPS = 6;
    %
    % Do not define anything further
    %
else
    %
    A = strings(175,1);
    %
    %   
    % Flag for reading the inputs from the file or material library
    % "0": material library in usermaterial.f will be used
    % "1": use the material parameters in excel file
    A(6)=1;   
    %
    %
    % Number variables in DEPVAR
    % Has a fixed size - including additional space for extra variables
    noPROPS = 175;
    %
    %
end
%
%
%
for ii=1:length(grain_order)
    %
    fprintf(inpFile, '\n*Material, name=MATERIAL-GRAIN%d',grain_order(ii));
    fprintf(inpFile, ['\n*Depvar\n', num2str(noDepvar), ',']);
    fprintf(inpFile, ['\n*User Material, constants=',num2str(noPROPS),'\n']);
    %
    % Euler angles
    A(1:3) = [euler_angle1(ii), euler_angle2(ii), euler_angle3(ii)];
    % Grain - ID
    A(4) = grain_order(ii);
    %
    % Phase - ID
    % IF DEFINED BY THE USER (>0)
    if matID>0
        %
        A(5) = material_order(ii);
        %
    % Use Dream3D output 
    else
        %
        A(5) = phase_order(ii);
        %
    end
    %
    %
    % Read the properties from the PROPS if desired
    if A16(6) ==1
        %
        % Loop through all different phases
        for iph = 1:numel(unique(phases))
            %
            % Column character (read next column for each phase)
            letter = char(iph+ 64);
            %
            xlRange = [letter,'1:', letter, '175']; % A1-A175
            [B]=readmatrix('PROPS.xlsx','Sheet','Material_parameters','DataRange',xlRange);
            %
            A(7:175) = B(7:175);
            %
        end
        %
    end
    %
    % Printing this information to file
    fprintf(inpFile, '%s, %s, %s, %s, %s, %s, %s, %s\n',A);
    %
end
%
fprintf(inpFile,'\n**');
fprintf(inpFile, '\n**\n** STEP: Loading\n**\n*Step, name=Loading, nlgeom=YES, inc=10000\n*Static\n0.01, 10., 1e-05, 10.');
%
% BOUNDARY CONDITIONS
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% NEEDS TO BE MODIFIED FOR 3D CASE
fprintf(inpFile,'\n** BOUNDARY CONDITIONS');
fprintf(inpFile,'\n**');
fprintf(inpFile,'\n** Name: BC-1 Type: Symmetry/Antisymmetry/Encastre');
fprintf(inpFile,'\n*Boundary');
fprintf(inpFile,'\n X0, XSYMM');
%
fprintf(inpFile,'\n**');
fprintf(inpFile,'\n** Name: BC-2 Type: Symmetry/Antisymmetry/Encastre');
fprintf(inpFile,'\n*Boundary');
fprintf(inpFile,'\n Y0, YSYMM');
%
fprintf(inpFile,'\n**');
fprintf(inpFile,'\n** Name: BC-3 Type: Displacement/Rotation');
fprintf(inpFile,'\n*Boundary');
fprintf(inpFile,'\n X1, 1, 1, 0.1');
%
% for 3D case
if numdim==3
    fprintf(inpFile,'\n**');
    fprintf(inpFile,'\n** Name: BC-4 Type: Symmetry/Antisymmetry/Encastre');
    fprintf(inpFile,'\n*Boundary');
    fprintf(inpFile,'\n Z0, ZSYMM');
end
%
fprintf(inpFile,'\n**');
%
fprintf(inpFile, '\n**\n** CONTROLS\n**');
fprintf(inpFile, '\n*Controls, reset\n**');
%
fprintf(inpFile, '\n*Controls, parameters=field, field=displacement\n**');
fprintf(inpFile, '\n , 1., , , , , ,  \n**');
%
fprintf(inpFile, '\n*Controls, parameters=field, field=hydrostatic fluid pressure\n**');
fprintf(inpFile, '\n , 1., , , , , ,  \n**');
%
fprintf(inpFile, '\n*Controls, parameters=field, field=rotation\n**');
fprintf(inpFile, '\n , 1., , , , , ,  \n**');
%
fprintf(inpFile, '\n*Controls, parameters=field, field=electrical potential\n**');
fprintf(inpFile, '\n , 1., , , , , ,  \n**');
%
fprintf(inpFile, '\n*Controls, parameters=time incrementation\n**');
fprintf(inpFile, '\n , , , , , , , 10, , , \n**');
fprintf(inpFile,'\n**');
%
fprintf(inpFile, '\n**\n** OUTPUT REQUESTS\n**');
fprintf(inpFile, '\n*Restart, write, frequency=0\n**');
fprintf(inpFile, '\n** FIELD OUTPUT: F-Output-1\n**\n*Output, field, variable=PRESELECT\n**');
if noDepvar>0
    fprintf(inpFile, '\n** FIELD OUTPUT: F-Output-2\n**\n*Element Output, directions=YES\nSDV,\n**');
end
fprintf(inpFile, '\n** HISTORY OUTPUT: H-Output-1\n**\n*Output, history, variable=PRESELECT\n**');
fprintf(inpFile, '\n*End Step');
%
% close the file
fclose(inpFile);
%
%
toc
%
return
%