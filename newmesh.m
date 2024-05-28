clc;
clear all;

%Inputs for fan mesh
Lf = 5; % length of fan base
Nlf = 30; % division of fan base
Df = 5; %Depth of fan mesh
Ndf = 25 ; % division of fan dept(left an right)

%All inputs
B = 2.5;
D = 8;
x = (D*(Lf-B))/(Lf*Df)+1;
L = x*Lf;
Nd = 30; %division for symmetric boundary in main mesh
Nb = 30; %division for fan mesh right of footing or rectangular division

%Function call

    %generateMesh(B,L,D,Nd,Nb);
    generateFanMesh(Lf,Df,Nlf,Ndf,B,Nb,D,L,Nd)

%function to genrate fan mesh
function generateFanMesh(Lf,Df,Nlf,Ndf,B,Nb,D,L,Nd)
    
    y_coor_of_symm_side = linspace(0,-Df,Ndf+1); 
    y_coor_of_footing_edge_side = zeros(size(y_coor_of_symm_side)); % footing edge y-coor to connect to the symm side to draw radial lines
    x_coor_of_symm_side = zeros(size(y_coor_of_symm_side)); %x-xoor of symmetric side
    x_coor_of_footing_edge_side = B*ones(size(y_coor_of_footing_edge_side)); %x-coor of footing edge point
    
    plot([x_coor_of_footing_edge_side;x_coor_of_symm_side],[y_coor_of_footing_edge_side;y_coor_of_symm_side],'k-');
    hold on;
    
    %for drawing lines below footing
    x_coor_of_fan_base = linspace(0,Lf,Nlf+1);
    y_coor_of_fan_base = -Df*ones(size(x_coor_of_fan_base));
    x_coor_of_footing_edge_for_below = B*ones(size(x_coor_of_fan_base));
    y_coor_of_footing_edge_for_below = zeros(size(x_coor_of_footing_edge_for_below));
    plot([x_coor_of_footing_edge_for_below;x_coor_of_fan_base],[y_coor_of_footing_edge_for_below;y_coor_of_fan_base],'k-');
    
    %for drawing inclined lines on right side
    y_coor_of_rightside_fan = linspace(0,-(Df+D),Ndf+1);
    x_coor_of_rightside_fan = L*ones(size(y_coor_of_rightside_fan));
    plot([x_coor_of_footing_edge_side;x_coor_of_rightside_fan],[y_coor_of_footing_edge_side;y_coor_of_rightside_fan],'k-');
    
    %for drawaing rectangular lines
    x_coor_of_footing_base = linspace(0,B,Nb+1);
    x_coor_of_footing_base([1,Nb+1])=[];
    y_coor_of_footing_base = zeros(size(x_coor_of_footing_base));
    
    eqn_of_left_incl_line = polyfit([B,0],[0,-Df],1);
    eqn_of_right_incl_line = polyfit([B,Lf],[0,-Df],1);
    m_left=eqn_of_left_incl_line(1); c_left=eqn_of_left_incl_line(2);
    m_right = eqn_of_right_incl_line(1); c_right = eqn_of_right_incl_line(2);
    
    y_coor_of_left_incl_line = m_left*x_coor_of_footing_base + c_left;
    
    x_coor_of_right_incl_line = (y_coor_of_left_incl_line - c_right)/m_right;
    y_coor_of_footing_right = zeros(size(x_coor_of_right_incl_line));
    
    plot([x_coor_of_footing_base;x_coor_of_footing_base;x_coor_of_right_incl_line;x_coor_of_right_incl_line],[y_coor_of_footing_base;y_coor_of_left_incl_line;y_coor_of_left_incl_line;y_coor_of_footing_right],'k-');
    plot([0;0;Lf;Lf],[0;-Df;-Df;0],'k-');
    grid on;
    
    GetYCoordinateLeftRight(Nb,B,Ndf,Df,x_coor_of_footing_base,x_coor_of_footing_edge_side,y_coor_of_footing_edge_side,y_coor_of_symm_side,x_coor_of_symm_side,y_coor_of_left_incl_line);
    GetXCoordinateMiddle(Nb,Df,Lf, Nlf, y_coor_of_left_incl_line, x_coor_of_right_incl_line,x_coor_of_footing_edge_for_below,x_coor_of_fan_base,y_coor_of_footing_edge_for_below,y_coor_of_fan_base);
    GenerateMainMesh(B,Nlf,Ndf,Lf,Df,D,L,x_coor_of_fan_base,y_coor_of_fan_base,Nd);
    GetYCoordinate_For_Right_mesh(Nb,Nd,B,L,D,Df,Ndf,x_coor_of_right_incl_line,x_coor_of_footing_edge_side,y_coor_of_footing_edge_side,x_coor_of_rightside_fan,y_coor_of_rightside_fan,y_coor_of_left_incl_line);
end

%%function to generate Main mesh
%% 
function GenerateMainMesh(B,Nlf,Ndf,Lf,Df,D,L,x_coor_of_fan_base,y_coor_of_fan_base,Nd)
    x_coor_of_boundary_base = linspace(0,L,Nlf+1);
    y_coor_of_boundary_base = -(Df+D)*ones(size(x_coor_of_boundary_base));
    plot([x_coor_of_fan_base;x_coor_of_boundary_base],[y_coor_of_fan_base;y_coor_of_boundary_base],'k-');
    
    %To draw the rectangular lines in main mesh
    y_coor_of_symm_side = linspace(-Df,-(Df+D),Nd+1);
    y_coor_of_symm_side(1)=[];
    x_coor_of_symm_side = zeros(size(y_coor_of_symm_side));
    y_coor_right_of_footing = zeros(size(y_coor_of_symm_side));

    eqn_of_right_incl_line=polyfit([B,L],[0,-(Df+D)],1);
    m=eqn_of_right_incl_line(1); c=eqn_of_right_incl_line(2);
    x_coor_of_right_symm_line = (y_coor_of_symm_side-c)/m;
    plot([x_coor_of_symm_side;x_coor_of_right_symm_line;x_coor_of_right_symm_line],[y_coor_of_symm_side;y_coor_of_symm_side;y_coor_right_of_footing],'k-');

     GetXCoordinateForMainMeshBelow(Lf,Nlf, x_coor_of_fan_base, x_coor_of_boundary_base, y_coor_of_symm_side, D,Df, x_coor_of_symm_side, Nd, x_coor_of_right_symm_line)
    %disp(x_int_belowbase);
end
%% 
function  GetXCoordinateForMainMeshBelow(Lf,Nlf, x_coor_of_fan_base, x_coor_of_boundary_base, y_coor_of_symm_side, D,Df, x_coor_of_symm_side, Nd, x_coor_of_right_symm_line)
    y_coor_of_fan_base = -Df*ones(size(x_coor_of_fan_base));
    y_coor_of_boundary_base = -(Df+D)*ones(size(y_coor_of_fan_base));
    for i=1:Nlf
        eqn = polyfit([x_coor_of_fan_base(i),x_coor_of_boundary_base(i)],[y_coor_of_fan_base(i),y_coor_of_boundary_base(i)],1);
        m=eqn(1);
        c=eqn(2);
        x_coor_int(i,:)=(y_coor_of_symm_side-c)/m;
    end
    x_coor_int(1,:)=[];
    x_coor_int = x_coor_int';
    x_coor_of_fan_base([1,Nlf+1])=[];
    x_coor_of_symm_side=[0,x_coor_of_symm_side];
    x_coor_of_right_symm_line=[Lf,x_coor_of_right_symm_line];
    x_coor_int = [x_coor_of_fan_base;x_coor_int];
    x_coor_int = [x_coor_of_symm_side',x_coor_int,x_coor_of_right_symm_line'];
    
    y_coor_of_symm_side=[-Df,y_coor_of_symm_side];
    
    y_coor_int = repmat(y_coor_of_symm_side',1,Nlf+1);
    
    
    for i=1:Nd
        for j=1:Nlf
            Divide_Quadri(x_coor_int(i,j),x_coor_int(i+1,j),x_coor_int(i+1,j+1),x_coor_int(i,j+1),y_coor_int(i,j),y_coor_int(i+1,j),y_coor_int(i+1,j+1),y_coor_int(i,j+1));
        end
    end

end

%% function to get the Y_coord of the right side mesh
function GetYCoordinate_For_Right_mesh(Nb,Nd,B,L,D,Df,Ndf,x_coor_of_right_incl_line,x_coor_of_footing_edge_side,y_coor_of_footing_edge_side,x_coor_of_rightside_fan,y_coor_of_rightside_fan,y_coor_of_left_incl_line)
    eqn_of_right_inclined_line = polyfit([B,L],[0,-(D+Df)],1);
    y_coor_of_left_incl_line = flip(y_coor_of_left_incl_line);
    y_sym_main_mesh = linspace(-Df,-(D+Df),Nd+1);
    y_coor_right_incl_line = [y_coor_of_left_incl_line,y_sym_main_mesh];
    m_inclined = eqn_of_right_inclined_line(1); c_inclined=eqn_of_right_inclined_line(2);
    x_coor_of_right_incl_line = (y_coor_right_incl_line - c_inclined)/m_inclined;
    

    for i=1:Ndf
        eqn=polyfit([x_coor_of_footing_edge_side(i),x_coor_of_rightside_fan(i)],[y_coor_of_footing_edge_side(i),y_coor_of_rightside_fan(i)],1);
        m=eqn(1);
        c=eqn(2);
        y_coor_int(i,:) = m*x_coor_of_right_incl_line + c;
    end
    y_coor_int = [y_coor_int;y_coor_right_incl_line];
    x_coor_int = repmat(x_coor_of_right_incl_line,Ndf+1,1);
   
    for i=1:Ndf
        for j=1:Nb+Nd-1
            Divide_Quadri(x_coor_int(i,j),x_coor_int(i+1,j),x_coor_int(i+1,j+1),x_coor_int(i,j+1),y_coor_int(i,j),y_coor_int(i+1,j),y_coor_int(i+1,j+1),y_coor_int(i,j+1));
        end
    end
end

 

%% 

%function to get Y-coordinate of the interior of the fan mesh on right side.

function GetYCoordinateLeftRight(Nb,B,Ndf,Df,x_coor_of_footing_base,x_coor_of_footing_edge_side,y_coor_of_footing_edge_side,y_coor_of_symm_side,x_coor_of_symm_side,y_coor_of_left_incl_line)
    x_coor_of_footing_base=[x_coor_of_footing_base,B];
    x_coor_of_interior = repmat(x_coor_of_footing_base,Ndf+1,1);
    x_coor_of_interior = [x_coor_of_symm_side',x_coor_of_interior];
    

    
    x_coor_of_footing_edge_side([1,Ndf+1])=[];
    y_coor_of_footing_edge_side([1,Ndf+1])=[];
    y_coor_of_symm_side([1,Ndf+1])=[];
    x_coor_of_symm_side([1,Ndf+1])=[];
    %disp(y_coor_of_symm_side);
    
    
    for i=1:Ndf-1
        eqn = polyfit([x_coor_of_footing_edge_side(i),x_coor_of_symm_side(i)],[y_coor_of_footing_edge_side(i),y_coor_of_symm_side(i)],1);
        m= eqn(1); c=eqn(2);
        y_coor_interior(i,:)=m*x_coor_of_footing_base + c;
        %disp(y_coor_interior);
    end
    
    y_coor_size = size(y_coor_interior);
    y_coor_of_footing_base = zeros(1,y_coor_size(2));
    
    y_coor_of_left_incl_line = [y_coor_of_left_incl_line,0];
    y_coor_interior=[y_coor_of_footing_base;y_coor_interior;y_coor_of_left_incl_line];
    y_coor_int_size = size(y_coor_interior);
    y_coor_symm_side = linspace(0,-Df,Ndf+1);
    
    y_coor_interior=[y_coor_symm_side',y_coor_interior];
    y_coor_interior(:,y_coor_size(2)+1)=[];
    
    x_coor_of_interior(:,y_coor_size(2)+1)=[];
    
    %generate quads mesh
    for i=1:Ndf
        for j=1:Nb-1
            Divide_Quadri(x_coor_of_interior(i,j),x_coor_of_interior(i+1,j),x_coor_of_interior(i+1,j+1),x_coor_of_interior(i,j+1),y_coor_interior(i,j),y_coor_interior(i+1,j),y_coor_interior(i+1,j+1),y_coor_interior(i,j+1));
        end
    end
    
end

%function to get X-coordinate of the interior of the fan mesh on middle side.

function GetXCoordinateMiddle(Nb,Df,Lf, Nlf, y_coor_of_left_incl_line, x_coor_of_right_incl_line,  x_coor_of_footing_edge_for_below, x_coor_of_fan_base, y_coor_of_footing_edge_for_below,y_coor_of_fan_base)
    y_coor_of_left_incl_line = [-Df,y_coor_of_left_incl_line];
    x_coor_of_right_incl_line = [Lf,x_coor_of_right_incl_line];
    %disp(y_coor_of_left_incl_line);
    for i=1:Nlf
        eqn = polyfit([x_coor_of_footing_edge_for_below(i),x_coor_of_fan_base(i)],[y_coor_of_footing_edge_for_below(i),y_coor_of_fan_base(i)],1);
        m=eqn(1);
        c=eqn(2);
        x_interior(i,:)=(y_coor_of_left_incl_line - c)/m;
    end
    x_interior=[x_interior;x_coor_of_right_incl_line];
   
    y_interior = repmat(y_coor_of_left_incl_line,Nlf+1,1);
    
    for i=1:Nlf
        for j=1:Nb-1
            Divide_Quadri(x_interior(i,j),x_interior(i+1,j),x_interior(i+1,j+1),x_interior(i,j+1),y_interior(i,j),y_interior(i+1,j),y_interior(i+1,j+1),y_interior(i,j+1));
        end

    end
    



end

function generateMesh (B,L,D,Nd,Nb)
    %Horizontal boundary
    x_Footing_base = linspace(0, B, Nd+1); %Dividing base of footing in Nd segments
    x_H_boundarybase = linspace(0, L, Nd+1); % Dividing base of horizontal boundary in Nd segments
   
    y_H_boundarybase = -D*ones(size(x_H_boundarybase));
    y_footingbase = zeros(size(x_Footing_base));
    plot([x_Footing_base;x_H_boundarybase], [y_footingbase; y_H_boundarybase], 'k-'); %Plots vertical lines
    hold on;
    
    %Vertical boundary
    V_boundarybase = linspace(0, -D, Nd+1); %This is the y-coor of the right side vertical boundary
    x_V_footingbase = B*ones(size(V_boundarybase)); %x-coor array of footing edge
    x_V_boundarybase = L*ones(size(V_boundarybase)); %x-coor array of right vertical boundary edge
    y_V_footingbase = zeros(size(V_boundarybase)); %y-coor array of footing edge.
    
    plot([x_V_footingbase;x_V_boundarybase], [y_V_footingbase;V_boundarybase], 'k-'); %Plots inclined lines
    
    %Vertical Lines
    start = log10(B);
    stop = log10(L);
    x_top_boundary = logspace(start, stop, Nb+1);
    m = D/(B-L);
    c = (-D*B)/(B-L);
    
    y_diagonal = x_top_boundary*m + c; %This is same as the symmetric side y-coordinate
    y_top_boundary = zeros(size(x_top_boundary));
    x_symmetric_boundary = zeros(size(x_top_boundary)); %For below footing space
    
    plot([x_top_boundary;x_top_boundary],[y_top_boundary;y_diagonal],'k-'); %Plots vertical lines
    plot([x_top_boundary;x_symmetric_boundary],[y_diagonal;y_diagonal],'k-'); %Plots horizontal lines below footing space
    
    grid on;
    
    
   % x-coord of intersection points below footing base
   [x_int_belowbase]= GetXCoordinatesBelowFooting(Nd, x_Footing_base, x_H_boundarybase, y_diagonal, D, x_symmetric_boundary, Nb, x_top_boundary);
   
   % y-coord of intersection points right of footing base
   [y_int_right_base] = GetYCoordinatesRightOfFooting(Nd,B,L, V_boundarybase, x_top_boundary, Nb);
   

end


function [y_int_right_base] = GetYCoordinatesRightOfFooting(Nd,B,L,V_boundarybase,x_top_boundary,Nb)
    
    %Trimming first and last values of arrays as we already have those
    %coordinates

    V_boundarybase([1,Nd+1])=[]; 
    x_top_boundary([1,Nb+1]) = [];

    for i = 1:Nd-1
        
        X = [B, L];
        Y = [0, V_boundarybase(i)];

        coefficient = polyfit(X, Y, 1);

        m = coefficient(1);
        c = coefficient(2);
        y_int_right_base(:, i) = m*x_top_boundary +c; % y-coord of intersection points right of footing base

        %Each column = one inclined line from top
        %Each row = Vertical line intersect from footing edge
        %Ex: 3,2 = 3rd vertical line intersecting with 2nd inclined line
        
    end
end

%%Function to divide the quads into triangles and store intersecting
%%point coordinates
function [x_coord_midpoint,y_coord_midpoint] = Divide_Quadri (x1,x2,x3,x4,y1,y2,y3,y4)
    plot([x1;x3], [y1;y3], 'k-');
    plot([x2;x4], [y2;y4], 'k-');
    
    %To find intersecting point coordinates
    X1=[x1,x3];
    X2 = [x2,x4];
    Y1=[y1,y3];
    Y2=[y2,y4];
    
    eqn1 = polyfit(X1,Y1,1);
    eqn2 = polyfit(X2,Y2,1);
    
    m1=eqn1(1);m2=eqn2(1);
    c1=eqn1(2);c2=eqn2(2);
    
    x_coord_midpoint = (c1-c2)/(m2-m1);
    y_coord_midpoint = m2*x_coord_midpoint + c2;
    %disp(x_coord_midpoint+" "+y_coord);
end









