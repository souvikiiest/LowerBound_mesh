clc;
clear all;
close all;

%% All inputs
B = 1;
Nlf = 5; % division of fan base
Df = B; %Depth of fan mesh
Ndf = Nlf ; % division of fan dept(left and right)
L = 8*B;

D = 5*B;
%x = (D*(Lf-B))/(Lf*Df)+1;
Lf=Df*(L-B)/(Df+D) + B; % length of fan base
Nb = Ndf; %division for fan mesh right of footing or rectangular division
Nd = 5*Nb; %division for symmetric boundary in main mesh
r=2; % ratio for G.P series

%% Function call

   total_node_table = generateFanMesh(Lf,Df,Nlf,Ndf,B,Nb,D,L,Nd,r);

%% function to genrate fan mesh
function [total_node_table] = generateFanMesh(Lf,Df,Nlf,Ndf,B,Nb,D,L,Nd,r)
    
    y_coor_of_symm_side = linspace(0, -Df, Ndf+1); 
    y_coor_of_footing_edge_side = zeros(size(y_coor_of_symm_side)); % footing edge y-coor to connect to the symm side to draw radial lines
    x_coor_of_symm_side = zeros(size(y_coor_of_symm_side)); %x-xoor of symmetric side
    x_coor_of_footing_edge_side = B*ones(size(y_coor_of_footing_edge_side)); %x-coor of footing edge point
    
    plot([x_coor_of_footing_edge_side;x_coor_of_symm_side],[y_coor_of_footing_edge_side;y_coor_of_symm_side],'k-');
    hold on;
    
    %for drawing lines below footing
    
    x_coor_of_fan_base = linspace(0, Lf, Nlf+1);
    y_coor_of_fan_base = -Df*ones(size(x_coor_of_fan_base));
    x_coor_of_footing_edge_for_below = B*ones(size(x_coor_of_fan_base));
    y_coor_of_footing_edge_for_below = zeros(size(x_coor_of_footing_edge_for_below));
    plot([x_coor_of_footing_edge_for_below;x_coor_of_fan_base],[y_coor_of_footing_edge_for_below;y_coor_of_fan_base],'k-');
    
    %for drawing inclined lines on right side
    y_coor_of_rightside_fan = linspace(0, -(Df+D), Ndf+1);
    x_coor_of_rightside_fan = L*ones(size(y_coor_of_rightside_fan));
    plot([x_coor_of_footing_edge_side;x_coor_of_rightside_fan],[y_coor_of_footing_edge_side;y_coor_of_rightside_fan],'k-');
    
    %for drawaing rectangular lines
    %[x_coords,~]= divide_line_gp(B, 0, 0, 0, Nb, r);
    x_coor_of_footing_base = linspace(0, B, Nb+1);
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
    
    [node_left_fan_total] = GetYCoordinateLeftRight(Nb,B,Ndf,Df,x_coor_of_footing_base,x_coor_of_footing_edge_side,y_coor_of_footing_edge_side,y_coor_of_symm_side,x_coor_of_symm_side,y_coor_of_left_incl_line,r);
    [node_middle_fan] = GetXCoordinateMiddle(B, Nb,Df,Lf, Nlf, y_coor_of_left_incl_line, x_coor_of_right_incl_line,x_coor_of_footing_edge_for_below,x_coor_of_fan_base,y_coor_of_footing_edge_for_below,y_coor_of_fan_base);
    [node_middle_main] = GenerateMainMesh(B,Nlf,Ndf,Lf,Df,D,L,x_coor_of_fan_base,y_coor_of_fan_base,Nd,r);
    [node_right_side] = GetYCoordinate_For_Right_mesh(Nb,Nd,B,L,D,Df,Ndf,x_coor_of_right_incl_line,x_coor_of_footing_edge_side,y_coor_of_footing_edge_side,x_coor_of_rightside_fan,y_coor_of_rightside_fan,y_coor_of_left_incl_line,r);

    total_node_table = [node_left_fan_total; node_middle_main;node_middle_fan; node_right_side];
end

%% function to generate Main mesh 
function [node_middle_main] = GenerateMainMesh(B,Nlf,Ndf,Lf,Df,D,L,x_coor_of_fan_base,y_coor_of_fan_base,Nd,r)
    %[x_coords,~]= divide_line_gp(0, -(Df+D), L, -(Df+D), Nlf, r);
    x_coor_of_boundary_base = linspace(0, L, Nlf+1);
    y_coor_of_boundary_base = -(Df+D)*ones(size(x_coor_of_boundary_base));
    plot([x_coor_of_fan_base;x_coor_of_boundary_base],[y_coor_of_fan_base;y_coor_of_boundary_base],'k-');
    
    %To draw the rectangular lines in main mesh
    %[~,y_coords]= divide_line_gp(0, -Df, 0, -(Df+D), Nd, r);
    y_coor_of_symm_side = linspace(-Df, -(Df+D), Nd+1);
    y_coor_of_symm_side(1)=[];
    x_coor_of_symm_side = zeros(size(y_coor_of_symm_side));
    y_coor_right_of_footing = zeros(size(y_coor_of_symm_side));

    eqn_of_right_incl_line=polyfit([B,L],[0,-(Df+D)],1);
    m=eqn_of_right_incl_line(1); c=eqn_of_right_incl_line(2);
    x_coor_of_right_symm_line = (y_coor_of_symm_side-c)/m;
    plot([x_coor_of_symm_side;x_coor_of_right_symm_line;x_coor_of_right_symm_line],[y_coor_of_symm_side;y_coor_of_symm_side;y_coor_right_of_footing],'k-');

    [node_middle_main] = GetXCoordinateForMainMeshBelow(B,Lf,Nlf, x_coor_of_fan_base, x_coor_of_boundary_base, y_coor_of_symm_side, D,Df, x_coor_of_symm_side, Nd, x_coor_of_right_symm_line);
    %disp(x_int_belowbase);
end

%% function to get the X_coord of the main mesh below
function [node_middle_main] = GetXCoordinateForMainMeshBelow(B,Lf,Nlf, x_coor_of_fan_base, x_coor_of_boundary_base, y_coor_of_symm_side, D,Df, x_coor_of_symm_side, Nd, x_coor_of_right_symm_line)
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
            [x_coord_midpoint,y_coord_midpoint] = Divide_Quadri(x_coor_int(i,j),x_coor_int(i+1,j),x_coor_int(i+1,j+1),x_coor_int(i,j+1),y_coor_int(i,j),y_coor_int(i+1,j),y_coor_int(i+1,j+1),y_coor_int(i,j+1));
            x_coord_midpoint_main_middle(i,j) = x_coord_midpoint;
            y_coord_midpoint_main_middle(i,j) = y_coord_midpoint;
        end
    end
    
    [node_middle_main] = getNodes(B,flipud(x_coor_int)', flipud(y_coor_int)', flipud(x_coord_midpoint_main_middle)', flipud(y_coord_midpoint_main_middle)', false);
 
end

%% Function to get the Y_coord of the right side main mesh
function [node_right_side] = GetYCoordinate_For_Right_mesh(Nb,Nd,B,L,D,Df,Ndf,x_coor_of_right_incl_line,x_coor_of_footing_edge_side,y_coor_of_footing_edge_side,x_coor_of_rightside_fan,y_coor_of_rightside_fan,y_coor_of_left_incl_line,r)
    eqn_of_right_inclined_line = polyfit([B,L],[0,-(D+Df)],1);
    y_coor_of_left_incl_line = flip(y_coor_of_left_incl_line);

    %[~,y_coords]= divide_line_gp(0, -Df, 0, -(Df+D), Nd, r);
    y_sym_main_mesh = linspace(-Df, -(Df+D), Nd+1);
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
            [x_coord_midpoint,y_coord_midpoint] = Divide_Quadri(x_coor_int(i,j),x_coor_int(i+1,j),x_coor_int(i+1,j+1),x_coor_int(i,j+1),y_coor_int(i,j),y_coor_int(i+1,j),y_coor_int(i+1,j+1),y_coor_int(i,j+1));
            x_coord_midpoint_fan_right(i,j) = x_coord_midpoint;
            y_coord_midpoint_fan_right(i,j) = y_coord_midpoint;
        end
    end
    
    % Store the node coordinates

    [node_right_side] = getNodes_right(B,fliplr(x_coor_int), fliplr(y_coor_int), fliplr(x_coord_midpoint_fan_right), fliplr(y_coord_midpoint_fan_right));
end
 

%% Function to get Y-coordinate of the interior of the fan mesh on left side

function [node_left_fan_total] = GetYCoordinateLeftRight(Nb,B,Ndf,Df,x_coor_of_footing_base,x_coor_of_footing_edge_side,y_coor_of_footing_edge_side,y_coor_of_symm_side,x_coor_of_symm_side,y_coor_of_left_incl_line,r)
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

    %[~,y_coords]= divide_line_gp(0, 0, 0, -Df, Ndf, r);
    y_coor_symm_side = linspace(0, -Df, Ndf+1);
    
    y_coor_interior=[y_coor_symm_side',y_coor_interior];
    y_coor_interior(:,y_coor_size(2)+1)=[];
    
    x_coor_of_interior(:,y_coor_size(2)+1)=[];
    
    %generate quads mesh
    for i=1:Ndf
        for j=1:Nb-1
           [x_coord_midpoint,y_coord_midpoint] = Divide_Quadri(x_coor_of_interior(i,j),x_coor_of_interior(i+1,j),x_coor_of_interior(i+1,j+1),x_coor_of_interior(i,j+1),y_coor_interior(i,j),y_coor_interior(i+1,j),y_coor_interior(i+1,j+1),y_coor_interior(i,j+1));
            x_coord_midpoint_fan_left(i,j) = x_coord_midpoint;
            y_coord_midpoint_fan_left(i,j) = y_coord_midpoint;
        end
    end
    
    % Store the node coordinates

    [node_left_fan_total] = getNodes(B,x_coor_of_interior, y_coor_interior, x_coord_midpoint_fan_left, y_coord_midpoint_fan_left);
   
end

%% Function to get X-coordinate of the interior of the fan mesh on middle side.

function [node_middle_fan] = GetXCoordinateMiddle(B, Nb,Df,Lf, Nlf, y_coor_of_left_incl_line, x_coor_of_right_incl_line,  x_coor_of_footing_edge_for_below, x_coor_of_fan_base, y_coor_of_footing_edge_for_below,y_coor_of_fan_base)
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
           [x_coord_midpoint,y_coord_midpoint] = Divide_Quadri(x_interior(i,j),x_interior(i+1,j),x_interior(i+1,j+1),x_interior(i,j+1),y_interior(i,j),y_interior(i+1,j),y_interior(i+1,j+1),y_interior(i,j+1));
           x_coord_midpoint_fan_base(i,j) = x_coord_midpoint;
           y_coord_midpoint_fan_base(i,j) = y_coord_midpoint;
        end

    end
    
    % Store the node coordinates

    [node_middle_fan] = getNodes(B,x_interior, y_interior, x_coord_midpoint_fan_base, y_coord_midpoint_fan_base);
   
end


%% Function to divide the quads into triangles and store intersecting
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


%% function to get Nodes in arranged way (for left and middle)

function [node_left_fan] = getNodes(B,x_coor_of_interior, y_coor_interior, x_coord_midpoint_fan_left, y_coord_midpoint_fan_left, include_last_triangle)
    mid_point_size = size(x_coord_midpoint_fan_left);

    if nargin < 6
        include_last_triangle = true; % Default behavior is to include the last extra triangle
    
    end

    node_left_fan=[];
    flag=1;
    
    for i=1:mid_point_size(1)
        for j=1:mid_point_size(2)
            %first element
            node_left_fan(flag,1)=x_coord_midpoint_fan_left(i,j); node_left_fan(flag,2)=y_coord_midpoint_fan_left(i,j);flag=flag+1;
            node_left_fan(flag,1) = x_coor_of_interior(i,j); node_left_fan(flag,2) = y_coor_interior(i,j); flag=flag+1;
            node_left_fan(flag,1) = x_coor_of_interior(i+1,j); node_left_fan(flag,2) = y_coor_interior(i+1,j); flag=flag+1;
            %Second element
            node_left_fan(flag,1)=x_coord_midpoint_fan_left(i,j); node_left_fan(flag,2)=y_coord_midpoint_fan_left(i,j);flag=flag+1;
            node_left_fan(flag,1) = x_coor_of_interior(i+1,j); node_left_fan(flag,2) = y_coor_interior(i+1,j); flag=flag+1;
            node_left_fan(flag,1) = x_coor_of_interior(i+1,j+1); node_left_fan(flag,2) = y_coor_interior(i+1,j+1); flag=flag+1;
            % Third element
            node_left_fan(flag,1)=x_coord_midpoint_fan_left(i,j); node_left_fan(flag,2)=y_coord_midpoint_fan_left(i,j);flag=flag+1;
            node_left_fan(flag,1) = x_coor_of_interior(i+1,j+1); node_left_fan(flag,2) = y_coor_interior(i+1,j+1); flag=flag+1;
            node_left_fan(flag,1) = x_coor_of_interior(i,j+1); node_left_fan(flag,2) = y_coor_interior(i,j+1); flag=flag+1;
            % Fourth element
            node_left_fan(flag,1)=x_coord_midpoint_fan_left(i,j); node_left_fan(flag,2)=y_coord_midpoint_fan_left(i,j);flag=flag+1;
            node_left_fan(flag,1) = x_coor_of_interior(i,j+1); node_left_fan(flag,2) = y_coor_interior(i,j+1); flag=flag+1;
            node_left_fan(flag,1) = x_coor_of_interior(i,j); node_left_fan(flag,2) = y_coor_interior(i,j); flag=flag+1;

            % For last extra triangle 
            if include_last_triangle && j==mid_point_size(2)
                node_left_fan(flag,1)=B; node_left_fan(flag,2)=0;flag=flag+1;
                node_left_fan(flag,1)=x_coor_of_interior(i,j+1); node_left_fan(flag,2)=y_coor_interior(i,j+1);flag=flag+1;
                node_left_fan(flag,1) = x_coor_of_interior(i+1,j+1); node_left_fan(flag,2) = y_coor_interior(i+1,j+1); flag=flag+1;
            end
        end
    end
end


%% function to get Nodes in arranged way (for right side)

function [node_right_fan] = getNodes_right(B,x_coor_of_interior, y_coor_interior, x_coord_midpoint_fan_left, y_coord_midpoint_fan_left, include_last_triangle)
    mid_point_size = size(x_coord_midpoint_fan_left);
    
    if nargin < 6
        include_last_triangle = true; % Default behavior is to include the last extra triangle
    
    end
    node_right_fan=[];
    flag=1;
    
    for i=1:mid_point_size(1)
        for j=1:mid_point_size(2)
            %first element
            node_right_fan(flag,1)=x_coord_midpoint_fan_left(i,j); node_right_fan(flag,2)=y_coord_midpoint_fan_left(i,j);flag=flag+1;
            node_right_fan(flag,1) = x_coor_of_interior(i+1,j); node_right_fan(flag,2) = y_coor_interior(i+1,j); flag=flag+1;
            node_right_fan(flag,1) = x_coor_of_interior(i,j); node_right_fan(flag,2) = y_coor_interior(i,j); flag=flag+1;
            %Second element
            node_right_fan(flag,1)=x_coord_midpoint_fan_left(i,j); node_right_fan(flag,2)=y_coord_midpoint_fan_left(i,j);flag=flag+1;
            node_right_fan(flag,1) = x_coor_of_interior(i,j); node_right_fan(flag,2) = y_coor_interior(i,j); flag=flag+1;
            node_right_fan(flag,1) = x_coor_of_interior(i,j+1); node_right_fan(flag,2) = y_coor_interior(i,j+1); flag=flag+1;
            % Third element
            node_right_fan(flag,1)=x_coord_midpoint_fan_left(i,j); node_right_fan(flag,2)=y_coord_midpoint_fan_left(i,j);flag=flag+1;
            node_right_fan(flag,1) = x_coor_of_interior(i,j+1); node_right_fan(flag,2) = y_coor_interior(i,j+1); flag=flag+1;
            node_right_fan(flag,1) = x_coor_of_interior(i+1,j+1); node_right_fan(flag,2) = y_coor_interior(i+1,j+1); flag=flag+1;
            % Fourth element
            node_right_fan(flag,1)=x_coord_midpoint_fan_left(i,j); node_right_fan(flag,2)=y_coord_midpoint_fan_left(i,j);flag=flag+1;
            node_right_fan(flag,1) = x_coor_of_interior(i+1,j+1); node_right_fan(flag,2) = y_coor_interior(i+1,j+1); flag=flag+1;
            node_right_fan(flag,1) = x_coor_of_interior(i+1,j); node_right_fan(flag,2) = y_coor_interior(i+1,j); flag=flag+1;

            % For last extra triangle 
            if include_last_triangle && j==mid_point_size(2)
                node_right_fan(flag,1)=B; node_right_fan(flag,2)=0;flag=flag+1;
                node_right_fan(flag,1)=x_coor_of_interior(i+1,j+1); node_right_fan(flag,2)=y_coor_interior(i+1,j+1);flag=flag+1;
                node_right_fan(flag,1) = x_coor_of_interior(i,j+1); node_right_fan(flag,2) = y_coor_interior(i,j+1); flag=flag+1;
            end
        end
    end
end

%% Function for G.P series
function [x_coords,y_coords]= divide_line_gp(x0, y0, x1, y1, n, r)
    % Length of the line segment
    L = sqrt((x1 - x0)^2 + (y1 - y0)^2);
    
    % Calculate the first term 'a' of the GP
    a = L * (r - 1) / (r^n - 1);
    
    % Initialize arrays to store the coordinates
    x_coords = zeros(1, n+1);
    y_coords = zeros(1, n+1);
    
    % Set the starting point
    x_coords(1) = x0;
    y_coords(1) = y0;
    
    % Calculate the cumulative distance along the line
    cumulative_distance = 0;
    for k = 1:n
        segment_length = a * r^(k-1);
        cumulative_distance = cumulative_distance + segment_length;
        
        % Calculate the coordinates of the k-th division point
        x_coords(k+1) = x0 + (cumulative_distance / L) * (x1 - x0);
        y_coords(k+1) = y0 + (cumulative_distance / L) * (y1 - y0);
    end
end


