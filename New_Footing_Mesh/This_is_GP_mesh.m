close all;
clear all;
warning('off','all');

%Inputs for fan mesh
Nlf = 7; % division of fan base
%Depth of fan mesh
Ndf = Nlf; % division of fan dept(left an right)

%All inputs
B = 1;
Df = 1.5*B;
D = 10*B;
L = 10*B; %Total extent of boundary
Lf = Df*(L-B)/(Df+D) + B; % length of fan base

Nd = 2*Nlf; %division for symmetric boundary in main mesh 2*Nlf %causing problem in mesh for 16,10
Nb = Nlf; %division for fan mesh right of footing or rectangular division Nlf

p=24; %for yield function
fi=25;
c=0.0001;
gamma=18;
tic;
%% Function call
[total_node_table,no_of_element]=generateFanMesh(Lf,Df,Nlf,Ndf,B,Nb,D,L,Nd);
total_node_table = round(total_node_table,3);

[result,new_node] = findMatchingEdges(total_node_table); %result has the discontinuity planes with angles
[left_bound_result,right_bound_result,middle_bound_result] = findConsecutiveNodes(total_node_table, B);
[A_element,B_element, A_Yield, B_yield,A_bound,B_Boundary,A_discontinuity,B_discontinuity,C_matrix] = ...
    EquilibriumParentFunction(total_node_table,no_of_element,left_bound_result,...
    right_bound_result,middle_bound_result,result,p,fi,c,gamma);

%% This are some checking function to verify nodes,boundary
 %labelTriangles(total_node_table);
 %plotNodesWithNumbersInside(total_node_table);
 %DrawDiscontinuity(total_node_table,result);
 %checkNodes(total_node_table);
 %checkElementEquilibirum(total_node_table);
 %DrawBoundary(left_bound_result,middle_bound_result,right_bound_result)



A2 = [A_element;A_discontinuity;A_bound];
B2 = [B_element;B_discontinuity;B_Boundary];
A1 = A_Yield;
B1 = B_yield;

f=full(C_matrix');
%% using gurobi
% A=[A1;A2];
% b=[B1;B2];
% model.A = sparse(A);
% model.obj = -f;
% model.rhs = full(b);
% model.sense = [repmat('<', size(A1,1),1); repmat('=', size(A2,1),1)];
% model.modelsense ='min';
% params.outputflag = 0;
% params.Threads = 4;
% params.Method = 1;
% result = gurobi(model,params);
% solution = result.x;
% fval = result.objval;

%% using linprog

 %options = optimset('MaxIter',10000000,'Display','final','TolFun',1e-3,'TolX',1e-3);
 [solution,fval,exitflag,output]=linprog(f,A1,B1,A2,B2); %,[],[],options

 %% Display results

N_c= -fval/(c*B);
N_g = -(fval)/(gamma*B*B);
N_q = -(fval);
fprintf("p = %d, c = %d, fi = %d \n",p,c,fi);
disp("Collapse Load: "+ -fval);
disp("Nc: "+N_c);
disp("Nq: "+N_q);
disp("Ng: "+N_g);
time = toc;
disp("Elapsed Time: "+time);

%% Equilibrium parent function
function [A_element,B_element, A_Yield, B_yield,A_bound,B_Boundary,A_discontinuity,...
    B_discontinuity,C_matrix] = EquilibriumParentFunction(total_node_table,...
    no_of_element,left_bound_result,right_bound_result,middle_bound_result,result,p,fi,c,gamma)

    [A_element,B_element]=ElementEquilibrium(total_node_table,no_of_element,gamma);
    [A_Yield, B_yield] = YieldEquilibrium(no_of_element,p,fi,c);
    [A_bound,B_Boundary] = BoundaryCondition(left_bound_result,middle_bound_result,right_bound_result,no_of_element);
    [A_discontinuity,B_discontinuity] = DiscontinuityEquilibrium(result,no_of_element);
     C_matrix = ObjectiveFunction(middle_bound_result,no_of_element);
end

%% function to genrate fan mesh
function [total_node_table,no_of_element] = generateFanMesh(Lf,Df,Nlf,Ndf,B,Nb,D,L,Nd)

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
plot([0;0;Lf;Lf],[0;-Df;-Df;0],'b-');
grid on;

[node_left_fan_total] = GetYCoordinateLeftRight(Nb,B,Ndf,Df,x_coor_of_footing_base,x_coor_of_footing_edge_side,y_coor_of_footing_edge_side,y_coor_of_symm_side,x_coor_of_symm_side,y_coor_of_left_incl_line);
[node_middle_fan] =GetXCoordinateMiddle(B,Nb,Df,Lf, Nlf, y_coor_of_left_incl_line, x_coor_of_right_incl_line,x_coor_of_footing_edge_for_below,x_coor_of_fan_base,y_coor_of_footing_edge_for_below,y_coor_of_fan_base);
[y_coor_of_symm_side,node_middle_main]= GenerateMainMesh(B,Nlf,Ndf,Lf,Df,D,L,x_coor_of_fan_base,y_coor_of_fan_base,Nd);
[node_right_side] =GetYCoordinate_For_Right_mesh(y_coor_of_symm_side,Nb,Nd,B,L,D,Df,Ndf,x_coor_of_footing_edge_side,y_coor_of_footing_edge_side,x_coor_of_rightside_fan,y_coor_of_rightside_fan,y_coor_of_left_incl_line);

total_node_table = [node_left_fan_total; node_middle_main;node_middle_fan; node_right_side];
sized=size(total_node_table);
disp("Total number of element: "+sized(1)/3);
no_of_element = sized(1)/3;
%disp(total_node_table);

end

%% function to generate Main mesh
function [y_coor_of_symm_side,node_middle_main]= GenerateMainMesh(B,Nlf,Ndf,Lf,Df,D,L,x_coor_of_fan_base,y_coor_of_fan_base,Nd)
x_coor_of_boundary_base = linspace(0,L,Nlf+1);
y_coor_of_boundary_base = -(Df+D)*ones(size(x_coor_of_boundary_base));
plot([x_coor_of_fan_base;x_coor_of_boundary_base],[y_coor_of_fan_base;y_coor_of_boundary_base],'k-');

%To draw the rectangular lines in main mesh
initial_length = Df/Ndf;
initial_increase_ratio = 1.1;
initial_length_of_main_mesh = initial_increase_ratio*initial_length; %this is 'a' for the gp series on the Depth part.

%gp-series
r = findRateOfIncrease((D-Df),initial_length_of_main_mesh,Nd);
% disp(r);
a=initial_length_of_main_mesh; sum=a;coord=zeros(1,Nd);coord(1)=a;
Length = D-Df;i=2;
while(sum<=Length)
    sum=sum+ (a*r^(i-1));
    coord(i) = sum;
    i=i+1;
end
coord = -coord - Df;
%coord(end) = -(D+Df);


%y_coor_of_symm_side = linspace(-Df,-(Df+D),Nd+1);
y_coor_of_symm_side=coord;
x_coor_of_symm_side = zeros(size(y_coor_of_symm_side));
y_coor_right_of_footing = zeros(size(y_coor_of_symm_side));

eqn_of_right_incl_line=polyfit([B,L],[0,-(Df+D)],1);
m=eqn_of_right_incl_line(1); c=eqn_of_right_incl_line(2);
x_coor_of_right_symm_line = (y_coor_of_symm_side-c)/m;
plot([x_coor_of_symm_side;x_coor_of_right_symm_line;x_coor_of_right_symm_line],[y_coor_of_symm_side;y_coor_of_symm_side;y_coor_right_of_footing],'k-');

[node_middle_main]=GetXCoordinateForMainMeshBelow(B,Lf,Nlf, x_coor_of_fan_base, x_coor_of_boundary_base, y_coor_of_symm_side, D,Df, x_coor_of_symm_side, Nd, x_coor_of_right_symm_line);
%disp(node_middle_main);
end

%%
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
%disp(x_coor_of_right_symm_line);
x_coor_of_right_symm_line=[Lf,x_coor_of_right_symm_line];
%disp(x_coor_of_right_symm_line);

x_coor_int = [x_coor_of_fan_base;x_coor_int];

x_coor_int = [x_coor_of_symm_side',x_coor_int,x_coor_of_right_symm_line'];
%disp(x_coor_int);
y_coor_of_symm_side=[-Df,y_coor_of_symm_side];

y_coor_int = repmat(y_coor_of_symm_side',1,Nlf+1);
%disp(y_coor_int);
count=0;
for i=1:Nd
    for j=1:Nlf
        [x_coord_midpoint,y_coord_midpoint] = Divide_Quadri(x_coor_int(i,j),x_coor_int(i+1,j),x_coor_int(i+1,j+1),x_coor_int(i,j+1),y_coor_int(i,j),y_coor_int(i+1,j),y_coor_int(i+1,j+1),y_coor_int(i,j+1));
        x_coord_midpoint_main_middle(i,j) = x_coord_midpoint;
        y_coord_midpoint_main_middle(i,j) = y_coord_midpoint;
        count=count+1;
    end
end
[node_middle_main] = getNodes(B,flipud(x_coor_int)', flipud(y_coor_int)', flipud(x_coord_midpoint_main_middle)', flipud(y_coord_midpoint_main_middle)', false);
%disp("the count of main mesh below:"+count);
end

%% function to get the Y_coord of the right side mesh
function [node_right_side] = GetYCoordinate_For_Right_mesh(y_coor_of_symm_side,Nb,Nd,B,L,D,Df,Ndf,x_coor_of_footing_edge_side,y_coor_of_footing_edge_side,x_coor_of_rightside_fan,y_coor_of_rightside_fan,y_coor_of_left_incl_line)


eqn_of_right_inclined_line = polyfit([B,L],[0,-(D+Df)],1);
y_coor_of_left_incl_line = flip(y_coor_of_left_incl_line);
y_coor_of_symm_side = [-Df,y_coor_of_symm_side];
% disp(y_coor_of_symm_side);
y_sym_main_mesh = y_coor_of_symm_side;
%disp(y_sym_main_mesh);
y_coor_right_incl_line = [y_coor_of_left_incl_line,y_sym_main_mesh];
% disp(y_coor_right_incl_line);
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
        [x_coord_midpoint,y_coord_midpoint] =Divide_Quadri(x_coor_int(i,j),x_coor_int(i+1,j),x_coor_int(i+1,j+1),x_coor_int(i,j+1),y_coor_int(i,j),y_coor_int(i+1,j),y_coor_int(i+1,j+1),y_coor_int(i,j+1));
        x_coord_midpoint_fan_right(i,j) = x_coord_midpoint;
        y_coord_midpoint_fan_right(i,j) = y_coord_midpoint;

    end
end
[node_right_side] = getNodes_right(B,fliplr(x_coor_int), fliplr(y_coor_int), fliplr(x_coord_midpoint_fan_right), fliplr(y_coord_midpoint_fan_right));
end

%% function to get Y-coordinate of the interior of the fan mesh on right side.
function [node_left_fan_total] = GetYCoordinateLeftRight(Nb,B,Ndf,Df,x_coor_of_footing_base,x_coor_of_footing_edge_side,y_coor_of_footing_edge_side,y_coor_of_symm_side,x_coor_of_symm_side,y_coor_of_left_incl_line)
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
        [x_coord_midpoint,y_coord_midpoint] = Divide_Quadri(x_coor_of_interior(i,j),x_coor_of_interior(i+1,j),x_coor_of_interior(i+1,j+1),x_coor_of_interior(i,j+1),y_coor_interior(i,j),y_coor_interior(i+1,j),y_coor_interior(i+1,j+1),y_coor_interior(i,j+1));
        x_coord_midpoint_fan_left(i,j) = x_coord_midpoint;
        y_coord_midpoint_fan_left(i,j) = y_coord_midpoint;
    end
end
[node_left_fan_total] = getNodes(B,x_coor_of_interior, y_coor_interior, x_coord_midpoint_fan_left, y_coord_midpoint_fan_left);

end

%% function to get X-coordinate of the interior of the fan mesh on middle side.
function [node_middle_fan] = GetXCoordinateMiddle(B,Nb,Df,Lf, Nlf, y_coor_of_left_incl_line, x_coor_of_right_incl_line,  x_coor_of_footing_edge_for_below, x_coor_of_fan_base, y_coor_of_footing_edge_for_below,y_coor_of_fan_base)
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
[node_middle_fan] = getNodes(B,x_interior, y_interior, x_coord_midpoint_fan_base, y_coord_midpoint_fan_base);

end

%% Function to divide the quads into triangles and store intersecting point coordinates
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

%% function to find gp-series r
function r = findRateOfIncrease(S, a, n)
% Define the function representing the sum of the geometric series
sumGP = @(r) a * (r^n - 1) / (r - 1) - S;

% Initial guess for r
r_initial_guess = 1.1; % This can be adjusted based on the problem context

% Use fsolve to find the rate of increase r
options = optimoptions('fsolve', 'Display', 'off');
r = fsolve(sumGP, r_initial_guess, options);
end

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

%% Element equilibrium
function [A_element,B_element]= ElementEquilibrium(total_node_table, no_of_element,gamma)
    
    A_element = sparse([]);
    B_element = repmat([0;gamma],no_of_element,1);
    for i = 1:no_of_element
       
        start_row = (i-1) * 2 + 1;
        start_col = (i-1) * 9 + 1;
        
        x1 = total_node_table((i-1)*3 + 1, 1); y1 = total_node_table((i-1)*3 + 1, 2);
        x2 = total_node_table((i-1)*3 + 2, 1); y2 = total_node_table((i-1)*3 + 2, 2);
        x3 = total_node_table((i-1)*3 + 3, 1); y3 = total_node_table((i-1)*3 + 3, 2);
        
        e1 = y2 - y3; e2 = y3 - y1; e3 = y1 - y2;
        z1 = x3 - x2; z2 = x1 - x3; z3 = x2 - x1;
        
        twice_area = abs(e1 * z2 - e2 * z1);
    
        A_individual = (1 / twice_area) * [
            e1, 0, z1, e2, 0, z2, e3, 0, z3;
            0, z1, e1, 0, z2, e2, 0, z3, e3
        ];
        
        A_element(start_row:start_row+1, start_col:start_col+8) = A_individual;
    end
   
end

%% Yield condition
function [A_Yield, B_yield] = YieldEquilibrium(no_of_element, p, fi,c)
    % Total number of nodes
    fi=(pi/180)*fi;
    total_nodes = no_of_element * 3;
    num_decimal_places = 3;

    A_yield = sparse(p, 3);

    for i = 1:p
        A_i = cos(2 * pi * i / p) + sin(fi) * cos(pi / p);
        B_i = sin(fi) * cos(pi / p) - cos(2 * pi * i / p);
        C_i = 2 * sin(2 * pi * i / p);
        % to reduce precision
        A_i = round(A_i, num_decimal_places);
        B_i = round(B_i, num_decimal_places);
        C_i = round(C_i, num_decimal_places);

        A_yield(i, 1) = A_i;
        A_yield(i, 2) = B_i;
        A_yield(i, 3) = C_i;
    end
    
    D = 2 *c * cos(fi) * cos(pi / p);

    A_Yield = kron(eye(total_nodes), A_yield);
    
    B_yield = D * ones(size(A_Yield, 1), 1);
end

%% Boundary condition
function [A_bound,B_Boundary] = BoundaryCondition(left_bound_result,middle_bound_result,right_bound_result,no_of_element)
    row = size(left_bound_result,1)+size(middle_bound_result,1)+2*size(right_bound_result,1); %+size(middle_bound_result,1)
    col = 9*no_of_element;
    A_bound =sparse(row,col);
    % for right boundary theta=0

    T_r = [0, 1, 0; 0, 0, 1];
    T_l = [ 0, 0, -1];
    T_m = [0,0,1];
    B_Boundary = sparse(row,1);
    
    row_start = 1;
    for i = 1:size(left_bound_result, 1)
        index = left_bound_result(i, 1);
        
        col_start = 3 * (index - 1) + 1;
        A_bound(row_start, col_start:col_start+2) = T_l;
        B_Boundary_left(row_start,1)=0;
        row_start = row_start + 1;
    end
    b_counter=1;
    for i=1:size(middle_bound_result,1)
        index = middle_bound_result(i,1);
        col_start = 3*(index-1)+1;
        A_bound(row_start,col_start:col_start+2)=T_m;
        B_Boundary_middle(b_counter,1)=0;
        row_start=row_start+1;
        b_counter=b_counter+1;
    end
    b_counter=1;
    for i = 1:size(right_bound_result, 1)
        index = right_bound_result(i, 1);
        col_start = 3 * (index - 1) + 1;
        A_bound(row_start, col_start:col_start+2) = T_r(1,:);
        A_bound(row_start+1, col_start:col_start+2) = T_r(2,:);
        B_Boundary_right(b_counter,1)=0;
        row_start = row_start + 2;
        b_counter=b_counter+2;
    end 
    
     B_Boundary_right(b_counter-1,1)=0;
     B_Boundary=[B_Boundary_left;B_Boundary_middle;B_Boundary_right]; % ;B_Boundary_middle
    
end

%% function to label the traingle number
function labelTriangles(nodes_table)
    % Number of triangles
    numTriangles = size(nodes_table, 1) / 3;
    
    
    % Create a figure for plotting
   
    
    for i = 1:numTriangles
        % Extract the coordinates of the current triangle
        idx = (i-1)*3 + 1;
        x = nodes_table(idx:idx+2, 1);
        y = nodes_table(idx:idx+2, 2);
       
        % Calculate the centroid
        centroid_x = mean(x);
        centroid_y = mean(y);
        
        % Label the centroid with the triangle number
        text(centroid_x, centroid_y, num2str(i), 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    end
 
    
end

%% for storing the boundary edge node number and coordinates.
function [left_bound_result,right_bound_result,middle_bound_result] = findConsecutiveNodes(nodes_table, B)
    right_bound_result = [];
    left_bound_result = [];
    middle_bound_result = [];
    % Get the number of rows in nodes_table
    numNodes = size(nodes_table, 1);
    
    for i = 1:numNodes-1
        if nodes_table(i, 2) == 0 && nodes_table(i+1, 2) == 0
            if (nodes_table(i, 1) - B) > 0
                % Store the index and coordinates in the result array
                right_bound_result = [right_bound_result; i, nodes_table(i, :)];
            end
            if (nodes_table(i+1, 1) - B) > 0
                right_bound_result = [right_bound_result;i+1,nodes_table(i+1,:)];
            end
           
            % storing below footing boudnary edges
            if (nodes_table(i, 1) - B) <= 0
                % Store the index and coordinates in the result array
                middle_bound_result = [middle_bound_result; i, nodes_table(i, :)];
            end
            if (nodes_table(i+1, 1) - B) <= 0
                middle_bound_result = [middle_bound_result;i+1,nodes_table(i+1,:)];
            end
            
        end
        if nodes_table(i, 1) == 0 && nodes_table(i+1, 1) == 0
            if (nodes_table(i, 2)) <= 0
                % Store the index and coordinates in the result array
                left_bound_result = [left_bound_result; i, nodes_table(i, :)];
            end
            if (nodes_table(i+1, 2)) <= 0
                left_bound_result = [left_bound_result;i+1,nodes_table(i+1,:)];
            end
        end
    end
    % This is done to add the last boudnary edge beside footing.

    size_right_bound = size(right_bound_result,1);
    last_row= right_bound_result(size_right_bound,:);
    right_bound_result=[right_bound_result;last_row(1)+7,last_row(2),last_row(3)];
    right_bound_result=[right_bound_result;last_row(1)+9,B,0];

    %this is done to remove a duplicate values while storing the mid_bound_result
    size_midd_bound = size(middle_bound_result,1);
   % disp(size_midd_bound);
    middle_bound_result(size_midd_bound-3,:)=[];
    middle_bound_result(size_midd_bound-2,:)=[];

    %reshape the middle bound result so that nodes come in order.
% nodeNumbers = middle_bound_result(:, 1);
% reshapedNodeNumbers = reshape(nodeNumbers, 2, []).';
% reversedPairsNodeNumbers = flipud(reshapedNodeNumbers);
% reorderedNodeNumbers = reshape(reversedPairsNodeNumbers.', [], 1);
% [~, newOrder] = ismember(reorderedNodeNumbers, nodeNumbers);
% reorderedArray = middle_bound_result(newOrder, :);
% middle_bound_result = reorderedArray;


end

%% function to find discontinuity edges
function [result,new_node] = findMatchingEdges(nodes_array) %nodes_array same as total_node_table
result = [];
n = size(nodes_array, 1);
new_node=[];
index_array = (1:n)';
new_nodes_array = [index_array,nodes_array];

i=1;
while i<n-1
    triangle = new_nodes_array(i:i+2,:);
    new_node = [new_node;triangle;triangle(1,:)]; %repeats first triangle node at the end
    i=i+3;
end

n = size(new_node, 1);
for i = 1:n
    if mod(i,4)==0  %skips first node repeatition (of first element) with next element
        continue;
    end
    x1 = new_node(i, 2);
    y1 = new_node(i, 3);
    x2 = new_node(i+1, 2);
    y2 = new_node(i+1, 3);
    %
    if ~(x1 == 0 && x2 == 0) || ~(y1 == 0 && y2 == 0) %skips boundary edges
        
        for j = i+2:n-1
            x3 = new_node(j, 2);
            y3 = new_node(j, 3);
            x4 = new_node(j+1, 2);
            y4 = new_node(j+1, 3);
            if new_node(j,1) ~= new_node(i,1) && new_node(j,1) ~= new_node(i+1,1) && ...
                   new_node(j+1,1) ~= new_node(i,1) && new_node(j+1,1) ~= new_node(i+1,1) %prevent same node number for comparison
                
                if (x1 ==x3 && y1 == y3 && x2==x4 && y2==y4) || ...
                    (x1==x4 && y1==y4 && x2==x3 && y2==y3)
                    
                    dy = (y2 - y1);
                    dx = (x2 - x1);
       
                    angle_radian = atan2(dy,dx);
                    angle_degree = rad2deg(angle_radian);
                    if(abs(angle_degree) < 10^-4)
                        angle_degree = 0;
                    end
                    result = [result; new_node(i,1),new_node(j+1,1), new_node(i+1,1), new_node(j,1),angle_radian, angle_degree];
                end
            end
        end
    end
end

end

function [A_discontinuity,B_discontinuity] = DiscontinuityEquilibrium(result,no_of_element)
    no_of_discont_plane = size(result,1);
    A_discontinuity = sparse(4*no_of_discont_plane,3*no_of_element);
    B_discontinuity = sparse(4*no_of_discont_plane,1);
    for i=1:no_of_discont_plane
        theta_radian = result(i,5);
        T=[sin(theta_radian).^2, cos(theta_radian).^2, -sin(2*theta_radian); -0.5*sin(2*theta_radian), 0.5*sin(2*theta_radian), cos(2*theta_radian)];
        col1_start = (result(i,1)-1)*3+1;
        col2_start = (result(i,2)-1)*3+1;
        col3_start = (result(i,3)-1)*3+1;
        col4_start = (result(i,4)-1)*3+1;

        row = (i-1)*4+1;
        A_discontinuity(row,col1_start:col1_start+2) = T(1,:);
        A_discontinuity(row,col2_start:col2_start+2) = -T(1,:);

        A_discontinuity(row+1,col1_start:col1_start+2) = T(2,:);
        A_discontinuity(row+1,col2_start:col2_start+2) = -T(2,:);

        A_discontinuity(row+2,col3_start:col3_start+2) = T(1,:);
        A_discontinuity(row+2,col4_start:col4_start+2) = -T(1,:);

        A_discontinuity(row+3,col3_start:col3_start+2) = T(2,:);
        A_discontinuity(row+3,col4_start:col4_start+2) = -T(2,:);
    end
end

%% function to plot node numbers
function plotNodesWithNumbersInside(total_node_table)
    % Get the number of nodes
    numNodes = size(total_node_table, 1);
 
    for i = 1:3:numNodes
        % Extract the coordinates of the three nodes forming the triangle
        x_coords = total_node_table(i:i+2, 1);
        y_coords = total_node_table(i:i+2, 2);
        

        area = abs(x_coords(1)*(y_coords(2)-y_coords(3)) + x_coords(2)*(y_coords(3)-y_coords(1)) + x_coords(3)*(y_coords(1)-y_coords(2))) / 2;
        
        % Adjust the offset based on the area
        offset = 0.1 * sqrt(area);  

        % Calculate the centroid of the triangle
        centroid_x = mean(x_coords);
        centroid_y = mean(y_coords);
      
        
        % Annotate the nodes inside the triangle
       text(centroid_x - offset, centroid_y, num2str(i), 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
        text(centroid_x + offset, centroid_y, num2str(i+1), 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
        text(centroid_x, centroid_y + offset, num2str(i+2), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
    end
    
    
end

%% check discontinuity by plotting
function DrawDiscontinuity(total_node_table,result)
    hold on;
    for i=1:size(result,1)
        X=[total_node_table(result(i,2),1),total_node_table(result(i,3),1)];
        Y=[total_node_table(result(i,2),2),total_node_table(result(i,3),2)];
        plot(X,Y,'r','LineWidth',2);
    end
    
    
end

%% Objective function
function C_matrix = ObjectiveFunction(middle_bound_result,no_of_element)

    size_middle = size(middle_bound_result,1);
    C_matrix = sparse(1,9*no_of_element);
    for i=1:2:size_middle
        x1 = middle_bound_result(i,2);
        y1 = middle_bound_result(i,3);
        x2 = middle_bound_result(i+1,2);
        y2 = middle_bound_result(i+1,3);
   
        Length = sqrt((x1-x2)^2+(y1-y2)^2);
        T = (Length/2)*[0,1,0];
        index_1 = middle_bound_result(i,1);
        index_2 = middle_bound_result(i+1,1);
        col_start_1 = (index_1-1)*3+1;
        col_start_2 = (index_2-1)*3+1;
        C_matrix(1,col_start_1:col_start_1+2)=T;
        C_matrix(1,col_start_2:col_start_2+2)=T;
    end
end

function checkNodes(total_node_table)
    x_coords = total_node_table(:, 1);
    y_coords = total_node_table(:, 2);
    scatter(x_coords, y_coords, 15, 'filled'); 
    
end

function checkElementEquilibirum(total_node_table)


for i = 1:3:size(total_node_table, 1)

    X=[total_node_table(i,1),total_node_table(i+1,1),total_node_table(i+2,1)];
    Y=[total_node_table(i,2),total_node_table(i+1,2),total_node_table(i+2,2)];
    plot(X,Y,"Color","#D95319");
end


end

function DrawBoundary(left_bound_result,middle_bound_result,right_bound_result)
    bound_all = [left_bound_result;middle_bound_result;right_bound_result];
    for i=1:2:size(bound_all,1)
        X=[bound_all(i,2),bound_all(i+1,2)];
        Y=[bound_all(i,3),bound_all(i+1,3)];
        plot(X,Y,"Color","#0000FF",'LineWidth',2);
    end
end

