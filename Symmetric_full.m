clc;
clear all;

%All inputs
B = 2;
L = 20;
D = 30;
Nd = 20;
Nb = 20;

%Function call
generateMesh(B,L,D,Nd,Nb);

%%
function generateMesh (B,L,D,Nd,Nb)
    %Horizontal boundary
    %Right of symmetry
    x_Footing_base_right = linspace(0, B, Nd+1); %Dividing base of footing in Nd segments
    x_H_boundarybase_right = linspace(0, L, Nd+1); % Dividing base of horizontal boundary in Nd segments
   
    y_H_boundarybase_right = -D*ones(size(x_H_boundarybase_right));
    y_footingbase_right = zeros(size(x_Footing_base_right));
    plot([x_Footing_base_right;x_H_boundarybase_right], [y_footingbase_right; y_H_boundarybase_right], 'k-'); %Plots vertical lines
    hold on;
    
    %Left of symmetry
    x_Footing_base_left = linspace(0, -B, Nd+1); %Dividing base of footing in Nd segments
    x_H_boundarybase_left = linspace(0, -L, Nd+1); % Dividing base of horizontal boundary in Nd segments
   
    y_H_boundarybase_right = -D*ones(size(x_H_boundarybase_right));
    y_footingbase_right = zeros(size(x_Footing_base_right));
    plot([x_Footing_base_left;x_H_boundarybase_left], [y_footingbase_right; y_H_boundarybase_right], 'k-'); %Plots vertical lines

    %Vertical boundary
    %Right side
    V_boundarybase_right = linspace(0, -D, Nd+1); %This is the y-coor of the right side vertical boundary
    x_V_footingbase_right = B*ones(size(V_boundarybase_right)); %x-coor array of footing edge
    x_V_boundarybase_right = L*ones(size(V_boundarybase_right)); %x-coor array of right vertical boundary edge
    y_V_footingbase_right = zeros(size(V_boundarybase_right)); %y-coor array of footing edge.
    
    plot([x_V_footingbase_right;x_V_boundarybase_right], [y_V_footingbase_right;V_boundarybase_right], 'k-'); %Plots inclined lines
    
    %Left side
    V_boundarybase_left = linspace(0, -D, Nd+1); %This is the y-coor of the left side vertical boundary
    x_V_footingbase_left = -B*ones(size(V_boundarybase_right)); %x-coor array of footing edge
    x_V_boundarybase_left = -L*ones(size(V_boundarybase_right)); %x-coor array of left vertical boundary edge
    y_V_footingbase_left = zeros(size(V_boundarybase_right)); %y-coor array of footing edge.
    
    plot([x_V_footingbase_left;x_V_boundarybase_left], [y_V_footingbase_left;V_boundarybase_left], 'k-'); %Plots inclined lines
    
    %Vertical lines (Right side)
    start = log10(B);
    stop = log10(L);
    x_top_boundary_right = logspace(start, stop, Nb+1);
    m = D/(B-L);
    c = (-D*B)/(B-L);
    
    y_diagonal_right = x_top_boundary_right*m + c; %This is same as the symmetric side y-coordinate
    y_top_boundary_right = zeros(size(x_top_boundary_right));
    x_symmetric_boundary_right = zeros(size(x_top_boundary_right)); %For below footing space
    
    plot([x_top_boundary_right;x_top_boundary_right],[y_top_boundary_right;y_diagonal_right],'k-'); %Plots vertical lines
    plot([x_top_boundary_right;x_symmetric_boundary_right],[y_diagonal_right;y_diagonal_right],'k-'); %Plots horizontal lines below footing space
    
    %Vertical lines (Left side)
    start = log10(-B);
    stop = log10(-L);
    x_top_boundary_left = logspace(start, stop, Nb+1);
    m = D/(-B+L);
    c = (-D*(-B))/(-B+L);
    
    y_diagonal_left = x_top_boundary_left*m + c; %This is same as the symmetric side y-coordinate
    y_top_boundary_left = zeros(size(x_top_boundary_left));
    x_symmetric_boundary_left = zeros(size(x_top_boundary_left)); %For below footing space
    
    plot([x_top_boundary_left;x_top_boundary_left],[y_top_boundary_left;y_diagonal_left],'k-'); %Plots vertical lines
    plot([x_top_boundary_left;x_symmetric_boundary_left],[y_diagonal_left;y_diagonal_left],'k-'); %Plots horizontal lines below footing space
    
    grid on;
    hold off;
    
   % x-coord of intersection points below footing base (right side)
   [x_int_belowbase_right]= GetXCoordinatesBelowFootingRight(Nd, x_Footing_base_right, x_H_boundarybase_right, y_diagonal_right, D);
   
   % y-coord of intersection points right of footing base (right side)
   [y_int_right_base_right] = GetYCoordinatesRightOfFootingRight(Nd,B,L, V_boundarybase_right, x_top_boundary_right, Nb);
   
   % x-coord of intersection points below footing base (left side)
   [x_int_belowbase_left]= GetXCoordinatesBelowFootingLeft(Nd, x_Footing_base_left, x_H_boundarybase_left, y_diagonal_left, D);
    
   % y-coord of intersection points right of footing base (left side)
   [y_int_right_base_left] = GetYCoordinatesRightOfFootingLeft(Nd,B,L, V_boundarybase_left, x_top_boundary_left, Nb);
    
   disp(y_top_boundary_right);
   disp(y_top_boundary_left);
end

%%Right Side
function [x_int_belowbase_right] = GetXCoordinatesBelowFootingRight(Nd, x_Footing_base, x_H_boundarybase, y_symmetric_boundary, D)

    %Trimming first and last values of arrays as we already have those
    %coordinates

    y_symmetric_boundary([1,Nd+1])=[];
    x_Footing_base([1,Nd+1])=[];
    x_H_boundarybase([1,Nd+1])=[];
    
    for i = 1:Nd-1
        
        X = [x_Footing_base(i), x_H_boundarybase(i)];
        Y = [0, -D];

        coefficient = polyfit(X, Y, 1);

        m = coefficient(1);
        c = coefficient(2);
        x_int_belowbase_right(:,i) = (y_symmetric_boundary - c)/m; % x-coord of intersection points below footing base
        
        %Each column = Inclined line from left
        %Each row = Horizontal line from top
        %Ex: 3,2 = 3rd horizontal line intersecting with 2nd inclined line
    end
    
end

function [y_int_right_base_right] = GetYCoordinatesRightOfFootingRight(Nd,B,L,V_boundarybase,x_top_boundary,Nb)
    
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
        y_int_right_base_right(:, i) = m*x_top_boundary +c; % y-coord of intersection points right of footing base

        %Each column = one inclined line from top
        %Each row = Vertical line intersect from footing edge
        %Ex: 3,2 = 3rd vertical line intersecting with 2nd inclined line
        
    end
end


%%Left Side
function [x_int_belowbase_left] = GetXCoordinatesBelowFootingLeft(Nd, x_Footing_base, x_H_boundarybase, y_symmetric_boundary, D)

    %Trimming first and last values of arrays as we already have those
    %coordinates

    y_symmetric_boundary([1,Nd+1])=[];
    x_Footing_base([1,Nd+1])=[];
    x_H_boundarybase([1,Nd+1])=[];
    
    for i = 1:Nd-1
        
        X = [x_Footing_base(i), x_H_boundarybase(i)];
        Y = [0, -D];

        coefficient = polyfit(X, Y, 1);

        m = coefficient(1);
        c = coefficient(2);
        x_int_belowbase_left(:,i) = (y_symmetric_boundary - c)/m; % x-coord of intersection points below footing base
        
        %Each column = Inclined line from left
        %Each row = Horizontal line from top
        %Ex: 3,2 = 3rd horizontal line intersecting with 2nd inclined line
    end
    
    x_symmetric_boundary = x_symmetric_boundary';
    all =[x_Footing_base;x_int_belowbase;x_H_boundarybase];
    All_x_below_base=[x_symmetric_boundary all x_top_boundary'];
    
    y_symmetric_boundary = [0; y_symmetric_boundary'; -D];
    All_y_below_base = repmat(y_symmetric_boundary, 1, Nd+1);
    % disp(All_y_below_base);
    % disp(All_x_below_base);
    
    for i=1:Nb
        for j=1:Nd
              [x_coord,y_coord]=  Divide_Quadri (All_x_below_base(i,j),All_x_below_base(i+1,j),All_x_below_base(i+1,j+1),All_x_below_base(i,j+1),All_y_below_base(i,j),All_y_below_base(i+1,j),All_y_below_base(i+1,j+1),All_y_below_base(i,j+1));
             
        end
    end
end

function [y_int_right_base_left] = GetYCoordinatesRightOfFootingLeft(Nd,B,L,V_boundarybase,x_top_boundary,Nb)
    
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
        y_int_right_base_left(:, i) = m*x_top_boundary +c; % y-coord of intersection points right of footing base

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
    %disp(x_coord_midpoint+" "+y_coord_midpoint);
end