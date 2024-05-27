syms sx1 sy1 txy1 sx2 sy2 txy2 sx sy txy q1 t1 q2 t2 q3 t3 theta
s= -txy*sin(2*theta) + sx*sin(theta)*sin(theta)+sy*cos(theta)*cos(theta);
t = -0.5*sx*sin(2*theta) + 0.5*sy*sin(2*theta) + txy*cos(2*theta);
 for i=1:4
    s_i(i) = subs(s, {txy, sx, sy}, {sym(['txy' num2str(i)]), sym(['sx' num2str(i)]), sym(['sy' num2str(i)])});
    t_i(i) = subs(t, {txy, sx, sy}, {sym(['txy' num2str(i)]), sym(['sx' num2str(i)]), sym(['sy' num2str(i)])});
 end

A = s_i(1);
B = t_i(1);
C = s_i(2);
D = t_i(2);

eqns = [A == q1, B == t1, C == q2, D == t2];
var_m = [sx1;sy1;txy1;sx2;sy2;txy2];

[A_bound, b_bound_temp] = equationsToMatrix(eqns,var_m);
theta_bound = [0,pi/2,pi];
elements = numTriangles + 2;
A_bound_cell = cell(1, 3);
for i = 1:3
        theta2 = theta_bound(i);
        A_bound_i = subs(A_bound, {theta}, theta2);
        A_bound_cell{i} = A_bound_i;    
end    


%disp(A_bound_cell{1});
%disp(A_bound_cell{2});
%disp(A_bound_cell{3});
col_bound=9*(numTriangles+2);
node_mat=[1,2,2,3,(numTriangles+2)*3,(numTriangles+2)*3-2];
Bound_mat = zeros(12,col_bound);
%Loaded_Edge
Bound_mat(1,2)=1;
Bound_mat(2,3)=1;
Bound_mat(3,5)=1;
Bound_mat(4,6)=1;
%Non_loaded_edges
Bound_mat(5,4)=1;
Bound_mat(6,6)=-1;
Bound_mat(7,7)=1;
Bound_mat(8,9)=-1;
node_col = 3*(3*(numTriangles+2)-1)+1;
Bound_mat(11,node_col+1)=1;
Bound_mat(12,node_col+2)=1;
Bound_mat(9,node_col-5)=1;
Bound_mat(10,node_col-4)=1;


%b_bound_const=[q1;0;q2;0;q2;0;q3;0;0;0;0;0];

%disp (b_bound);
b_bound = [0;0;0;0;0;0;0;0];
Bound_mat([1, 3, 5, 7], :)=[];

