clc;
clear all;

% B = 5;
% p=24;
% numTriangles = 6;

Mesh
Stress_Equlibrium
Yield_Condition
Discontinuity_EQ
Boundary_cond

Global_A=[A_mat;Bound_mat;Disc_mat];
Global_A=double(Global_A);

%disp(Global_A);
Global_B=[b_equil;b_bound;b_disc];
Global_B=double(Global_B);
%disp (Global_B);
Yield_matrix=double(Yield_matrix);
%sigma matrix
syms sx sy txy;
Global_s=sym([]);
for i=1:3*E
    Global_s=[Global_s;sym(['sx' num2str(i)]);sym(['sy' num2str(i)]);sym(['txy' num2str(i)])];
end
%disp(s_mat);
c_matrix
