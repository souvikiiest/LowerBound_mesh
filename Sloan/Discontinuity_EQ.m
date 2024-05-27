
syms  theta sx1 sx2 sx3 sy1 sy2 sy3 tx1 tx2 tx3 txy  sx sy txy1 txy2 txy3 txy4 sx4 sy4 ;
s= -txy*sin(2*theta) + sx*sin(theta)*sin(theta)+sy*cos(theta)*cos(theta);
t = -0.5*sx*sin(2*theta) + 0.5*sy*sin(2*theta) + txy*cos(2*theta);
 for i=1:4
    s_i(i) = subs(s, {txy, sx, sy}, {sym(['txy' num2str(i)]), sym(['sx' num2str(i)]), sym(['sy' num2str(i)])});
    t_i(i) = subs(t, {txy, sx, sy}, {sym(['txy' num2str(i)]), sym(['sx' num2str(i)]), sym(['sy' num2str(i)])});
 end
 %Disont equil eqn starts
 A = s_i(1)-s_i(2); 
 B = s_i(3)-s_i(4);
 C = t_i(1)-t_i(2);
 D = t_i(3)-t_i(4);
 %Discont equil eqn ends

 eqns = [A==0,C==0,B==0,D==0];
 var_m = [sx1;sy1;txy1;sx2;sy2;txy2;sx3;sy3;txy3;sx4;sy4;txy4;];
 [Aeq_d,b_d] = equationsToMatrix(eqns,var_m); %Aeq_d = A equlibrium matrix
 size_d = size(theta1);
 A_discont_cell = cell(1,size_d(2));
 
 for i=1:size_d(2)    
     A_discont_i = subs(Aeq_d, {theta}, theta1(i));
     A_discont_cell{i} = A_discont_i; 
 end
%disp(A_discont_cell{1});
E=(numTriangles+2);
disc_rows = 4*(E + (E+2) - 1 - 3);
Disc_mat = zeros(disc_rows, 9*E);
flag = 0;
node_mat=[1,4,3,5];
for i=1:size_d(2)
    p=node_mat(1)+3*i;
     q=node_mat(2)+3*i;
      r=node_mat(3)+3*i;
    s=node_mat(4)+3*i;

    node_mat=[node_mat,p,q,r,s];
end
%disp(node_mat);
row_num=1;
index=1;
for i=1:size_d(2)
    
    p=node_mat(index);
    q=node_mat(index+1);
    for j=1:2
        for z=1:3
            Disc_mat(row_num,3*p-3+z)=A_discont_cell{i}(j,z);
            Disc_mat(row_num,3*q-3+z)=A_discont_cell{i}(j,z+3);
            Disc_mat(row_num+2,3*p+3+z)=A_discont_cell{i}(j+2,z+6);
            Disc_mat(row_num+2,3*q+z)=A_discont_cell{i}(j+2,z+9);
        end
        row_num=row_num+1;
    end
    index=index+4;
    row_num=row_num+2;
end
%disp(A_discont_cell{1}(1,1));

size_b = size(Disc_mat);
b_disc = zeros(size_b(1),1);