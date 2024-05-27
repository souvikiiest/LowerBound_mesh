%load('Mesh.mat', 'finalArray');
syms x y x1 y1 x2 y2 x3 y3 sx1 sx2 sx3 sy1 sy2 sy3 tx1 tx2 tx3 g;
P = [1, x, y];
C = [1, x1, y1; 1, x2, y2; 1, x3, y3];
%eqn = N*P/C;
N = P * inv(C); 
%Shape function (1X3)
N1 = N(1);
N2 = N(2);
N3 = N(3);
%End of Shape function
%Stress eq condition starts
Sx=N1*sx1 + N2*sx2 +N3*sx3;%(sigma_x)
Sy=N1*sy1 + N2*sy2 +N3*sy3;%(sigma_y)
Txy=N1*tx1 + N2*tx2 +N3*tx3;%(tau_xy)

dsx_dx = diff(Sx,x); %d(sigma_x)/dx;
dsy_dy = diff(Sy,y); %d(sigma_y)/dy;
dtxy_dx = diff(Txy,x);%d(Tau_xy)/dx;
dtxy_dy = diff(Txy,y);%d(Tau_xy)/dy;

Eq1 = dsx_dx + dtxy_dy; %Eq1 = 0
Eq2 = dsy_dy + dtxy_dx; % Eq2 = gamma
sigma_m = [sx1;sy1;tx1;sx2;sy2;tx2;sx3;sy3;tx3;];
eqns = [Eq1 == 0, Eq2 == g];
% Stress equilibrium matrix (A)
[A_stress,b_e] = equationsToMatrix(eqns, sigma_m);
elements = numTriangles + 2;
A_stress_cell = cell(1, elements);

for i=1:elements
        x1_val = finalArray(1,1,i);
        x2_val = finalArray(2,1,i);
        x3_val = finalArray(3,1,i);
        y1_val = 0;
        y2_val = finalArray(2,2,i);
        y3_val = finalArray(3,2,i);
        %disp(['x1:', num2str(x1_val), ', y1:', num2str(y1_val), ', x2:', num2str(x2_val), ', y2:', num2str(y2_val), ', x3:', num2str(x3_val), ', y3:', num2str(y3_val)]);
    A_stress_i = subs(A_stress, {x1, x2, x3, y1, y2, y3}, {x1_val, x2_val, x3_val, y1_val, y2_val, y3_val});
    A_stress_cell{i} = A_stress_i;
end
%disp(A_stress_cell{3});
rows = 2*(numTriangles+2);
col = 9*(numTriangles+2);
A_mat = zeros(rows,col);
flag=0;
temp=1;
for i=1:numTriangles+2
    
    for j=1:2
        flag =flag+1;
        z=temp;
        for k=1:9
            A_mat(flag,z) = A_stress_cell{i}(j,k);
            z=z+1;
        end
        

    end
    temp=9*i+1;

end
size_a = size(A_mat);
z=1;
b_e=[0;0]; %gamma = 18
for i=1:size_a(1)/2
    b_equil(z,1) = b_e(1,1);
    b_equil(z+1,1) = b_e(2,1);
    z=z+2;
end
%disp(b_equil);