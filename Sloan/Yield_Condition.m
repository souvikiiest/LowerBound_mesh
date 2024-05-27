
syms A sx B sy txy C D c k sx1 sx2 sx3 sy1 sy2 sy3 txy1 txy2 txy3 c1 c2 c3 c

fi= phi*pi/180;
p=24;
s = (cos(2*pi*k/p)+sin(fi)*cos(pi/p))*sx;
s2 = (sin(fi)*cos(pi/p)-cos(2*pi*k/p))*sy;
s3 = 2*sin(2*pi*k/p)*txy; 
s4 = 2*c*cos(fi)*cos(pi/p);
Yield_arr = sym(zeros(3, 4, 3*numTriangles));
s_yield3 = subs(s3, k, 3);

flag=1;
for z = 1:p
    s_yield1 = s; %subs(s, k, z);
    s_yield2 = s2;%subs(s2, k, z);
    s_yield3 = s3;%subs(s3, k, z);
    s_yield4 = s4;
    for i = 1:3
        s_yield1_1 = subs(s_yield1, {sx}, {sym(['sx' num2str(i)])});
        Yield_arr(i,1, z) = s_yield1_1;
        
        s_yield2_1 = subs(s_yield2, {sy}, {sym(['sy' num2str(i)])});
        Yield_arr(i,2, z) = s_yield2_1;
        
            s_yield3_1 = subs(s_yield3, {txy}, {sym(['txy' num2str(i)])});
            Yield_arr(i,3, z) = s_yield3_1;

        s_yield4_1 = subs(s_yield4, {c}, {sym(['c' num2str(i)])});
        Yield_arr(i,4, z) = s_yield4_1;
        flag=flag+1;
    end
    
end
%disp(Yield_arr);
temp=1;
Yield_2d = sym(zeros(p*3,4));
for i=1:p
    for j=1:3
        Yield_2d(temp,:)=Yield_arr(j,:,i);
        temp=temp+1;
    end
end
%disp(Yield_2d);

coefficients = sym(zeros(p, 3));
for i=1:3:p
    for j=1:3
        z=j-1;

        coefficients(i+z,1)=coeffs(Yield_2d(i+z,1), sym(['sx' num2str(j)]));
        coefficients(i+z,2)=coeffs(Yield_2d(i+z,2), sym(['sy' num2str(j)]));
        coefficients(i+z,3)=coeffs(Yield_2d(i+z,3), sym(['txy' num2str(j)]));

    end
end
%disp(coefficients);
i = 1;
coefficients_val = sym(zeros(p, 3));
    for j=1:p
        z=j-1;
        coefficients_val(j,1)=subs(coefficients(j,1),k,j);
        coefficients_val(j,2)=subs(coefficients(j,2),k,j);
        coefficients_val(j,3)=subs(coefficients(j,3),k,j);

        %disp(i+z);
    end
    
coeff_simpl = simplify(coefficients_val);
A_yield = vpa(coeff_simpl,4);

%Yield_matrix = zeros(3*p*(numTriangles+2),9*(numTriangles+2));
E=(numTriangles+2);
Yield_matrix = kron(eye(E*3), A_yield);
%disp(Yield_matrix);
%disp(coefficients);
c= 0.0001;
b_yield_0 = 2*c*cos(fi)*cos(pi/p);
for i = 1:size(Yield_matrix)
    b_yield(i,1)=b_yield_0;
end    

%disp(b_yield);
