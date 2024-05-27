
syms A sx B sy C txy D c k sx1 sx2 sx3 sy1 sy2 sy3 txy1 txy2 txy3 c1 c2 c3 
fi= phi*pi/180;
p=3;
s = (cos(2*pi*k/p)+sin(fi)*cos(pi/p))*sx;
s2 = (sin(fi)*cos(pi/p)-cos(2*pi*k/p))*sy;
s3 = 2*sin(2*pi*k/p)*txy; 
s4 = 2*c*cos(fi)*cos(pi/p);
Yield_arr = sym(zeros(3, 4, 3*numTriangles));

flag=1;
for z = 1:numTriangles*3
    s_yield1 = subs(s, k, z);
    s_yield2 = subs(s2, k, z);
    s_yield3 = subs(s3, k, z);
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
temp=1;
Yield_2d = sym(zeros(numTriangles*9,4));
for i=1:numTriangles*3
    for j=1:3
        Yield_2d(temp,:)=Yield_arr(j,:,i);
        temp=temp+1;
    end
end

coefficients = sym(zeros(9*numTriangles, 3));
for i=1:3:9*numTriangles
    for j=1:3
        z=j-1;
     
        coefficients(i+z,1)=coeffs(Yield_2d(i+z,1), sym(['sx' num2str(j)]));
        coefficients(i+z,2)=coeffs(Yield_2d(i+z,2), sym(['sy' num2str(j)]));
        coefficients(i+z,3)=coeffs(Yield_2d(i+z,3), sym(['txy' num2str(j)]));
        
    end
end
% fi=fi;
% newcoef=subs(coefficients,fi,fi)
coeff_simpl = simplify(coefficients);
A_yield = vpa(coeff_simpl,4);
disp(A_yield);
