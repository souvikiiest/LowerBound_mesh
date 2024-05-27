theta_s=3.14;
size_c = size(Global_s);
c_mat_zeros = zeros (size_c(1)-6,1);
c_mat_load = (width/2)*[(sin(theta_s))^2;(cos(theta_s))^2;-sin(2*theta_s);(sin(theta_s))^2;(cos(theta_s))^2;-sin(2*theta_s)];
c_mat=[c_mat_load;c_mat_zeros];
% Minimize -[c]'*[sigma_c_mat]
% SUbjected to [Global_A]*[Global_s]<=Global_b
f=-c_mat;
k=Global_s;
% B_matrix = zeros(129,1);
[solution, fval]=linprog(f,Yield_matrix,b_yield,Global_A,Global_B);
disp(fval);
qf = -fval/width;
if fi == 0
    Nc = qf/c;
    disp("Nc:"+Nc);
else
    % Nq = 1+(qf*tan(fi))/(c);
    % Nc = (Nq-1)*(cot(fi));
    % disp("Nc:"+Nc);
    % disp("Nq:" + Nq);
    disp(qf);
    disp(qf/c);
end

