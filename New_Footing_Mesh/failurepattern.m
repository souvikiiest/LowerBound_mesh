This_is_GP_mesh;

X1 = X;
Elem_X = ElemX;
Elem_Y = ElemY;

% X=X1;
% ElemX=Elem_X;
% ElemY=Elem_Y;

B = B;
ElemX_temp = [];
ElemY_temp = [];

row = 1;
for i = 1:3:size(ElemX)
    ElemX_temp(row, :) = ElemX(i:i+2, 1)';
    ElemY_temp(row, :) = ElemY(i:i+2, 1)';
    row = row + 1;
end
ElemX = ElemX_temp;
ElemY = ElemY_temp;
EX = ElemX';
EY = ElemY';

phi = fi * (pi / 180);
c = c;
for i = 1:3*size(ElemX, 1)
    sigma1(i) = X(3*i-2);
    sigma2(i) = X(3*i-1);
    sigma3(i) = X(3*i);
    ac(i) = (sigma1(i) - sigma2(i))^2 + (2 * sigma3(i))^2;
    dc(i) = (2 * c * cos(phi) - (sigma1(i) + sigma2(i)) * sin(phi))^2;
    kk(i) = ac(i) / dc(i);
end

for i = 1:3*size(ElemX, 1)
    af(i) = EX(i);
    bf(i) = EY(i);
end

% ag(:, 1) = af;
% ag(:, 2) = bf;
% ag(:, 3) = kk;
ag(:,1:2) = total_node_table(:,1:2);
ag(:,3) = kk;
x = ag(:, 1);
Y = ag(:, 2);
z = ag(:, 3);

% x = x / B;
% Y = Y / B;

x = round(x * 1e8) / 1e8;
Y = round(Y * 1e8) / 1e8;
xmin = min(x);
ymin = min(Y);
xmax = max(x);
ymax = max(Y);
xres = 500;
yres = 500;
figure;
xv = linspace(xmin, xmax, xres);
yv = linspace(ymin, ymax, yres);
[Xinterp, Yinterp] = meshgrid(xv, yv);
Zinterp = griddata(x, Y, z, Xinterp, Yinterp);

surface(Xinterp, Yinterp, Zinterp);
hold on;
grid off
axis equal
shading interp

% nc = 64;
% temp = 0.96;
% for i = 1:nc
%     ncol(i, 1:3) = temp;
%     temp = temp - 1/nc;
%     if temp <= 0
%         temp = 0;
%     end
% end

axis equal
map = colormap("jet");
colormap(map(5:end, :));
hcol = colorbar;
cpos = get(hcol, 'Position');
cpos(4) = cpos(4) / 2;
cpos(3) = cpos(3) / 2;
cpos(1) = cpos(1) + 0.105;
cpos(2) = cpos(2) + 0.075;
set(hcol, 'Position', cpos)
late = max(x);
clear x y;
hold on
xlabel('x/B', 'Fontweight', 'bold', 'Fontsize', 15);
ylabel('y/B', 'Fontweight', 'bold', 'Fontsize', 15);
set(gca, 'Fontsize', 16);
axis([xmin - 0.15 xmax ymin - 0.15 ymax])
set(gcf, 'PaperType', 'A4')
box off;
%colormap parula
print -r600 -dtiff hb2.5as1ab1;
