

numTriangles = 10; % Number of triangles
B = 1; % Width of foting
width=B;
phi= 10;
b=0.1;
apex = [B, 0];
angle_m = (45-phi/2)*(pi/180);
angle_p = (45 + phi/2)*(pi/180);
angleIncrement = pi /(2* numTriangles);
theta1{numTriangles+1} = zeros;

triangleVertices = zeros(3, 2, numTriangles); % 3-D array
r0 = (B)/cos(angle_p);
% Active zone
plot([B,0,0,B],[0,0,-B*tan(angle_p),0],'-b', 'LineWidth', 2);
hold on;

%Radial shear zone
for i = 1:numTriangles+1
    angle = (i-1)*angleIncrement;

    r = r0 * exp(b*(angle));
    vertex1 = [apex(1)-r*cos(angle_p+angle),apex(2)-r*sin(angle_p+angle)];
    vertex2 = [apex(1) - r*cos(angle_p+angle+angleIncrement),apex(2) - r*sin(angle_p+angle+angleIncrement)];
    
    theta1{i} = angle_p+angle;

    triangleVertices(:,:,i) = [apex;vertex1;vertex2];
    plot(triangleVertices(1:2, 1, i), triangleVertices(1:2, 2, i), '-b', 'LineWidth', 2);
    
    if i > 1 
    plot([triangleVertices(2,1,i),triangleVertices(2,1,i-1)],[triangleVertices(2,2,i),triangleVertices(2,2,i-1)],'-b', 'LineWidth', 2);
    end
   
    hold on;
    
end
%footing
x_d = B+2*r*cos(angle_m);
plot([triangleVertices(2,1,i),apex(1),x_d,triangleVertices(2,1,i)],[triangleVertices(2,2,i),apex(2),0,triangleVertices(2,2,i)],'-b', 'LineWidth', 2);
 x=[0,B,B,0];
 y=[0,0,1,1];
 fill(x,y,'r');

 %passive zone
finalArray=zeros(3,2,numTriangles+3);
finalArray(:,:,2:end-1) = triangleVertices;
finalArray(:,:,1) = [B,0;0,0;0,0];
finalArray(:,:,end)=[B,0;x_d,0;0,0];
save('data.mat','theta1');
%theta1{numTriangles+2} = angle_p+angle+angle_m;
 % for i=1:numTriangles+1
 %     theta1{i}=theta1{i}+3.14;
 % end

for i=1:numTriangles+2
    finalArray(3,:,i) = finalArray(2,:,i+1);
end

axis equal;
hold off;
