clc
clearvars

%% Initial (A) and Final (B) Configurations

R_A = rot(10*pi/180,'x');
p_A = [0 0 0];

R_B = rot(20*pi/180,'x') * rot(70*pi/180,'y') * rot(40*pi/180,'z');
p_B = [3 3 3];

%% Configurations --> Unit Dual Quaternians

R_A_UQ = rot2quat(R_A);
p_A_Q = [0 p_A];
A_UDQ = [R_A_UQ, 1/2*quatProduct(p_A_Q , R_A_UQ)];

R_B_UQ = rot2quat(R_B);
p_B_Q = [0 p_B];
B_UDQ = [R_B_UQ, 1/2*quatProduct(p_B_Q , R_B_UQ)];

%% Screw Linear Interpolation (ScLERP)

Transf = qualQuatProduct(dualQuatConjucate(A_UDQ), B_UDQ);
[l, m, theta, d] = screwParameters(Transf);

tau = 0:0.1:1;
C = zeros(length(tau),8);
R = zeros(3,3,length(tau));
p = zeros(length(tau),3);

for i = 1:length(tau)
    P = [cos(tau(i)*theta/2) sin(tau(i)*theta/2)*l];
    Q = [-tau(i)*d/2*sin(tau(i)*theta/2) tau(i)*d/2*cos(tau(i)*theta/2)*l + sin(tau(i)*theta/2)*m];
    C(i,:) = qualQuatProduct(A_UDQ, [P Q]);
    [R(:,:,i), p(i,:)] = dualQuat2transformation(C(i,:));
end

%% 3D Simulation
% Cube Dimensions
length = 3;
width = 2;
height = 1;
plot3D(length,width,height,R_A,p_A,R_B,p_B,tau,R,p)


%% Functions

function R = rot(angle, axis)
if axis == 'z'
    R = [...
        cos(angle) -sin(angle) 0
        sin(angle) cos(angle)  0
        0          0           1];
    
elseif axis == 'y'
    R = [...
        cos(angle)  0  sin(angle)
        0           1 0
        -sin(angle) 0 cos(angle)];
    
elseif axis == 'x'
    R = [...
        1 0           0
        0 cos(angle) -sin(angle)
        0 sin(angle)  cos(angle)];
end
end

function Q = rot2quat(R)
q_0 = 1/2*sqrt(trace(R)+1);
q_r = 1/2*[...
    sign(R(3,2) - R(2,3))*sqrt(R(1,1) - R(2,2) - R(3,3) + 1)
    sign(R(1,3) - R(3,1))*sqrt(R(2,2) - R(3,3) - R(1,1) + 1)
    sign(R(2,1) - R(1,2))*sqrt(R(3,3) - R(1,1) - R(2,2) + 1)];
Q = [q_0 q_r.'];
end

function R = quat2rot(Q)
q_0 = Q(1);
q_1 = Q(2);
q_2 = Q(3);
q_3 = Q(4);
R(1,1) = q_0*q_0 + q_1*q_1 - q_2*q_2 - q_3*q_3;
R(1,2) = 2*(q_1*q_2 - q_0*q_3);
R(1,3) = 2*(q_0*q_2 + q_1*q_3);
R(2,1) = 2*(q_0*q_3 + q_1*q_2);
R(2,2) = q_0*q_0 - q_1*q_1 + q_2*q_2 - q_3*q_3;
R(2,3) = 2*(q_2*q_3 - q_0*q_1);
R(3,1) = 2*(q_1*q_3 - q_0*q_2);
R(3,2) = 2*(q_0*q_1 + q_2*q_3);
R(3,3) = q_0*q_0 - q_1*q_1 - q_2*q_2 + q_3*q_3;
end

function PQ = quatProduct(P, Q)
p_0 = P(1);
q_0 = Q(1);
p_r = P(2:4);
q_r = Q(2:4);
scalarPart = p_0*q_0 - dot(p_r,q_r);
vectorPart = p_0*q_r + q_0*p_r + cross(p_r,q_r);
PQ = [scalarPart vectorPart];
end

function quatConjucate = quatConjucate(Q)
q_0 = Q(1);
q_r = Q(2:4);
quatConjucate = [q_0 -q_r];
end

function dualQuatConjucate = dualQuatConjucate(A)
P = A(1:4);
Q = A(5:8);
dualQuatConjucate = [quatConjucate(P) quatConjucate(Q)];
end

function AB = qualQuatProduct(A, B)
P_A = A(1:4);
Q_A = A(5:8);
P_B = B(1:4);
Q_B = B(5:8);
AB = [quatProduct(P_A,P_B) quatProduct(P_A,Q_B)+quatProduct(Q_A,P_B)];
end

function [u, theta] = quat2AxisAngle(Q)
q_0 = Q(1);
q_r = Q(2:4);
if (norm(q_r) <= 1e-12)
    % Null or full rotation, the angle is 0 (modulo 2*pi) --> singularity: The unit vector u is indeterminate.
    % By convention, set u to the default value [0, 0, 1].
    u = [0, 0, 1];
    theta = 0;
else
    u = q_r/norm(q_r);
    theta = 2*atan2(norm(q_r), q_0);
end
end

function [R, p] = dualQuat2transformation(A)
P = A(1:4);
Q = A(5:8);
R = quat2rot(P);
p_Q = 2 * quatProduct(Q, quatConjucate(P));
p = p_Q(2:4);
end

function [l, m, theta, d] = screwParameters(A)
P = A(1:4);
Q = A(5:8);
[l, theta] = quat2AxisAngle(P);
if theta == 0 || theta == pi
    disp("Screw axis is at infinity!")
else
    [~, t] = dualQuat2transformation(A);
    d = dot(t,l);
    m = 1/2 * (cross(t,l) + (t - d * l) * cot(theta/2));
end
end

function plot3D(x,y,z,R_A,p_A,R_B,p_B,tau,R,p)
Rec_Cuboid = [...
    x/2 -x/2 -x/2  x/2 x/2  x/2  x/2  x/2  x/2 -x/2 -x/2  x/2 -x/2 -x/2 -x/2 -x/2
    y/2  y/2  y/2  y/2 y/2 -y/2 -y/2  y/2 -y/2 -y/2 -y/2 -y/2 -y/2  y/2  y/2 -y/2
    z/2  z/2 -z/2 -z/2 z/2  z/2 -z/2 -z/2 -z/2 -z/2  z/2  z/2  z/2  z/2 -z/2 -z/2];
Vertices = [-x/2 -y/2 -z/2; -x/2 y/2 -z/2; x/2 y/2 -z/2; x/2 -y/2 -z/2; -x/2 -y/2 z/2; -x/2 y/2 z/2; x/2 y/2 z/2; x/2 -y/2 z/2].';
Faces = [1 2 3 4;5 6 7 8;3 4 8 7;1 2 6 5;2 3 7 6;1 4 8 5];
figure
for i = 1:length(tau)-1
    verts_A = R_A * Vertices + p_A.';
    patch('Faces',Faces,'Vertices',verts_A.','FaceColor',[0.7 0.8 1],'EdgeColor','k','LineWidth',1);
    
    verts_B = R_B * Vertices + p_B.';
    patch('Faces',Faces,'Vertices',verts_B.','FaceColor',[0.7 0.8 1],'EdgeColor','k','LineWidth',1);
    
    hold on
    Rec_Cuboid_i = R(:,:,i)*Rec_Cuboid + p(i,:).';
    plot3(Rec_Cuboid_i(1,:),Rec_Cuboid_i(2,:),Rec_Cuboid_i(3,:),'Color',[0,0,0] + 0.6,'LineWidth',0.5)
    hold on
    plot3(Rec_Cuboid_i(1,2),Rec_Cuboid_i(2,2),Rec_Cuboid_i(3,2),'o','Color',[0,0,0] + 0.6,'LineWidth',2,'MarkerSize',2)
    
    axis equal
    set(gcf,'Color',[1 1 1]);
    axis off
    view(30,20)
    pause(0.1)
end
hold on
text(p_A(1),p_A(2),p_A(3),'(A)')
plot3(p(:,1), p(:,2),p(:,3),'k')
end
