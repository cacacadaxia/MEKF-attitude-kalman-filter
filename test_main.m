% % ǰ�����нǶȲ����ķ�ʽ����ô�ͱȽ���Ч
clear all;
close all;
N = 30;
sf = 200;

[qTrue,gyro] = genTrig(N,sf);
qMea = genMea(qTrue);
q0 = mulQua1(qTrue(1,:),expQua([pi,0,0]));


%%
dt = 1/sf;
Nt = length(gyro);

% rotation or vector measurement
if isempty(qMea)
    meaType = 'V';
else
    meaType = 'A';
end

% default parameters
U.P0 = 100;
U.GyroAngleRW = 0.1*pi/180;
U.rvstd = 0.1;
U.vMeaStd = [0.1,0.1];

% U = parseVar(varargin,U);

% Q and R
Q = eye(3)*U.GyroAngleRW^2*dt^2;
R = eye(3)*U.rvstd^2;%%��Ҫ


% pre-allocate memory
qEst = zeros(Nt,4);
P = zeros(3,3,Nt);
if strcmp(meaType,'V')
    qMea = zeros(Nt,4);
end

% initialize
qEst(1,:) = q0;
P(:,:,1) = eye(3)*U.P0^2;

% filter iteration
for nt = 2:Nt
    % integration
    av = 0.5*(gyro(nt-1,:) + gyro(nt,:));%%����λ��
    qEst(nt,:) = mulQua1(qEst(nt-1,:),expQua(av*dt));%%Ϊ�β�����ת����ģ�
    
    % uncertainty propagation
%     ��һ��
    F1 = expRot(av*dt)';
%     �ڶ���
    S = omegaMatrix(av*dt);
    normV  = sqrt(S(1,2)^2+S(1,3)^2+S(1,3)^2);
    F2 = eye(3)+sin(normV)/normV*S(:,:)+...
            (1-cos(normV))/normV^2*S(:,:)^2;
%     ������    
%%���޵��﹫ʽһ����
    rotation_part = av*dt;
    rotation_part = rotation_part';
    %%�����ָ��ӳ��Ĵ�С
    theta = sqrt(rotation_part'*rotation_part);         % �����������ģ
    %%ʵ�ֵ�ʱ���ȼ����1/theta��Ȼ��ʹ�ó˷��滻�������еĳ�����
    a = rotation_part/theta;            % �����ʲô��
    c = cos(theta);
    s = sin(theta);
    F3 = c*eye(3)+(1-c)*a*a'+s*omegaMatrix(a);
    
    %%
    F = F1;
    P(:,:,nt) = F*P(:,:,nt-1)*F'+Q;
    % update
    K = P(:,:,nt)*(P(:,:,nt)+R)^-1;     %%Kʼ����3ά�ģ������Э�����������ص�
    dv = K * logQua1(mulQua1(invQua1(qEst(nt,:)),qMea(nt,:)))';
    qEst(nt,:) = mulQua1(qEst(nt,:),expQua(dv)');
    P(:,:,nt) = (eye(3)-K)*P(:,:,nt);
end

%% ����ĳһ�����Ǵ�������ģ����Ҳ���
l1 = sum(qEst - qTrue,2);
l2 = sum(qEst + qTrue,2);
if abs(sum(l1(end-5:end)))<0.1
   qEstMerr =  qEst - qTrue;
   disp(1)
elseif abs(sum(l2(end-5:end)))<0.1
   qEstMerr =  qEst + qTrue;
   disp(2)
else
    fprintf('fuck , whats wrong\n');
end

figure;
subplot(2,1,1)
plot(qEstMerr);
subplot(2,1,2)
plot(qMea - qTrue);

%% 
function [ q ] = mulQua1( q1, q2)
% ����������λ��Ԫ���ĳ˷�
q1 = q1';
q2 = q2';
% multiplicate
q(1,:) = q1(1,:).*q2(1,:)-q1(2,:).*q2(2,:)-q1(3,:).*q2(3,:)-q1(4,:).*q2(4,:);
q(2,:) = q1(1,:).*q2(2,:)+q1(2,:).*q2(1,:)+q1(3,:).*q2(4,:)-q1(4,:).*q2(3,:);
q(3,:) = q1(1,:).*q2(3,:)-q1(2,:).*q2(4,:)+q1(3,:).*q2(1,:)+q1(4,:).*q2(2,:);
q(4,:) = q1(1,:).*q2(4,:)+q1(2,:).*q2(3,:)-q1(3,:).*q2(2,:)+q1(4,:).*q2(1,:);
q = q';
end

function [ v, u, theta ] = logQua1( q )
% ��Ԫ���Ķ���ӳ��
% check size and unitness
q = q';
% calculate log
normQv = sqrt(sum(q(2:4,:).^2));
u = q(2:4,:)./normQv;
theta = wrapToPi(atan2(normQv,q(1,:))*2);

% identity
indi = find(q(1,:)==1);
if ~isempty(indi)
    u(:,indi) = [1;0;0];
    theta(indi) = 0;
end
% format result
v = u.*theta;
v = v';
end

function [ invq ] = invQua1( q, checku )
q = q';
invq = [q(1,:);-q(2:4,:)];
invq = invq';
end



function [ q ] = expQua( v, check0 )
% calculate the exponential map from pure quaternion to unit quaternion
% If v is a 3-by-n or n-by-3 matrix, treat v as n 3-vectors
% If v is a 4-by-n or n-by-4 matrix, treat v as n pure quaternions
% in this case, if check0==true (default), check if v is pure quaternions
% q returns n unit quaternions in the same dimension as v

if ~exist('check0','var') || isempty(check0)
    check0 = true;
end

% check dimensions and pure quaternion
if size(v,1)==3
    tran = false;
    v = v/2;
elseif size(v,2)==3
    v = v'/2;
    tran = true;
elseif size(v,1)==4
    if check0
        if ~isempty(find(v(1,:)~=0,1))
            error('v must be pure quaternions');
        end
    end
    tran = false;
    v = v(2:4,:);
elseif size(v,2)==4
    if check0
        if ~isempty(find(v(:,1)~=0,1))
            error('v must be pure quaternions');
        end
    end
    tran = true;
    v = v(:,2:4)';
else
    error('v must be of size 4-n, n-4 for pure quaternions or 3-n, n-3 for vectors')
end

% calculate exponential map
theta = sqrt(sum(v.^2));
q = [cos(theta);v./theta.*sin(theta)];

if ~isempty(find(theta==0,1))
    q(:,theta==0) = [1;0;0;0];
end

% format result
if tran
    q = q';
end

end

function [omega]=omegaMatrix(data)

% wx=data(1)*pi/180;
% wy=data(2)*pi/180;
% wz=data(3)*pi/180;
wx=data(1);
wy=data(2);
wz=data(3);

omega=[
    0,-wz,wy;
    wz,0,-wx;
    -wy,wx,0
    ];

end