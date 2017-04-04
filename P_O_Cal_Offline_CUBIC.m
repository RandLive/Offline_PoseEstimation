clear;close;clc;

global Anz_Punkte Anz_Spulen ...
    m C r funktionscount zeit U Ia ...
    max_abstand ...
    StartX StartY StartZ StartPhi StartTheta

cd DataSave
cd Cube
% % % -80dBW
load Um_Noise_80dBW
load POR
% 
% SigLength = max(size(Um_NoNoise));
% SigNoiseLevel = -60; 
% 
% 
% Um = zeros(SigLength,8);
% 
% 
% for i=1:SigLength
% Um(i,:)=Um_NoNoise(i,:)+wgn(1,8,SigNoiseLevel);
% end

Um = Um_Noise_80dBW;
% POR=POR;

% % -40dBW
% load Um_Noise_40dBW
% load POR_40dBW
% Um=Um_Noise_40dBW;
% POR=POR_40dBW;

cd ..
cd ..

StartX = 0.1;
StartY = 0.1;
StartZ = 0.1;
StartPhi = 45;
StartTheta = 45;

f = 1000;
Ia = 0.660156294774090;
I_ = Ia;
Verstaerkung = 1.443984589506881e+02;
C = - 0.010268772800074;
N_sende = 207;
% R_sende = 0.02;
R_sende = 0.015;

funktionscount=2000;                    % stop option for fsolve functon
Anz_Punkte=1;                           % Number of sensor coils
Anz_Spulen=8;                           % Number of transmitter coils
zeit=3;                               % time limitation for fsolve
max_abstand=0.2;

K_ =pi*N_sende*I_*R_sende^2;            % pre-factor for magnetic dipole
m=zeros(Anz_Spulen,3);

faktor= [1 1 1 1 1 1 1 1]';            % Without Calibration
% % magnetisches Dipol
phi=-90;    theta=90;       m(1,:)=faktor(1)*K_*[sind(theta)*cosd(phi)      sind(theta)*sind(phi)        cosd(theta)];
phi=30;     theta=60;       m(2,:)=faktor(2)*K_*[sind(theta)*cosd(phi)      sind(theta)*sind(phi)        cosd(theta)];
phi=90;     theta=30;    	m(3,:)=faktor(3)*K_*[sind(theta)*cosd(phi)      sind(theta)*sind(phi)        cosd(theta)];
phi=120;    theta=60;   	m(4,:)=faktor(4)*K_*[sind(theta)*cosd(phi)      sind(theta)*sind(phi)        cosd(theta)];
phi=120;    theta=90;       m(5,:)=faktor(5)*K_*[sind(theta)*cosd(phi)      sind(theta)*sind(phi)        cosd(theta)];
phi=240;	theta=30;    	m(6,:)=faktor(6)*K_*[sind(theta)*cosd(phi)      sind(theta)*sind(phi)        cosd(theta)];
phi=0;      theta=0;     	m(7,:)=faktor(7)*K_*[sind(theta)*cosd(phi)      sind(theta)*sind(phi)        cosd(theta)];
phi=0;      theta=90 ;      m(8,:)=faktor(8)*K_*[sind(theta)*cosd(phi)      sind(theta)*sind(phi)        cosd(theta)];
% Ortskoordinaten der Spulen
r(1,:)=[-0.05 0.10 0];
r(2,:)=[-0.15 0 0];
r(3,:)=[-0.05 -0.10 0];
r(4,:)=[0.05 -0.10 0];
r(5,:)=[0.15 0 0];
r(6,:)=[0.05 0.10 0];
r(7,:)=[0.05 0 0];
r(8,:)=[-0.05 0 0];

length_U = max(size(Um));
POM = zeros(length_U,5);
POM_KF = zeros(length_U,5);

mode = input(' 1-no update \n 2-update with previous value \n 3-update with filtered value \n 4-update with ground truth \n Enter the working mode:');

for i=1:length_U
    
    switch mode
        case 1 % No Update
            StartX = 0;
            StartY = 0;
            StartZ = 0.1;
            StartPhi = 45;
            StartTheta = 45;
            U=Um(i,:);
            [x,y,z,phi,theta]=POS(U,StartX,StartY,StartZ,StartPhi,StartTheta);
            POM(i,:)=[x,y,z,phi,theta];
            disp(i);
            pause(0);
            if i>=15
                POM_KF(i,:) = fun_KF_CV(POM(i-14:i,:))';
            else
                POM_KF(i,:) = POM(i,:);
            end
            
        case 2 % Update with last Value
            [x,y,z,phi,theta]=POS(U,StartX,StartY,StartZ,StartPhi,StartTheta);
            POM(i,:)=[x,y,z,phi,theta];
            disp(i);
            pause(0);
            StartX = x;
            StartY = y;
            StartZ = z;
            StartPhi = phi;
            StartTheta = theta;
            U=Um(i,:);
            if i>=15
                POM_KF(i,:) = fun_KF_CV(POM(i-14:i,:))';
            else
                POM_KF(i,:) = POM(i,:);
            end
         
        case 3 % Update with filtered Value
            U=Um(i,:);
            [x,y,z,phi,theta]=POS(U,StartX,StartY,StartZ,StartPhi,StartTheta);
            POM(i,:)=[x,y,z,phi,theta];
            if i>=15
                POM_KF(i,:) = fun_KF_CV(POM(i-14:i,:))';
            else
                POM_KF(i,:) = POM(i,:);
            end
            disp(i);
            pause(0);
            StartX = POM_KF(i,1);
            StartY = POM_KF(i,2);
            StartZ = POM_KF(i,3);
            StartPhi = POM_KF(i,4);
            StartTheta = POM_KF(i,5);
            if i>=15
                POM_KF(i,:) = fun_KF_CV(POM(i-14:i,:))';
            else
                POM_KF(i,:) = POM(i,:);
            end
            
        case 4 % Update with Real Value
            StartX = POR(i,1);
            StartY = POR(i,2);
            StartZ = POR(i,3);
            StartPhi = POR(i,4);
            StartTheta = POR(i,5);
            U=Um(i,:);
            [x,y,z,phi,theta]=POS(U,StartX,StartY,StartZ,StartPhi,StartTheta);
            POM(i,:)=[x,y,z,phi,theta];
            disp(i);
            pause(0);
            if i>=15
                POM_KF(i,:) = fun_KF_CV(POM(i-14:i,:))';
            else
                POM_KF(i,:) = POM(i,:);
            end
    end
    
end

POM = [1000*POM(:,1:3), POM(:,4:5)];
POM_KF = [1000*POM_KF(:,1:3), POM_KF(:,4:5)];
POR = [1000*POR(:,1:3), POR(:,4:5)];


scatter3(POM(:,1),POM(:,2),POM(:,3),'ro');hold on;scatter3(POM_KF(:,1),POM_KF(:,2),POM_KF(:,3),'b*');hold on;scatter3(POR(:,1),POR(:,2),POR(:,3),'g.');legend('Measued','Filtered','Truth');
% axis(1000*[-.02 .12 -.01 .01 0 .2]);

all_P_Error = POM(:,1:3)-POR(:,1:3);
rms_P_Error = rms(sqrt(sum((all_P_Error).^2)));

all_O_Error = POM(:,4:5)-POR(:,4:5);
rms_O_Error = rms(abs(all_O_Error(:,1))+abs(all_O_Error(:,2)));

all_P_Error_KF = POM_KF(:,1:3)-POR(:,1:3);
rms_P_Error_KF = rms(sqrt(sum((all_P_Error_KF).^2)));

all_O_Error_KF = POM_KF(:,4:5)-POR(:,4:5);
rms_O_Error_KF = rms(abs(all_O_Error_KF(:,1))+abs(all_O_Error_KF(:,2)));

fprintf('The Position RMSE is = %f mm \n',rms_P_Error);
fprintf('The Orientation RMSE is = %f degree \n',rms_O_Error);
fprintf('The Position RMSE (After KF) is = %f mm \n',rms_P_Error_KF);
fprintf('The Orientation RMSE (After KF) is = %f degree \n',rms_O_Error_KF);


