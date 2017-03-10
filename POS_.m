function [x,y,z,phi,theta]=POS_(U,x,y,z,phi,theta)

[x,y,z,phi,theta,flag]=POS_fast(U,x,y,z,phi,theta);
if flag==0
    StartPhi = 0;
    [x,y,z,phi,theta,flag]=POS_fast(U,x,y,z,StartPhi,theta);
end
if flag==0
    StartPhi = 90;
    [x,y,z,phi,theta,flag]=POS_fast(U,x,y,z,StartPhi,theta);
end
if flag==0
    StartPhi = 180;
    [x,y,z,phi,theta,flag]=POS_fast(U,x,y,z,StartPhi,theta);
end
if flag==0
    StartPhi = 270;
    [x,y,z,phi,theta,flag]=POS_fast(U,x,y,z,StartPhi,theta);
end
if flag==0
    StartPhi = 360;
    [x,y,z,phi,theta,flag]=POS_fast(U,x,y,z,StartPhi,theta);
end
if flag==0
    disp('5')
    [x,y,z,phi,theta,flag]=POS_fast(U,x,y,z,StartPhi,theta);
end

