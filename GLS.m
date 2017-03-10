% Nicht lineares GLS

function F=GLS(x)

global U m r C

xv(1)=x(1);
xv(2)=x(2);
xv(3)=x(3);
n=[sind(x(5))*cosd(x(4)) sind(x(5))*sind(x(4)) cosd(x(5))];
F=zeros(1,5);

for i=1:8
    F(i)= U(i)-C*dot(Bfeld(m(i,:),xv-r(i,:)),n);
end