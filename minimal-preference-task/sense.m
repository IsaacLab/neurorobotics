function [Id,Ii]=sense(L,x,y,p,R,ang)
   
   Th=100;
   sA=1;
   sB=0.03;
      
   xi=x + R*cos(p+ang);
   yi=y + R*sin(p+ang);
    
   xd=x + R*cos(p-ang);
   yd=y + R*sin(p-ang);
   
   di=sqrt((xi-L(1))^2 + (yi-L(2))^2);
      
   piz=atan((L(2)-yi)/(L(1)-xi))+pi*(L(1)-xi<0) - (p+ang);
   piz=mod(piz,2*pi);
   Ii=(1+cos(piz*2))*(cos(piz)>cos(pi/2))*0.5./(1+exp(sB*(di-Th)));
      
   dd=sqrt((xd-L(1))^2 + (yd-L(2))^2);
      
   pd=atan((L(2)-yd)/(L(1)-xd))+pi*(L(1)-xd<0) - (p-ang);
   pd=mod(pd,2*pi);
   Id=(1+cos(pd*2))*(cos(pd)>cos(pi/2))*0.5./(1+exp(sB*(dd-Th)));

end
