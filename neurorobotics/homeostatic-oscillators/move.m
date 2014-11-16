function [x,y,p,vt]=move(x,y,p,Mi,Mr,R,T)
  
  vt=(Mr+Mi)/2;
  vr=(Mr-Mi)/2/R;
    
  xant=x;
  yant=y;
  pant=p;
    
  x= x + vt*T*cos(p);
  y= y + vt*T*sin(p);
    
  p= p + vr*T;
  % p= mod(p,2*pi);
    
  if(x==NaN)
    display([x,y,p,xant,yant,pant,Mr,Mi,R,T,vt,vr])
  end

end
