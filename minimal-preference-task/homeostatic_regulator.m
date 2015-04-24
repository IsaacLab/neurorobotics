function y=homeostatic_regulator(x,H1,H2)

    x=mod(x+pi,2*pi)-pi;
    
    y=(abs(x)-H1*pi)/H2/pi;
    y(abs(x)<=(H1*pi))=0;
    y(abs(x)>(H1*pi+H2*pi))=1;

end 
