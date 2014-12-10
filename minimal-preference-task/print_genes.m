function [w,s1,s2,alpha,H1,H2,eps_K,phi_r,phi_l,phi0] = print_genes(genes,Nsize)

   %sensors
    s1(1,1)=bin2dec(genes(1,:))/31*16-8;
    s1(2,1)=bin2dec(genes(2,:))/31*16-8;
    s2(1,1)=bin2dec(genes(3,:))/31*16-8;
    s2(2,1)=bin2dec(genes(4,:))/31*16-8;
    
    %controller
    alpha=(bin2dec(genes(5,:))+1)/32*5; 
 
    H1=bin2dec(genes(6,:))/31*0.20; 
    H2=bin2dec(genes(7,:))/31*0.25;

    phi_r=(bin2dec(genes(8,:))+1)/32*2*pi;
    phi_l=(bin2dec(genes(9,:))+1)/32*2*pi;
    
    w=zeros(Nsize,1);
    for ii=1:Nsize
        w(ii)=bin2dec(genes(ii+9,:))/31*5;
    end
    
    phi0=zeros(Nsize,1);
    for ii=1:Nsize
        phi0(ii)=bin2dec(genes(ii+9+Nsize*1,:))/31*pi-pi/2;
    end
    
    iK=1;
    eps_K=zeros(Nsize,Nsize);
    for ii=1:Nsize
        for jj=1:Nsize
            if(ii~=jj)
                eps_K(ii,jj) = bin2dec(genes(iK+9+Nsize*2,:))/31*0.9;
                iK=iK+1;
            end
        end
    end
end
