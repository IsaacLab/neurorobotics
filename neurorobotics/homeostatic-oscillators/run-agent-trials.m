clear all
close all

Nrep=1;
Trials=1000;

noise=0;
sigma=1;
sensory_inversion=0;

plasticity=1;


for rep=1:Nrep
    

mode1={'P','NP'};
mode2={'situated','passively_coupled','decoupled'};

task=0;

for iM1=1
%     for iM2=1:2

        
%         name=['dfa_',mode1{iM1},'_',mode2{iM2}]
if noise ==1
        name=['dfa1_noise_',mode1{iM1}]
else
        name=['dfa1_',mode1{iM1}]
end
        
        
        
%         decoupled=iM2==3;
%         passively_coupled=iM2==2;
%         
        

save_on=0;

C=zeros(1,Nrep);
IX=zeros(1,Nrep);


load network
% 

Ntrials=Trials;

% N=1250;
% TimeDuration=N*T;
% T=1/Fc;
% N=TimeDuration/T;
% t=0:T:N*T-T;


genes=individuos{1};

[w,s1,s2,alpha,H1,H2,eps_K,phi_r,phi_l,phi0] = print_genes(genes,Nsize);
        
alpha=alpha

xx=[];
xx1=[];

xx_pas=[];
xx1_pas=[];


       
dist=[];
RT=[];
            
            %Inicializa Osciladores

            K=zeros(Nsize,Nsize,1);
            dK=zeros(Nsize,Nsize,1);
            dK=(eps_K~=0).*rand(Nsize,Nsize)*4*pi;
            dK1=zeros(Nsize,Nsize,1);
            
            
            if plasticity==0
               refK=randi(100,1,1);
               load(['record_1/record_trial_',num2str(refK)],'K') 
               Kf=K(:,:,N/2);
            end
            
            K_pas=zeros(Nsize,Nsize,1);
            dK_pas=zeros(Nsize,Nsize,1);
            dK_pas=(eps_K~=0).*rand(Nsize,Nsize)*4*pi;
            dK1_pas=zeros(Nsize,Nsize,1);

            z=rand(Nsize,N*Ntrials+1)*2*pi;
            phi13=zeros(1,1);
            phi23=zeros(1,1);
            z1=zeros(Nsize,1);
            phi=zeros(Nsize,N*Ntrials); 
            z_h=zeros(Nsize,1); 
            
 
            z_pas=rand(Nsize,N*Ntrials+1)*2*pi;
            phi13_pas=zeros(1,1);
            phi23_pas=zeros(1,1);
            z1_pas=zeros(Nsize,1);
            phi_pas=zeros(Nsize,N*Ntrials); 
            z_h_pas=zeros(Nsize,1);       
            
            for trial=1:1
%             clc
%             display(name)
%             display([rep,trial/Ntrials])        

            if trial~=1
                z(:,1) = z(:,end);
                dK(:,:,1) = dK(:,:,end);
            end

                %Inicializa Robot
                x=zeros(1,1);
                y=zeros(1,1);
%                 p=zeros(1,1);
                p=rand(1,1)*2*pi;

                Mr=zeros(1,1);
                Ml=zeros(1,1);
                Mr_pas=zeros(1,1);
                Ml_pas=zeros(1,1);
                vt=zeros(1,1);

                Ir1=zeros(1,1);
                Il1=zeros(1,1);

                Ir2=zeros(1,1);
                Il2=zeros(1,1);
                
                D=zeros(1,1);
                
%                 L1=(rand(1,2)-0.5)*4;
%                 offset=rand(1,2);
%                 offset=offset./repmat(sum(offset,2),1,2)*4;
%                 L1=L1+offset.*sign(L1);

                r1=rand(1,1)*2*pi;
                D1=100 + 50*rand(1,1);
                L1=D1*[cos(r1),sin(r1)];

                r2= r1 + pi/2 + rand(1,1)*pi;
                D2=100 + 50*rand(1,1);
                L2=D2*[cos(r2),sin(r2)];
              
                i0=1;

                for i=1:(N*Ntrials)
                    
                    i0=i0+1;
%                     
%                 if mod(i,N*Ntrials/1000)==0
% %                     clc
%                     display([i/N/Ntrials])
%                 end
% %                 

                  if sensory_inversion==1
                    [Il1,Ir1]=sense(L1,x,y,p,R,S_ang);
                    [Il2,Ir2]=sense(L2,x,y,p,R,S_ang);
                  elseif noise ==1 
                    Ir1=rand(1,1)*sigma;
                    Ir2=rand(1,1)*sigma;
                    Il1=rand(1,1)*sigma;
                    Il2=rand(1,1)*sigma;
%                     
%                     Ir1=sin(i*T*2*pi*1/50);
%                     Ir2=sin(i*T*2*pi*1/50+pi/3);
%                     Il1=sin(i*T*2*pi*1/50+pi);
%                     Il2=sin(i*T*2*pi*1/50+pi+pi/3);               
%                     
                  else
                    [Ir1,Il1]=sense(L1,x,y,p,R,S_ang);
                    [Ir2,Il2]=sense(L2,x,y,p,R,S_ang);
                   end
                   
                    Ilr=[[Ir1,Il1]*s1;[Ir2,Il2]*s2;zeros(Nsize-2,1)];

                 
                    switch task
                        case 1
                            I= Ilr.*[1;0;zeros(Nsize-2,1)] +nn*randn(Nsize,1);
                        case 2
                            I= Ilr.*[0;1;zeros(Nsize-2,1)] +nn*randn(Nsize,1);
                        case 3
                            I= Ilr.*[1;(rand(1,1)<0.15);zeros(Nsize-2,1)] +nn*randn(Nsize,1);
                        case 4
                            I= Ilr.*[(rand(1,1)<0.15);1;zeros(Nsize-2,1)] +nn*randn(Nsize,1);
                        case 0
                            I=Ilr;
                    end 
                    
%                     if decoupled==1
%                         if i==1 
%                             I0=[rand(1,2)*s1;rand(1,2)*s2;zeros(Nsize-2,1)];
%                             I=I0;   
%                         end
%                     end
  
                    zz=repmat(z(:,i),1,Nsize);
                    zz=zz-zz';
                    zz=(mod(zz+pi,2*pi)-pi); 
                    
                    zz_pas=repmat(z_pas(:,i),1,Nsize);
                    zz_pas=zz_pas-zz_pas';
                    zz_pas=(mod(zz_pas+pi,2*pi)-pi); 

                    if plasticity == 1
                        K=(1-cos(dK))/2*alpha.*(mod(dK,4*pi)<2*pi);
                    else
                        K=Kf;
                    end
                    K=(1-cos(dK))/2*alpha.*(mod(dK,4*pi)<2*pi);
                    K_pas=(1-cos(dK_pas))/2*alpha.*(mod(dK_pas,4*pi)<2*pi);

                    
                    z1 = w + I + diag(K*sin(zz));

                    phi(:,i) = angle(K*exp(1j*z(:,i))) - z(:,i);
                    phi(:,i) = (mod(phi(:,i)+pi,2*pi)-pi); 
                    

                    z1_pas = w + I + diag(K_pas*sin(zz_pas));

                    phi_pas(:,i) = angle(K_pas*exp(1j*z_pas(:,i))) - z_pas(:,i);
                    phi_pas(:,i) = (mod(phi_pas(:,i)+pi,2*pi)-pi); 
                    

                    p_phi = eps_K.*repmat(homeostatic_regulator((phi(:,i)-phi0),H1,H2),1,Nsize);

                    diff_z= zz' - repmat(phi0,1,Nsize);
                    
                    p_phi_pas = eps_K.*repmat(homeostatic_regulator((phi_pas(:,i)-phi0),H1,H2),1,Nsize);

                    diff_z_pas= zz_pas' - repmat(phi0,1,Nsize);



                    if plasticity == 1
                        dK1 = (p_phi).*(1-cos(diff_z))/2;
                        dK1_pas = (p_phi_pas).*(1-cos(diff_z_pas))/2;
                    else
                        dK1 = zeros(Nsize);
                        dK1_pas = zeros(Nsize);
                    end


                    z(:,i+1) = z(:,i) + z1*T;
                    dK = dK + dK1*T;
                    
                    z_pas(:,i+1) = z_pas(:,i) + z1_pas*T;
                    dK_pas = dK_pas + dK1_pas*T;

                    
                    Mr = (sin(phi(min([Nsize,4]),i)-phi_r))*Mov;
                    Ml = (sin(phi(3,i)-phi_l))*Mov;
                    
                    Mr_pas = (sin(phi_pas(min([Nsize,4]),i)-phi_r))*Mov;
                    Ml_pas = (sin(phi_pas(3,i)-phi_l))*Mov;

                    [x,y,p]=move(x,y,p,Ml,Mr,R,T);

                    vt=(Mr+Ml)/2;
                    
                    d_min=min([sqrt((x-L1(1)).^2 + (y-L1(2)).^2);sqrt((x-L2(1)).^2 + (y-L2(2)).^2)],[],1);
                
                    
            
                    
                    if d_min<25 || i0>(1250*10)
%                         display(d_min)

                        display([name,'_',num2str(rep)])
                        display([i/(N*Ntrials)])
                        display([sqrt((x-L1(1)).^2 + (y-L1(2)).^2),sqrt((x-L2(1)).^2 + (y-L2(2)).^2)])
                        
                        dist=[dist,[sqrt((x(end)-L1(1)).^2 + (y(end)-L1(2)).^2);sqrt((x(end)-L2(1)).^2 + (y(end)-L2(2)).^2)]];
                     
                        
                        x=0;
                        y=0;
                        p=rand(1,1)*2*pi;
                        
                        r1=rand(1,1)*2*pi;
                        D1=100 + 50*rand(1,1);
                        L1=D1*[cos(r1),sin(r1)];

                        r2= r1 + pi/2 + rand(1,1)*pi;
                        D2=100 + 50*rand(1,1);
                        L2=D2*[cos(r2),sin(r2)];
                        
                        RT=[RT,i0*T];
                        i0=0;


                    end

                end
                
                

                
              
                xx=[xx,cos(phi)];
                xx_pas=[xx_pas,cos(phi_pas)];
           
                
                xx1=[xx1,cos(z)];
                xx1_pas=[xx1_pas,cos(z_pas)];
              

                if save_on ==1
                    savefile=['record_trial_',num2str(trial),'.mat'];
                    save(savefile);
                end
            
            end

figure
plot(dist')


figure
plot(xx')

% H=hilbert(mean(xx1,1));
% x=abs(H);
% 
x=mean(xx,1);
% x=xx(3,:);


save([name,'_',num2str(rep)],'xx1','xx1_pas','RT','dist','T')


    end
% end

end

