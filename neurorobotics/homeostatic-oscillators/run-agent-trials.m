clear all
close all

Nrep=1;
Trials=1000;

% Choose in the agent is feeded with noise
noise=0; % 0 - siguated agent; 1 - agent feeded with noise
sigma=1; 

% Choose if the agent has synaptic plasticity
plasticity=1; % 0 - frozen weights agent; 1 - plastic agents

mode1={'P','NP'};
mode2={'situated','passively_coupled','decoupled'};

task=0; %task 0 for two ligths
save_on=0; % save_on=1 for saving results into file

if plasticity==1
    iM1=1
else
    iM1=2
    
for rep=1:Nrep

    if noise ==1
            name=['dfa1_noise_',mode1{iM1}]
    else
            name=['dfa1_',mode1{iM1}]
    end

    load network
    
    Ntrials=Trials;
    genes=individuos{1};
    
    [w,s1,s2,alpha,H1,H2,eps_K,phi_r,phi_l,phi0] = print_genes(genes,Nsize);

    xx=[];
    xx1=[];
           
    dist=[];
    RT=[];
        
    %Initialize oscillators

    K=zeros(Nsize,Nsize,1);
    dK=zeros(Nsize,Nsize,1);
    dK=(eps_K~=0).*rand(Nsize,Nsize)*4*pi;
    dK1=zeros(Nsize,Nsize,1);
            
    if plasticity==0
       refK=randi(100,1,1);
       load(['record/record_trial_',num2str(refK)],'K') % Values previously generated for a plastic robot 
       Kf=K(:,:,N/2);
    end

    z=rand(Nsize,N*Ntrials+1)*2*pi;
    phi13=zeros(1,1);
    phi23=zeros(1,1);
    z1=zeros(Nsize,1);
    phi=zeros(Nsize,N*Ntrials); 
    z_h=zeros(Nsize,1); 
            
    for trial=1:1

        if trial~=1
            z(:,1) = z(:,end);
            dK(:,:,1) = dK(:,:,end);
        end

        %Initialize Robot
        x=0;
        y=0;
        p=rand(1,1)*2*pi;

        Mr=0;
        Ml=0;
        vt=0;
        % Sensors
        Ir1=0;
        Il1=0;
        Ir2=0;
        Il2=0;
        % Distance
        D=0;
        % Light 1
        r1=rand(1,1)*2*pi;
        D1=100 + 50*rand(1,1);
        L1=D1*[cos(r1),sin(r1)];
        % Light 2
        r2= r1 + pi/2 + rand(1,1)*pi;
        D2=100 + 50*rand(1,1);
        L2=D2*[cos(r2),sin(r2)];
              
        i0=1;

        for i=1:(N*Ntrials)
                    
            i0=i0+1;


            % Sense the environment
            if  noise ==1 
                Ir1=rand(1,1)*sigma;
                Ir2=rand(1,1)*sigma;
                Il1=rand(1,1)*sigma;
                Il2=rand(1,1)*sigma;
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
                    
            % Compute next state of the oscillator network
            zz=repmat(z(:,i),1,Nsize);
            zz=zz-zz';
            zz=(mod(zz+pi,2*pi)-pi); 

            if plasticity == 1
                K=(1-cos(dK))/2*alpha.*(mod(dK,4*pi)<2*pi); % Weight mapping function
            else
                K=Kf;
            end
                    
            z1 = w + I + diag(K*sin(zz));

            phi(:,i) = angle(K*exp(1j*z(:,i))) - z(:,i);
            phi(:,i) = (mod(phi(:,i)+pi,2*pi)-pi); 
                    

            p_phi = eps_K.*repmat(homeostatic_regulator((phi(:,i)-phi0),H1,H2),1,Nsize);
            diff_z= zz' - repmat(phi0,1,Nsize);
                    

            dK1 = (p_phi).*(1-cos(diff_z))/2;

            z(:,i+1) = z(:,i) + z1*T;
            dK = dK + dK1*T;

            % Compute motor states and move robot            
            Mr = (sin(phi(min([Nsize,4]),i)-phi_r))*Mov;
            Ml = (sin(phi(3,i)-phi_l))*Mov;

            [x,y,p]=move(x,y,p,Ml,Mr,R,T);

            vt=(Mr+Ml)/2;
                    
            d_min=min([sqrt((x-L1(1)).^2 + (y-L1(2)).^2);sqrt((x-L2(1)).^2 + (y-L2(2)).^2)],[],1);
        
            % Re-generate ligths if reached by robot        
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
                
                
            xx=[xx,cos(z)];
      
           if save_on ==1
                savefile=['record/record_trial_',num2str(trial),'.mat'];
            save(savefile);
        end
    end
    
    figure
    plot(dist')
    
    figure
    plot(xx')
    save([name,'_',num2str(rep)],'xx','RT','dist','T')

end

