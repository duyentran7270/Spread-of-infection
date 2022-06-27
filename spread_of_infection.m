%% Case 1: We suppose that the exposed students are not infectious immediately 
clc
clear all
close all 

nsim = 1e3; % number of simulations

sq_side = 10; % (meters) side length of a square 
npers = 51; % total number of students (one of them is infected initially) 
std = 0.3; % parameter of the normal distribution
min_dist = 0.5; % (meters) if the infected student comes closer than 0.5 m to others, they risk to be infected
prob = 0.8; % probability to get infected
Tend = 4 * 60; % time when the simulation ends (4 hours = 240 minutes)
tstep = 5; % (minutes) time step

t=0:tstep:Tend;
nstep=length(t)-1;

for k=1:nsim
    i=1;    
    
    % x=10*rand(50,1);x(51,1)=0;
    % y=10*rand(50,1);y(51,1)=0;
    x=10*rand(51,1);
    y=10*rand(51,1);
    
    infect=[51];
    
    for i=2:nstep 
        x(:,i)=x(:,i-1)+std*randn(51,1);
        y(:,i)=y(:,i-1)+std*randn(51,1);
        x(find(x>10))=10;x(find(x<0))=0;
        y(find(y>10))=10;y(find(y<0))=0;
        add_in=[];
        tak=[];
        for j=1:50;
            dis=sqrt((x(end,i)-x(j,i)).^2+(y(end,i)-y(j,i))^2);
            if dis<=0.5 
                add_in=[add_in j];
            end
        end
        ni=length(add_in);
        tak=add_in(randperm(ni,round(0.8*ni)));  
        infect=[infect tak];
        infect=unique(infect);
    end
    
    ninfe(k)=length(infect);
    
    if k==1, 
        infe=infect;
        no_infe=1:51;
        no_infe(infe)=[];
        xp=x(:,end);yp=y(:,end);
    end
end 

mean_infected=mean(ninfe)

figure
plot(xp(infe,end),yp(infe,end),'r*');
hold on
plot(xp(no_infe,end),yp(no_infe,end),'ko');
grid on 
title('Final locations of students');
legend('infected students','healthy students','location','bestoutside');


%% Case 2: We suppose that the infected students immediately infect others

clc
clear all
close all 

nsim = 1e3; % number of simulations

sq_side = 10; % (meters) side length of a square 
npers = 51; % total number of students (one of them is infected initially) 
std = 0.3; % parameter of the normal distribution
min_dist = 0.5; % (meters) if the infected student comes closer than 0.5 m to others, they risk to be infected
prob = 0.8; % probability to get infected
Tend = 4 * 60; % time when the simulation ends (4 hours = 240 minutes)
tstep = 5; % (minutes) time step

t=0:tstep:Tend; % time
nstep=length(t)-1;

for k=1:nsim
    i=1;    
    
    % x=10*rand(50,1);x(51,1)=0;
    % y=10*rand(50,1);y(51,1)=0;
    x=10*rand(51,1);    %generate random walk
    y=10*rand(51,1);
    
    infect=[51];
    
    % figure (i)
    % scatter(z(1,i),z(2,i),'ko','filled'); hold on
    % plot(x(1:50,i),y(1:50,i),'ro'); hold off
    
    for i=2:nstep 
        x(:,i)=x(:,i-1)+std*randn(51,1);
        y(:,i)=y(:,i-1)+std*randn(51,1);
        x(find(x>10))=10;x(find(x<0))=0;
        y(find(y>10))=10;y(find(y<0))=0;
        add_in=[];
        tak=[];
        for j=1:50;
            dis=sqrt((x(end,i)-x(j,i))^2+(y(end,i)-y(j,i))^2);
            if dis<=0.5 
                add_in=[add_in j];
            end
            if length(infect)>1
            for m=1:length(infect)-1
                dism=sqrt((x(infect(m),i)-x(j,i))^2+(y(infect(m),i)-y(j,i))^2);
                if dism<=0.5
                    add_in=[add_in j];
                end 
            end
            end
        end
        add_in=unique(add_in);
        ni=length(add_in);
        tak=add_in(randperm(ni,round(0.8*ni)));  
        infect=[infect tak];
        infect=unique(infect);
    end
    
    ninfe(k)=length(infect);
    
    if k==1, 
        infe=infect;
        no_infe=1:51;
        no_infe(infe)=[];
        xp=x(:,end);yp=y(:,end);
    end 
end

mean_infected=mean(ninfe)

figure
plot(xp(infe,end),yp(infe,end),'r*');
hold on
plot(xp(no_infe,end),yp(no_infe,end),'ko');
grid on 
title('Final locations of students');
legend('infected students','healthy students','location','bestoutside');