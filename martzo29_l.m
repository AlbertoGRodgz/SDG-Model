%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Barrio's Code %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Definicion de variables%%%
clear all;

%%   parameters
N =6; % number of countries
dt=.0005;
gi=0.6;
gk=20;
neighbors=0;
goal_l = 50;
mem = 5;
ngoals = 17;     %number of goals
iterations = 500;
lower_limit=30;   % wealth distribution limits
upper_limit=100;

r_anual = .1; % porcentaje de recururso por a?o
xn=zeros(ngoals,N,iterations);

%% initial condityions
g = randi(N,ngoals,N);
x = g;           %  interaction strngth vetween countris and goals

%A =  importdata('countriesM.txt') ;%  Countries vonnections matrix 
%A=spones(N);
Z=6;   % Z mero de vecinos
A=sprandsym(N,Z/N);
A=A/(max(max(A))-min(min(A)));
A=A-diag(diag(A));

A=[0,1,0,0,0,0;1,0,0,0,0,0;0,0,0,1,0,0;0,0,1,0,0,0;0,0,0,0,0,1;0,0,0,0,1,0];
n = find(A ~= 0);
[I,J] = ind2sub(size(A),n);  %      first neighbors A;
A2 = A.^2;
A2=A2-diag(diag(A2));
nn = find(A2 ~= 0);
[II,JJ] = ind2sub(size(A2),nn);%   second neighbors A2;

%init=randi([0 12],ngoals,N);
ggg=load('scores_2000.csv');
ggg=ggg';
init(1:17,1:6)=ggg(1:17,1:6)/2;
for ii=1:mem
    xn(:,:,ii)=sort(init,2,'descend');    % initial values for memory
end
xn=sort(xn);
glinks =  importdata('corelacionM.txt') ; %hoals correlation matriz

%%   Generate a vector of values sampled from a Pareto distribution within a specified range.

% U = rand(1, max(N));
% X = (1 ./ U).^(1.^(1. / 2)) ;
% pareto_samples_normalized = (X-min(X))/(max(X)-min(X)); %Normalize the Pareto samples to [0, 1]
% 
% pareto_vec = lower_limit + pareto_samples_normalized * (upper_limit - lower_limit);
% r = pareto_vec;
r= 500*[1 .95 .2 02 .007 .001];
r=sort(r,'ascend');
gg=ones(ngoals,N)/2;
qq=sum(gg);
for i=1:ngoals
    for j=1:N
        e(i,j)=gg(i,j)/qq(j);
    end
end

%r_anual=sort(r,'ascend')*.1;
%%   figures
figure(11)
clf
surface(glinks)
figure(12)
clf
surface(A)

figure(13)
clf
surface(A2)

figure(14)
clf
rpp=[
0.067004382
0.009353274
0.091714877
0.430166578
0.034276889
0.407325773];
loglog(sort(r,'descend'))
hold on
loglog(sort(rpp*1000,'descend'))
r=sort(rpp*1000,'descend');

%%
%load test1
r=r*2.6;

%% Dynamical loop
for cont=mem+1:iterations
    
    for i=1:N
        sources(:,cont)=r;
        
        for j=1:ngoals
            
            
            for k=length(n)
                neigh(I(k),J(k))=g(j,i);
            end
             neighbors =sum(sum(neigh));
            for k=length(nn)
                fneigh(II(k),JJ(k))=g(j,i);
            end
             fneighbors =sum(sum(fneigh));
            
            ginteractions = 0;
            fginteractions = 0;
            for k=1:ngoals
                ginteractions = (glinks(k,j) * x(j,i)) *gi + ginteractions;
                %fginteractions =  x(j,i) + fginteractions;
            end
            
            memory = 0;
            for m=cont-mem:cont-1
                memory = memory + xn(j,i,m);
            end
            xn(j,i,cont) =xn(j,i,cont-1)+ (r(i)*(r_anual)*((1-e(j,i))*(goal_l-x(j,i))/goal_l)+...
                gk/2*(fneighbors/length(nn))+gk*(neighbors/length(n))+ginteractions+abs((memory/mem)))*dt ;
            
            if xn(j,i,cont)>=goal_l
                xn(j,i,cont)=goal_l;
            end
            
            %neighbors=0;
        end
    end
    x = xn(:,:,cont);
    %%
    sxn=sum(xn)/ngoals;
    %  for i=1:ngoals
    xg(:,:)=goal_l-sxn(1,:,:);
    %end
    figure(1)
    clf

    % Graficar cada grupo con su respectivo color
    
    plot(xg(1, :)', 'color', [0 0 .7], 'linewidth',2)  % Rojo
    hold on
    box on
    plot(xg(2, :)', 'color',[.5 0 .5],  'linewidth', 2) % Azul
    
    plot(xg(3, :)', 'color',[1 0.1 0.1],  'linewidth', 2) % Verde
    plot(xg(4, :)', 'color',[ 1 .8 .6], 'linewidth', 2) % Azul
    plot(xg(5, :)', 'color',[0 .8 .0], 'linewidth', 2) % Verde
    plot(xg(6, :)', 'color',[.6 .5 .05], 'linewidth',2) % Verde
    ylabel('Distance to goal achievement')
    xlabel('Time in iterations (arb units)')
    set(gca, 'fontsize', 16, 'fontweight', 'bold', 'linewidth', 2)
    %axis([mem iterations 0 goal_l])
    title('Average behaviour of countries')
    
    % Definir etiquetas en la leyenda
    cl = char('USA-EU', 'UKC empire','Rusia','China-India','LL America','Agrica');
    legend(cl,'location', 'southoutside', 'orientation', 'vertical')
        set(gca, 'fontsize', 16, 'fontweight', 'bold', 'linewidth', 2)

    pause(0.01)
    hold off
end

%%
%      figure(2)
% clf
%
%
%     xg1(:,:)=goal_l-xn(1,:,:);
%     xg2(:,:)=goal_l-xn(2,:,:);
%     xg3(:,:)=goal_l-xn(3,:,:);
%
% for ii=1:ngoals
%     for jj=1:N
%         xg1(ii,:)=goal_l-xn(ii,jj,:);
%         hold on
%      subplot(ngoals,1,ii),plot(xg1(ii,:)','linewidth',2)
%     set(gca,'fontsize',16,'fontweight','bold','linewidth',2)
%     ylabel(['oal ', num2str(ii)])
%     pause(0.01)
%     end
% end
%     subplot(3,1,2), plot(xg2','linewidth',2)
%     set(gca,'fontsize',16,'fontweight','bold','linewidth',2)
%     %xlabel('tiem (iterations x 50') ylabel('goal 2')%axis([0 t -1.5 1.5])
%
%     subplot(3,1,3),plot(xg3','linewidth',2)
%     set(gca,'fontsize',16,'fontweight','bold','linewidth',2)
%      xlabel('time in iterations (arb units)')
%    ylabel('goal 3')%axis([0 t -1.5 1.5])
%%

figure(3)
clf
%     xg1(:,:)=goal_l-xn(1,1:N,:);
%     pl1=plot(xg1','linewidth',4,'color','r');
%     %legend('goal 1')
%
% hold on
%cc=char(['r','y','b','m','k','c']);
% for i=1:ngoals
%         xg1(:,:)=goal_l-xn(i,1:N,:);
%             %pip=['pl',num2str(i)];
%             plot(xg1','linewidth',4,'color',cc(i));
%          pause(0.01)
% end
tt=goal_l-sum(xn,2)/N;
for i=1:ngoals
    hold on
    xg1(i,:)=tt(i,:);%xn(i,N,:);
    pip=[];
    plot(xg1(i,:)','linewidth',4);
    pip=([pip 'country',num2str(i)]);
    
    pause(0.01)
end
%
%sl=['goal 1','goal 2','goal 3','goal 4','goal 5'];
set(gca,'fontsize',16,'fontweight','bold','linewidth',2)
ylabel('Distance to goal achivement')%
xlabel('time in iterations (arb units)')
title('average  behaviour of goals')
%       axis([mem iterations 0 goal_l])

%     sl=['goal ' ,'1',]
%  for i=2:9
%      sl=[sl;['goal ' num2str(i) ,]];
%  end
sg=1:ngoals;
sl=num2str(sg');
legend(sl,'location','eastoutside')
box on
hold off


%%
%
% figure(5)
% clf
% sg1=sum((xg(1:4, :)))/4;
% sg2=sum((xg(5:8, :)))/4; %
% sg3=sum((xg(9:11, :)))/3; %
% sg4=sum((xg(12:14, :)))/3; %
% sg5=sum((xg(15:20, :)))/6; %
% 
% 
% plot(sg1', 'r', 'linewidth',2)  % Rojo
% set(gca,'fontsize',16,'fontweight','bold','linewidth',2)
% 
% hold on
% plot(sg2', 'b',  'linewidth', 2) % Azul
% 
% plot(sg3, 'g',  'linewidth', 2) % Verde
% plot(sg4, 'k', 'linewidth', 2) % Azul
% plot(sg5, 'y', 'linewidth', 2) % Verde
% legend('gluster1','gcuster2','cluster3','cluuster4','cluster5')
% xlabel('time in iterations (arb units)')
% box on
% title('average  behaviour of clusters')
% ylabel('Distance to goal achivement')%
% 
% 
% hold off


