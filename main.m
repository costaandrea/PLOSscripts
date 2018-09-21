%%%Scripts to replicate the figures in "On the Calculation of Betweenness 
%%%Centrality in MarineConnectivity Studies Using Transfer Probabilities"
%%%Costa et al., (2017); PlosOne
%%%
%%%This script requires the m_map package
%%%(https://www.eoas.ubc.ca/~rich/map.html)
%%%
%%%Author: Andrea Costa 2017/01/19 - andrea.costa@univ-amu.fr
%%%
%%%Credit for the connectivity matrices to Katell Guizien: "Vulnerability of 
%%%marine benthic metapopulations: implications of spatially structured 
%%%connectivity for conservation practice." Divers. Distrib. 2014;20(12):1392-1402
%%%

clear all
close all

%add useful stuff
addpath('F:\ANDREA\NUMERICAL_METHODS\m_map')


%load data
connectivity = importdata('connectivity_10_30_composite_2004_2006notext.dat');

B=reshape(connectivity',32,32,20);
%change names. I'll use it later
Ci=B;
adjme=B;


%choose what to do
metric = 'old'; %old or new %%you should use old-low and new-high to reproduce the paper's results
range = 'low'; %low or high


Q=[7]; %which matrices you wanna analyze
%Q=1:20;

for iii=1:numel(Q)
kkk=Q(iii);

if strcmp(metric,'new')
 Cu=ones(32,32);
 xx=Ci(:,:,kkk);
 
vals=xx(:);
[v, in] = sort(vals);%sort ascendent

indx=in(v>0);

hm = floor(numel(indx)/20);%index of 5 per cent 

if strcmp(range,'high')
 erase = indx(1:end-hm);%erase lower 95 per cent
elseif strcmp(range,'low')
 erase = indx(hm+1:end);%erase top 95 per cent
end

xx(erase)=0;

adj = Cu ./ xx;  %1/x

 adj1=zeros(32,32);
  for i=1:32
   for j=1:32
     adj1(i,j) = log(adj(i,j)); %log(1/x)
      
   end
  end


elseif strcmp(metric,'old')
    
  xx=Ci(:,:,kkk);
 
vals=xx(:);
[v, in] = sort(vals);%sort ascendent

indx=in(v>0);

hm = floor(numel(indx)/20);%index of 5 per cent 

if strcmp(range,'high')
 erase = indx(1:end-hm);%erase lower 95 per cent
elseif strcmp(range,'low')
 erase = indx(hm+1:end);%erase top 95 per cent
end

xx(erase)=0;

adj1=xx;   
    
end

adjme=adj1;
if strcmp(metric,'new')
    adjme(find(isinf(adjme)))=0; %do not plot Inf values  
end


%define coordinates on a circle
theta=linspace(0,2*pi,33);
[LONG,LAT]=pol2cart(theta,1);

    
N=length(adjme); %number of sites

adjm=adjme-eye(N).*adjme; %remove diagonal

weight=adjm;  %define edges' weight


minno=min(min(adjm(adjm>0)));
maxxo=max(max(adjm));


figure('Name',['matrix_',num2str(Q(iii)),'_',num2str(minno),'_',num2str(maxxo)])

m_proj('miller','lon',[-1.3 1.3],'lat',[-1.3 1.3]);
set(gca,'ytick',[],'xtick',[])
%m_gshhs_h('patch',[.5 .5 .5]);
m_grid('ytick',[],'yticklabels',[],'xtick',[],...
       'xticklabels',[],'linestyle','none','ticklen',.02)
%m_usercoast('mm','patch',[.5 .5 .5]);
hold on


%%%Choose the first colormap
cm = winter;
colormap(cm);

aa=0;
bb=0;
cc=0;


for i=N:-1:1   %%m_quiver

  C = find(adjm(i,:)); %find connections spanning the lines of adjm

    for ii=1:numel(C) %calculates components of vectors starting from node i
    
      k=C(ii); %more manageble index
         
         if (adjm(i,k)>=minno && adjm(i,k)<=maxxo) || (adjm(k,i)>=minno && adjm(k,i)<=maxxo)            
   
        LONGI=[LONG(i) LONG(k)]; %geographical coordinates of the two Areas
        LATI=[LAT(i) LAT(k)];

           indCol = findnearest((1:length(cm))*(maxxo-minno)/length(cm),(adjm(i,k)-minno));
           col = cm(indCol,:);
           
       if (adjm(i,k)>=minno && adjm(i,k)<=maxxo) && ~(adjm(k,i)>=minno && adjm(k,i)<=maxxo) %plot only the first one
           aa=aa+1;
          [su,a12,a21] = m_idist(LONGI(1), LATI(1), LONGI(2), LATI(1));
          [lonu,latu,a21] = m_fdist(LONGI(1), LATI(1), a12, 4*su/5);
           
          longs=[LONGI(1),lonu];
          lats=[LATI(1),latu];
          
          U=m_lldist(longs,lats);
          
          if LONGI(1)>LONGI(2)  %change sign of U-component if necessary
              U= -U;
          end
              
          [sv,a12,a21] = m_idist(LONGI(1), LATI(1), LONGI(1), LATI(2));
          [lonv,latv,a21] = m_fdist(LONGI(1), LATI(1), a12, 4*sv/5);
         
          longs=[LONGI(1),lonv];
          lats=[LATI(1),latv];
          
          V=m_lldist(longs,lats);
          
          if LATI(1)>LATI(2)  %change sign of V-component if necessary
              V= -V;
          end
          
            m_quiver(LONGI(1), LATI(1), U/100, V/100, 'Color', col,'LineWidth', 2);

            hold on
            m_line([LONGI(1) LONGI(2)], [LATI(1) LATI(2)], 'Color',col,'LineWidth', 2);
            hold on
            
       elseif (adjm(i,k)>=minno && adjm(i,k)<=maxxo) && (adjm(k,i)>=minno && adjm(k,i)<=maxxo)   %plot both
           bb=bb+1;
          [su,a12,a21] = m_idist(LONGI(1), LATI(2), LONGI(2), LATI(2));
          [lonu,latu,a21] = m_fdist(LONGI(1), LATI(1), a12, su/2);
              
          longs=[LONGI(1),lonu];
          lats=[LATI(1),latu];
          
          U=m_lldist(longs,lats);
          
          if LONGI(1)>LONGI(2)  %change sign of U-component if necessary
              U= -U;
          end
              
          [sv,a12,a21] = m_idist(LONGI(1), LATI(1), LONGI(1), LATI(2));
          [lonv,latv,a21] = m_fdist(LONGI(1), LATI(1), a12, sv/2);
           
          longs=[LONGI(1),lonv];
          lats=[LATI(1),latv];
          
          V=m_lldist(longs,lats);
          
          if LATI(1)>LATI(2)  %change sign of V-component if necessary
              V= -V;
          end
          
          indCol2 = findnearest((1:length(cm))*(maxxo-minno)/length(cm),(adjm(k,i)-minno));
           col2 = cm(indCol2,:);
           
            m_quiver(LONGI(1), LATI(1), U/100, V/100, 'Color', col, 'LineWidth', 2);
              hold on

            m_quiver(LONGI(2), LATI(2), -U/100, -V/100, 'Color', col2, 'LineWidth', 2);
            
             hold on
          
            
      elseif ~(adjm(i,k)>=minno && adjm(i,k)<=maxxo) && (adjm(k,i)>=minno && adjm(k,i)<=maxxo) %plot only the second one
          cc=cc+1;
          [su,a12,a21] = m_idist(LONGI(1), LATI(1), LONGI(2), LATI(1));
          [lonu,latu,a21] = m_fdist(LONGI(1), LATI(1), a12, 4*su/5);
           
          longs=[LONGI(1),lonu];
          lats=[LATI(1),latu];
          
          U=m_lldist(longs,lats);
          
          if LONGI(1)>LONGI(2)  %change sign of U-component if necessary
              U= -U;
          end
              
          [sv,a12,a21] = m_idist(LONGI(1), LATI(1), LONGI(1), LATI(2));
          [lonv,latv,a21] = m_fdist(LONGI(1), LATI(1), a12, 4*sv/5);
         
          longs=[LONGI(1),lonv];
          lats=[LATI(1),latv];
          
          V=m_lldist(longs,lats);
          
          
          if LATI(1)>LATI(2)  %change sign of V-component if necessary
              V= -V;
          end
          
          indCol2 = findnearest((1:length(cm))*(maxxo-minno)/length(cm),(adjm(k,i)-minno));
           col2 = cm(indCol2,:);

            m_quiver(LONGI(2), LATI(2), -U/100, -V/100, 'Color', col2, 'LineWidth', 2);
              hold on

            m_line([LONGI(2) LONGI(1)], [LATI(2) LATI(1)], 'Color',col2,'LineWidth', 2);
            hold on
      end
      
        end %$%
        
    end
          
end

hold on

hcb=colorbar('EastOutside');
set(hcb,'position',[.85 .11 .05 .815])
caxis([minno, maxxo]);

if strcmp(metric,'new')
    
set(hcb,'YTick',linspace(minno,maxxo,6));
 lbl=linspace(minno,maxxo,6)';
lbl2 = num2str(lbl,'%0.1f');

 set(hcb,'yticklabel',lbl2,'FontSize',16)
 
elseif strcmp(metric,'old')
    
 set(hcb,'YTick',linspace(minno,maxxo,6));

 lbl=linspace(minno,maxxo,6)';

 if strcmp(range,'low')
  set(hcb,'yticklabel',num2str(lbl*1e5,'%0.1f'),'FontSize',16)
 
  annotation('textbox', [0.93, 0.99, 0, 0], 'string', 'x10^{-5}','FontSize',14)
  
 elseif strcmp(range,'high')
  set(hcb,'yticklabel',num2str(lbl*1e2,'%0.1f'),'FontSize',16)
 
  annotation('textbox', [0.93, 0.99, 0, 0], 'string', 'x10^{-2}','FontSize',14)
 end
 
end
 
%freeze colorbar
cbfreeze(hcb);



%%%NODE BETWEENNESS
betw1=zeros(N,1);
for kkk=1:numel(Q)
 Cu=ones(32,32);
 adj = Cu ./ Ci(:,:,Q(kkk));
 adj1=zeros(32,32);
  for i=1:32
   for j=1:32
     adj1(i,j) = log(adj(i,j)); %log(1/x)
      
   end
  end
  
if strcmp(metric,'old')  
  adj1=Ci(:,:,kkk);
end

adjme=adj1;  
adjme(find(adjme==Inf))=0;


betw1 = betweenness_wei(adjme);

weight = betw1;


%%% New colorbar
cm = autumn;
colormap(cm);

 for j=1:N
        
        indCol = findnearest((1:length(cm))*max(weight)/length(cm),weight(j));
        col = cm(indCol,:);  
        

     m_plot(LONG(j),LAT(j),'.','MarkerSize',42,'Color', col)  %draw markers
     m_plot(LONG(j),LAT(j),'o','MarkerSize',14,'Color','k')
 end

hold on

 hC = colorbar('WestOutside');
 set(hC,'position',[.12 .11 .05 .815])


set(hC,'XTick',[0,0.4]);
set(hC,'XTick',[0,1]);
 caxis([0, max(weight)]);

theta=linspace(0,2*pi,33);
[LONG2,LAT2]=pol2cart(theta,1.2);
LONG2=LONG2+0.05;
LAT2=LAT2-0.07;
 for i=1:N
  [X,Y]=m_ll2xy(LONG2(i),LAT2(i));

  Mstr=num2str(i);

  text(X,Y,num2str(i),'Color','k','FontWeight','bold',...
  'FontSize',14, 'HorizontalAlignment','right', 'VerticalAlignment','bottom')

 end
end
    set(gca,'fontsize',18)
end

          

%%%All betweennessess
for kkk=1:20
    betwenn(:,kkk) = betweenness_wei(newproba(B(:,:,kkk)));
end

figure
colormap(hot)
imagesc(betwenn)
h=colorbar;
ylabel(h,'Betweenness','fontsize',14)
set(gca,'fontsize',14)
ylabel('Sites','fontsize',14)
xlabel('Matrices','fontsize',14)
