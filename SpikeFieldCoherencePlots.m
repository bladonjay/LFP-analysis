


load('C:\Users\jhbladon\Desktop\Rapid Analysis\MEC LEC\Designmatrices1-28-15.mat')
load('C:\Users\jhbladon\Desktop\Rapid Analysis\MEC LEC\Cellnames1-26-15.mat')
clear MSmatrix MDmatrix LSmatrix LDmatrix;
%load('C:\Users\jhbladon\Desktop\Rapid Analysis\MEC 2-24-15.mat')
load('C:\Users\jhbladon\Desktop\Rapid Analysis\SIall5-12-15.mat')


Mmatrix2=Mmatrix;
Lmatrix2=Lmatrix;
[Lmatrix,Legend] = easySI2Design(SIalllec);

%[Mmatrix,Legend] = easySI2Design(SIallmec);


stats=LECtheta;
%%

thetaadd=[stats.sess];
%Mmatrix=[Mmatrix thetaadd'];
newmatrix=[Lmatrix thetaadd'];

thetapref=newmatrix(:,end-1);
%thetapref=thetapref(~isnan(thetapref));
thetalock=newmatrix(:,end);
%thetalock=thetalock(~isnan(thetalock));
pluss=thetapref>0; negs=thetapref<0;
% plot just the theta phases in histo
plot(thetapref(pluss),thetalock(pluss),'r*')
hold on; plot(abs(thetapref(negs)),thetalock(negs),'b*');


figure; hist(thetapref,20);
% plot in a polar plot
figure; [y,x]=hist(thetapref,20);
polar(x,y);

figure; plot(thetapref,thetalock,'k*'); hold on
plot(-3:.1:3,cos(-3:.1:3)*.3+.5);
% from the histogram it looks like three separate cell groups, so lets just
% pretend theyre three

% PICK BIN sizes for degree increments and show the regression of the stats
%%
%peakcells=newmatrix(:,end)>.1 & newmatrix(:,end-1)<2.1 & newmatrix(:,end-1)>-1.7 ;
%valleycells=newmatrix(:,end)>.1 & ~peakcells;

figure; plot(thetapref,thetalock,'k*'); hold on
plot(-3:.1:3,cos(-3:.1:3)*.3+.5,'k');

lockthresh=.1;
peakcells=newmatrix(:,end-1)>-1 & newmatrix(:,end-1)<2 & newmatrix(:,end)>lockthresh;
valleycells=abs(newmatrix(:,end-1))>1.5 & newmatrix(:,end)>lockthresh;
middlecells=newmatrix(:,end-1)>.5 & newmatrix(:,end-1)<2.1 & newmatrix(:,end)>lockthresh;
%
plot(newmatrix(peakcells,end-1),newmatrix(peakcells,end),'r*');
hold on;
plot(newmatrix(valleycells,end-1),newmatrix(valleycells,end),'b*');

plot(newmatrix(middlecells,end-1),newmatrix(middlecells,end),'g*');
%%
pcrit=.01;
contextcells=newmatrix(:,3)>.2 & newmatrix(:,10)<pcrit;
itemcells=newmatrix(:,2)>.2 & newmatrix(:,9)<pcrit;
plot(newmatrix(contextcells,end-1),newmatrix(contextcells,end),'r*');
hold on
plot(newmatrix(itemcells,end-1),newmatrix(itemcells,end),'b*');

%%  The average SI and the percent sig for theta locked cells

% this has lows mids and highs

figure;
inds1=[1,4,2,5,3]; %[8,12,9,11,10]; %
MECsis=nanmean(newmatrix(:,inds1));
Mids=nanmean(newmatrix(middlecells,inds1));
Valleys=nanmean(newmatrix(valleycells,inds1));
peaks=nanmean(newmatrix(peakcells,inds1));

% not error of mean though
MECsis2=SEM(newmatrix(:,inds1));
Mids2=SEM(newmatrix(middlecells,inds1));
Valleys2=SEM(newmatrix(valleycells,inds1));
Peaks2=SEM(newmatrix(peakcells,inds1));

bars=[MECsis' Mids' Valleys' peaks'];
errors=[MECsis2' Mids2' Valleys2' Peaks2'];
hbs=bar(bars);
set(hbs(1),'FaceColor','k'); set(hbs(2),'FaceColor','g');
set(hbs(3),'FaceColor','b'); set(hbs(4),'FaceColor','r');
hold on;
exes=[[1:length(inds1)]' [1:length(inds1)]' [1:length(inds1)]' [1:length(inds1)]'];
exes(:,1)=exes(:,1)-.25; exes(:,2)=exes(:,2)-.1;
exes(:,3)=exes(:,3)+.1; exes(:,4)=exes(:,4)+.25;
errorbar(exes,bars,errors,'k.');
set(gca,'XTickLabel',{'Context','Item*Cxt','Item','Item*Pos','Pos','Valence'})
legend('MEC','mids','valley','Peak Locked');
ylabel('Average SI'); xlabel('SI Somparison');
title('SI Trends by Cell Type');
%ylim([0 1]);


%%
figure;
pcrit=.01;
inds2=[8,12,9,11,10];

MECsis=nanmean(newmatrix(:,inds2)<pcrit);
Mids=nanmean(newmatrix(middlecells,inds2)<pcrit);
Valleys=nanmean(newmatrix(valleycells,inds2)<pcrit);
peaks=nanmean(newmatrix(peakcells,inds2)<pcrit);



bars=[MECsis' Mids' Valleys' peaks'];
hbs=bar(bars);
set(hbs(1),'FaceColor','k'); set(hbs(2),'FaceColor','g');
set(hbs(3),'FaceColor','b'); set(hbs(4),'FaceColor','r');
set(gca,'XTickLabel',{'Context','Item*Cxt','Item','Item*Pos','Pos'})
legend('MEC','mid','valley','Peak Locked');
ylabel('Proportion of Cells'); xlabel('SI Somparison');
title('SI Trends by Cell Type');
ylim([0 1]);

%% the nueman way is to do it just in half









%% just for grins lets see if they split by average firing rate:
figure;
inds1=  [1,4,2,5,3]; % [8,12,9,11,10]; % 
ratethresh=10;

ratethresh1=7; ratethresh2=12;

tops=newmatrix(:,7)>ratethresh2;
bottoms=newmatrix(:,7)<ratethresh1;

Highs=nanmean(newmatrix(tops,inds1));
Lows=nanmean(newmatrix(bottoms,inds1));
Highs2=SEM(newmatrix(tops,inds1));
Lows2=SEM(newmatrix(bottoms,inds1));
exes=[[1:length(inds1)]' [1:length(inds1)]'];
exes(:,1)=exes(:,1)-.15; exes(:,2)=exes(:,2)+.15;

bars=[Highs' Lows'];
errors=[Highs2' Lows2'];
hbs=bar(bars); hold on
errorbar(exes,bars,errors,'k.');
set(gca,'XTickLabel',{'Context','Item*Cxt','Item','Item*Pos','Pos'})
legend('Highs','Lows');
titlename=['rate cutoff ' num2str(ratethresh) ' Hz '];
title(titlename);
%%
figure;
inds1=[8,12,9,11,10]; %[1,4,2,5,3]; %
ratethresh=12; pcrit=.05;

%ratethresh1=7; ratethresh2=12;

tops=newmatrix(:,7)>ratethresh;
bottoms=newmatrix(:,7)<ratethresh;

Highs=nanmean(newmatrix(tops,inds1)<pcrit);
Lows=nanmean(newmatrix(bottoms,inds1)<pcrit);

bars=[Highs' Lows'];

hbs=bar(bars); hold on
%errorbar(exes,bars,errors,'k.');
set(gca,'XTickLabel',{'Context','Item*Cxt','Item','Item*Pos','Pos'})
legend('Highs','Lows');
titlename=['rate cutoff ' num2str(ratethresh) ' Hz '];
title(titlename);


