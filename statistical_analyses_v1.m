%cd('/Users/scooley2/Dropbox (University of Oregon)/Reservoir_Review/datasets');
cd('/Users/sarahcooley/Dropbox (University of Oregon)/Reservoir_Review/datasets');
load('complete_dataset_mar22.mat');

%clear out_data
c = compare_data;
num_obs = [c.num_obs]';
c(num_obs < 2) = [];
num_obs = [c.num_obs]';

for i = 1:length(c)
    type(i,1) = get_type(c(i).valid);
    sv(i,:) = c(i).sea_var_abs;
    svv(i,:) = c(i).sea_var;
    tv(i,:) = c(i).tot_var_abs;
    tvv(i,:) = c(i).tot_var;
    rv(i,:) = c(i).rmse;
    av(i,:) = c(i).ann_err_abs;
    avv(i,:) = c(i).ann_err;
    vol(i,1) = c(i).rv_mcm;
    rr(i,:) = c(i).rout;
    
end

%what is the average error (by group)
%type(type < 4) = 4;
t = 1;
clear o
var = 100*rv(type == t,:)./vol(type == t);
o(1,:) = median(var);
var = 100*svv(type == t,:)./vol(type == t);
o(2,:) = median(var);
var = 100*sv(type == t,:)./vol(type == t);
o(3,:) = median(var);
var = 100*avv(type == t,:)./vol(type== t);
o(4,:) = median(var);
var = 100*av(type == t,:)./vol(type == t);
o(5,:) = median(var);
var = 100*tvv(type == t,:)./vol(type == t);
o(6,:) = median(var);
var = 100*tv(type == t,:)./vol(type == t);
o(7,:) = median(var);
var = rr(type == t,:);
o(8,:) = median(var);

%rv = 100*rv./vol;
%rv = 100*rv./vol;


%% boxplots of variability
rp = 100*rv./vol;
%RMSE
figure(1)
subplot(2,2,1)
var = rp(type == 1,:);
boxchart(var,'MarkerStyle','.')
ylim([0 100])
xticklabels({'GLWS','GRS','GloLakes','GRDL-Y','GRDL-L'})
ylabel('Mean RMSE (% Capacity)')
set(gca,'FontSize',14);
title(['n = ' num2str(sum(type == 1))])
box on

subplot(2,2,2)
var = rp(type <= 2,:);
var(:,3) = [];
boxchart(var,'MarkerStyle','.')
ylim([0 100])
xticklabels({'GLWS','GRS','GRDL-Y','GRDL-L'})
ylabel('Mean RMSE (% Capacity)')
set(gca,'FontSize',14);
title(['n = ' num2str(sum(type <= 2))])
box on


subplot(2,2,3)
type1 = type;
type1(type == 1) = 3;
var = rp(type1 == 3,:);
var(:,[1 4]) = [];
boxchart(var,'MarkerStyle','.')
ylim([0 100])
xticklabels({'GRS','GloLakes','GRDL-L'})
ylabel('Mean RMSE (% Capacity)')
set(gca,'FontSize',14);
title(['n = ' num2str(sum(type1 == 3))])
box on


subplot(2,2,4)
var = rp(type <= 4,:);
var(:,[1 3 4]) = [];
boxchart(var,'MarkerStyle','.')
ylim([0 100])
xticklabels({'GRS','GRDL-L'})
set(gca,'FontSize',14);
ylabel('Mean RMSE (% Capacity)')
title(['n = ' num2str(sum(type <= 4))])
box on

%% figure 2 - boxplot by type
figure(2)
subplot(2,3,1)
hold off
var = 100*svv(type == 1,:)./vol(type == 1,:);
boxchart(var,'MarkerStyle','.')
xticklabels({'GLWS','GRS','GloLakes','GRDL-Y','GRDL-L'})
ylabel('Seasonal Variability Dif (% Capacity)')
ylim([ -100 100])
set(gca,'FontSize',14)
box on

subplot(2,3,2)
hold off

var = 100*avv(type == 1,:)./vol(type == 1,:);
boxchart(var,'MarkerStyle','.')
xticklabels({'GLWS','GRS','GloLakes','GRDL-Y','GRDL-L'})
ylabel('Annual Variability Dif (% Capacity)')
ylim([ -50 50])
set(gca,'FontSize',14)
box on



subplot(2,3,3)
hold off

var = 100*tvv(type == 1,:)./vol(type == 1,:);
boxchart(var,'MarkerStyle','.')
xticklabels({'GLWS','GRS','GloLakes','GRDL-Y','GRDL-L'})
ylabel('Total Variability Dif (% Capacity)')
ylim([ -225 225])
set(gca,'FontSize',14)
box on


subplot(2,3,4)
var = 100*sv(type == 1,:)./vol(type == 1,:);
boxchart(var,'MarkerStyle','.')
xticklabels({'GLWS','GRS','GloLakes','GRDL-Y','GRDL-L'})
ylabel('Seasonal Variability Dif (% Capacity)')
ylim([ 0 100])
set(gca,'FontSize',14)
box on


subplot(2,3,5)
var = 100*av(type == 1,:)./vol(type == 1,:);
boxchart(var,'MarkerStyle','.')
xticklabels({'GLWS','GRS','GloLakes','GRDL-Y','GRDL-L'})
ylabel('Annual Variability Dif (% Capacity)')
ylim([ 0 50])
set(gca,'FontSize',14)
box on


subplot(2,3,6)
var = 100*tv(type == 1,:)./vol(type == 1,:);
boxchart(var,'MarkerStyle','.')
xticklabels({'GLWS','GRS','GloLakes','GRDL-Y','GRDL-L'})
ylabel('Total Variability Dif (% Capacity)')
ylim([ 0 225])
set(gca,'FontSize',14)
box on


%% statistics for each individual value
%var = rv;
%var = nanmean(var,2);
clear o
for p = 1:4
    if p == 1; v = rv; end
    if p == 2; v = sv; end
    if p == 3; v = av; end
    if p == 4; v = tv; end
    var = 100*nanmean(v,2)./vol;

    o(1,p) = mean(var);
    o(2,p) = median(var);
    o(3,p) = prctile(var,25);
    o(4,p) = prctile(var,75);
    o(5,p) = 100*sum(abs(var) < 5)./length(var);
    o(6,p) = 100*sum(abs(var) < 10)./length(var);
    o(7,p) = 100*sum(abs(var) < 25)./length(var);
    o(8,p) = 100*sum(abs(var) < 50)./length(var);
end

%% what predicts error?
area = [c.area]';
lat = [c.lat]';
sdi = [c.sdibc]';
cloud = [c.cloud]';
tri = [c.tri_elv]';
year = [c.year]';
vol = [c.rv_mcm]';
rm = 100*rv./vol;
var = nanmean(rm,2);
num_obs = [c.num_obs]';
gvar = 100*[c.mean_sea_var]'./vol;

clear o
%var = gvar;
[r,p] = corrcoef(var,log(area));
o.arear = r(2);
o.areap = p(2);
[r,p] = corrcoef(var,log(vol));
o.volr = r(2);
o.volp = p(2);
[r,p] = corrcoef(var,cloud);
o.cloudr = r(2);
o.cloudp = p(2);
[r,p] = corrcoef(var(~isnan(tri)),tri(~isnan(tri)));
o.trir = r(2);
o.trip = p(2);
[r,p] = corrcoef(var,lat);
o.latr = r(2);
o.latp = p(2);
[r,p] = corrcoef(var,sdi);
o.sdir = r(2);
o.sdip = p(2);
[r,p] = corrcoef(var,num_obs);
o.nor = r(2);
o.nop = p(2);
[r,p] = corrcoef(var(year > 0),year(year > 0));
o.yearr = r(2);
o.yearp = p(2);
[r,p] = corrcoef(var,gvar);
o.gvarr = r(2);
o.gvarp = p(2);


figure(3)
subplot(2,3,1)
cmap = lines(5);
cc = cmap(1,:);
semilogx(vol,var,'o','MarkerFaceColor',cc,'MarkerEdgeColor',cc);
alpha(0.5);
box on
ylim([0 200])

xlabel('Capacity (MCM)');
ylabel('RMSE (% Capacity)')
set(gca,'FontSize',14)
subplot(2,3,2)
semilogx(area,var,'o','MarkerFaceColor',cc,'MarkerEdgeColor',cc);
box on
xlabel('Area (km^2)');
ylabel('RMSE (% Capacity)')
set(gca,'FontSize',14)
ylim([0 200])

subplot(2,3,3)
scatter(gvar,var,'o','MarkerFaceColor',cc,'MarkerEdgeColor',cc);
alpha(0.5)
box on
xlabel('Mean Seasonal Variability (% Capacity)');
ylabel('RMSE (% Capacity)')
set(gca,'FontSize',14)
ylim([0 200])

subplot(2,3,4)
scatter(sdi,var,'MarkerFaceColor',cc,'MarkerEdgeColor',cc);
alpha(0.5)
box on
xlabel('Shoreline Development Index');
ylabel('RMSE (% Capacity)');
set(gca,'FontSize',14)
axis([0 50 0 200])

subplot(2,3,5)
scatter(tri,var,'MarkerFaceColor',cc,'MarkerEdgeColor',cc);
alpha(0.5)

box on
xlabel('Terrain Ruggedness Index');
ylabel('RMSE (% Capacity)');
set(gca,'FontSize',14)

ylim([0 200])

subplot(2,3,6)
scatter(cloud,var,'MarkerFaceColor',cc,'MarkerEdgeColor',cc);
alpha(0.5)

box on
xlabel('Cloudiness (%)');
ylabel('RMSE (% Capacity)');
set(gca,'FontSize',14)
ylim([0 200])


x = [var gvar  tri sdi cloud log(vol) log(area) num_obs abs(lat) year];
%z = [log(vol) gvar];
II = ~isnan(tri) & year > 0;
x = x(II,:);
%z = z(II,:);
[rho,pval] = partialcorr(x);
rho = array2table(rho, ...
    'VariableNames',{'RMSE','Sea Var','TRI','SDI','Cloud','Log(Vol)','Log(Area)','Num Obs','Abs(Lat)','Year'},...
    'RowNames',{'RMSE','Sea Var','TRI','SDI','Cloud','Log(Vol)','Log(Area)','Num Obs','Abs(Lat)','Year'});
pval = array2table(pval, ...
    'VariableNames',{'RMSE','Sea Var','TRI','SDI','Cloud','Log(Vol)','Log(Area)','Num Obs','Abs(Lat)','Year'},...
    'RowNames',{'RMSE','Sea Var','TRI','SDI','Cloud','Log(Vol)','Log(Area)','Num Obs','Abs(Lat)','Year'});

%% group by capacity and variability
vv = vol;
t1 = prctile(vv,20);
t2 = prctile(vv,40);
t3 = prctile(vv,60);
t4 = prctile(vv,80);
t5 = prctile(vv,100);
G = zeros(length(var),1);
thresh = [0 t1 t2 t3 t4 t5];
for i = 1:5
    G(vv > thresh(i) & vv <= thresh(i+1))= i;
end
figure(4)
subplot(2,4,2)
boxchart(G,var,'MarkerStyle','.');
ylim([0 50])
xticklabels({'Q1','Q2','Q3','Q4','Q5'})
xlabel('Reservoir Capacity')
xticks([1 2 3 4 5])
box on
set(gca,'FontSize',14)
ylabel('NRMSE (% Capacity)')


vv = gvar;
t1 = prctile(vv,20);
t2 = prctile(vv,40);
t3 = prctile(vv,60);
t4 = prctile(vv,80);
t5 = prctile(vv,100);
G = zeros(length(var),1);
thresh = [0 t1 t2 t3 t4 t5];
for i = 1:5
    G(vv > thresh(i) & vv <= thresh(i+1))= i;
end
figure(4)
subplot(2,4,1)
boxchart(G,var,'MarkerStyle','.');
ylim([0 100])
xticklabels({'Q1','Q2','Q3','Q4','Q5'})
xlabel('Mean Seasonal Variability')
xticks([1 2 3 4 5])
box on
set(gca,'FontSize',14)
ylabel('NRMSE (% Capacity)')


vv = area;
t1 = prctile(vv,20);
t2 = prctile(vv,40);
t3 = prctile(vv,60);
t4 = prctile(vv,80);
t5 = prctile(vv,100);
G = zeros(length(var),1);
thresh = [0 t1 t2 t3 t4 t5];
for i = 1:5
    G(vv > thresh(i) & vv <= thresh(i+1))= i;
end
figure(4)
subplot(2,4,3)
boxchart(G,var,'MarkerStyle','.');
ylim([0 50])
xticklabels({'Q1','Q2','Q3','Q4','Q5'})

xlabel('Reservoir Surface Area')
xticks([1 2 3 4 5])
box on
set(gca,'FontSize',14)
ylabel('NRMSE (% Capacity)')

vv = abs(lat);
t1 = prctile(vv,20);
t2 = prctile(vv,40);
t3 = prctile(vv,60);
t4 = prctile(vv,80);
t5 = prctile(vv,100);
G = zeros(length(var),1);
thresh = [0 t1 t2 t3 t4 t5];
for i = 1:5
    G(vv > thresh(i) & vv <= thresh(i+1))= i;
end
figure(4)
subplot(2,4,4)
boxchart(G,var,'MarkerStyle','.');
ylim([0 50])
xticklabels({'Q1','Q2','Q3','Q4','Q5'})
xlabel('Latitude (Absolute Value)')
xticks([1 2 3 4 5])
box on
set(gca,'FontSize',14)
ylabel('NRMSE (% Capacity)')

vv = sdi;
t1 = prctile(vv,20);
t2 = prctile(vv,40);
t3 = prctile(vv,60);
t4 = prctile(vv,80);
t5 = prctile(vv,100);
G = zeros(length(var),1);
thresh = [0 t1 t2 t3 t4 t5];
for i = 1:5
    G(vv > thresh(i) & vv <= thresh(i+1))= i;
end
figure(4)
subplot(2,4,5)
boxchart(G,var,'MarkerStyle','.');
ylim([0 50])
xticklabels({'Q1','Q2','Q3','Q4','Q5'})
xticks([1 2 3 4 5])
xlabel('Shoreline Development Index')
box on
set(gca,'FontSize',14)
ylabel('NRMSE (% Capacity)')

vv = cloud;
t1 = prctile(vv,20);
t2 = prctile(vv,40);
t3 = prctile(vv,60);
t4 = prctile(vv,80);
t5 = prctile(vv,100);
G = zeros(length(var),1);
thresh = [0 t1 t2 t3 t4 t5];
for i = 1:5
    G(vv > thresh(i) & vv <= thresh(i+1))= i;
end
figure(4)
subplot(2,4,6)
boxchart(G,var,'MarkerStyle','.');
ylim([0 50])
xticklabels({'Q1','Q2','Q3','Q4','Q5'})

xlabel('Mean Cloudiness')
xticks([1 2 3 4 5])
box on
set(gca,'FontSize',14)
ylabel('NRMSE (% Capacity)')

vv = tri(~isnan(tri));
t1 = prctile(vv,20);
t2 = prctile(vv,40);
t3 = prctile(vv,60);
t4 = prctile(vv,80);
t5 = prctile(vv,100);
G = zeros(length(vv),1);
thresh = [0 t1 t2 t3 t4 t5];
for i = 1:5
    G(vv > thresh(i) & vv <= thresh(i+1))= i;
end
figure(4)
subplot(2,4,7)
boxchart(G,var(~isnan(tri)),'MarkerStyle','.');
ylim([0 50])
xticklabels({'Q1','Q2','Q3','Q4','Q5'})

xlabel('Terrain Ruggedness Index')
xticks([1 2 3 4 5])
box on
set(gca,'FontSize',14)
ylabel('NRMSE (% Capacity)')



vv = year(year > 0);
t1 = prctile(vv,20);
t2 = prctile(vv,40);
t3 = prctile(vv,60);
t4 = prctile(vv,80);
t5 = prctile(vv,100);
G = zeros(length(vv),1);
thresh = [0 t1 t2 t3 t4 t5];
for i = 1:5
    G(vv > thresh(i) & vv <= thresh(i+1))= i;
end
figure(4)
subplot(2,4,8)
boxchart(G,var(year > 0),'MarkerStyle','.');
ylim([0 50])
xticklabels({'Q1','Q2','Q3','Q4','Q5'})

xlabel('Year Built')
xticks([1 2 3 4 5])
box on
set(gca,'FontSize',14)
ylabel('NRMSE (% Capacity)')





%% zscores
%calculate z-scores - based on mean variability 

rm = 100*[compare_data.rmse]'./[compare_data.rv_mcm]';
rm = nanmean(rm,2);
seavar = 100*[compare_data.mean_sea_var]'./[compare_data.rv_mcm]';
thresh(1,1) = 0;
for j = 1:10
    thresh(j+1,1) = prctile(seavar,j*10);
end
for i = 1:length(rm)
    for j = 1:10
        if seavar(i) > thresh(j) && seavar(i) <= thresh(j+1)
            group(i,1) = j;
        end
    end
end

grpmeans = splitapply(@median,rm,group);
grpstds = splitapply(@std,rm,group);

for i = 1:length(rm)
    zsvar(i,1) = (rm(i) - grpmeans(group(i)))./grpstds(group(i));
    zmvar(i,1) = (rm(i) - grpmeans(group(i)))./grpmeans(group(i));
end
var = zsvar;


vv = vol;
t1 = prctile(vv,20);
t2 = prctile(vv,40);
t3 = prctile(vv,60);
t4 = prctile(vv,80);
t5 = prctile(vv,100);
G = zeros(length(var),1);
thresh = [0 t1 t2 t3 t4 t5];
for i = 1:5
    G(vv > thresh(i) & vv <= thresh(i+1))= i;
end
figure(5)
subplot(2,4,2)
boxchart(G,var,'MarkerStyle','.');
xticklabels({'Q1','Q2','Q3','Q4','Q5'})
xlabel('Reservoir Capacity')
xticks([1 2 3 4 5])
box on
set(gca,'FontSize',14)
ylabel('NRMSE Z-Score (by Var)')
ylim([-1.5 2])


vv = gvar;
t1 = prctile(vv,20);
t2 = prctile(vv,40);
t3 = prctile(vv,60);
t4 = prctile(vv,80);
t5 = prctile(vv,100);
G = zeros(length(var),1);
thresh = [0 t1 t2 t3 t4 t5];
for i = 1:5
    G(vv > thresh(i) & vv <= thresh(i+1))= i;
end
figure(5)
subplot(2,4,1)
boxchart(G,var,'MarkerStyle','.');
xticklabels({'Q1','Q2','Q3','Q4','Q5'})
xlabel('Mean Seasonal Variability')
xticks([1 2 3 4 5])
box on
set(gca,'FontSize',14)
ylabel('NRMSE Z-Score (by Var)')
ylim([-1.5 2])


vv = area;
t1 = prctile(vv,20);
t2 = prctile(vv,40);
t3 = prctile(vv,60);
t4 = prctile(vv,80);
t5 = prctile(vv,100);
G = zeros(length(var),1);
thresh = [0 t1 t2 t3 t4 t5];
for i = 1:5
    G(vv > thresh(i) & vv <= thresh(i+1))= i;
end
figure(5)
subplot(2,4,3)
boxchart(G,var,'MarkerStyle','.');
xticklabels({'Q1','Q2','Q3','Q4','Q5'})
ylabel('NRMSE Z-Score (by Var)')
xticks([1 2 3 4 5])
box on
set(gca,'FontSize',14)
xlabel('Reservoir Surface Area')
ylim([-1.5 2])

vv = abs(lat);
t1 = prctile(vv,20);
t2 = prctile(vv,40);
t3 = prctile(vv,60);
t4 = prctile(vv,80);
t5 = prctile(vv,100);
G = zeros(length(var),1);
thresh = [0 t1 t2 t3 t4 t5];
for i = 1:5
    G(vv > thresh(i) & vv <= thresh(i+1))= i;
end
figure(5)
subplot(2,4,4)
boxchart(G,var,'MarkerStyle','.');
xticklabels({'Q1','Q2','Q3','Q4','Q5'})
xlabel('Latitude (Absolute Value)')
xticks([1 2 3 4 5])
box on
set(gca,'FontSize',14)
ylabel('NRMSE Z-Score (by Var)')
ylim([-1.5 2])



vv = sdi;
t1 = prctile(vv,20);
t2 = prctile(vv,40);
t3 = prctile(vv,60);
t4 = prctile(vv,80);
t5 = prctile(vv,100);
G = zeros(length(var),1);
thresh = [0 t1 t2 t3 t4 t5];
for i = 1:5
    G(vv > thresh(i) & vv <= thresh(i+1))= i;
end
figure(5)
subplot(2,4,5)
boxchart(G,var,'MarkerStyle','.');
xticklabels({'Q1','Q2','Q3','Q4','Q5'})
xticks([1 2 3 4 5])
xlabel('Shoreline Development Index')
box on
set(gca,'FontSize',14)
ylabel('NRMSE Z-Score (by Var)')
ylim([-1.5 2])

vv = cloud;
t1 = prctile(vv,20);
t2 = prctile(vv,40);
t3 = prctile(vv,60);
t4 = prctile(vv,80);
t5 = prctile(vv,100);
G = zeros(length(var),1);
thresh = [0 t1 t2 t3 t4 t5];
for i = 1:5
    G(vv > thresh(i) & vv <= thresh(i+1))= i;
end
figure(5)
subplot(2,4,6)
boxchart(G,var,'MarkerStyle','.');
xticklabels({'Q1','Q2','Q3','Q4','Q5'})

xlabel('Mean Cloudiness')
xticks([1 2 3 4 5])
box on
set(gca,'FontSize',14)
ylabel('NRMSE Z-Score (by Var)')
ylim([-1.5 2])

vv = tri(~isnan(tri));
t1 = prctile(vv,20);
t2 = prctile(vv,40);
t3 = prctile(vv,60);
t4 = prctile(vv,80);
t5 = prctile(vv,100);
G = zeros(length(vv),1);
thresh = [0 t1 t2 t3 t4 t5];
for i = 1:5
    G(vv > thresh(i) & vv <= thresh(i+1))= i;
end
figure(5)
subplot(2,4,7)
boxchart(G,var(~isnan(tri)),'MarkerStyle','.');
xticklabels({'Q1','Q2','Q3','Q4','Q5'})

xlabel('Terrain Ruggedness Index')
xticks([1 2 3 4 5])
box on
set(gca,'FontSize',14)
ylabel('NRMSE Z-Score (by Var)')
ylim([-1.5 2])



vv = year(year > 0);
t1 = prctile(vv,20);
t2 = prctile(vv,40);
t3 = prctile(vv,60);
t4 = prctile(vv,80);
t5 = prctile(vv,100);
G = zeros(length(vv),1);
thresh = [0 t1 t2 t3 t4 t5];
for i = 1:5
    G(vv > thresh(i) & vv <= thresh(i+1))= i;
end
figure(5)
subplot(2,4,8)
boxchart(G,var(year > 0),'MarkerStyle','.');
xticklabels({'Q1','Q2','Q3','Q4','Q5'})

xlabel('Year Built')
xticks([1 2 3 4 5])
box on
set(gca,'FontSize',14)
ylabel('NRMSE Z-Score (by Var)')
ylim([-1.5 2])



































