cd('/Users/scooley2/Dropbox (University of Oregon)/Reservoir_Review/gauge_data/');
load('complete_comp_v2.mat');
ids = [all_data.final_id]';


cd('/Users/scooley2/Dropbox (University of Oregon)/Reservoir_Review/datasets');
load('complete_dataset_mar22.mat');
for i=  1:length(out_data)
    out_data(i).ats(:,1) = out_data(i).ts.grs;
    if out_data(i).valid(3) == 1
    out_data(i).ats(:,2) = out_data(i).ts.glolakes(138:377);
    end
end
comp_data = out_data;

grids = [comp_data.grand_id]';
lia = ismember(grids,ids);
comp_out = comp_data(lia);
count = 1;
for i = 1:length(comp_out)
    type(i,1) = get_type(comp_out(i).valid);
    ind = comp_out(i).grand_id;
    idgauge = find(ids == ind);
    gauge_ts = all_data(idgauge).data;
    months = all_data(idgauge).months;
    cts = comp_out(i).ats;
    [rmse(i,:), sea_err(i,:), sea_err_abs(i,:), tot_err(i,:), ann_err(i,:), ann_err_abs(i,:),r2val(i,:),pval(i,:)] = get_gauge_error(gauge_ts, cts, months,0);
    vol(i,1) = comp_out(i).rv_mcm;
    latout = comp_out(i).lat;
    lonout = comp_out(i).lon;
    tout(i,1) = all_data(idgauge).type;
end
rmse_norm_abs = 100*rmse./vol;
sea_erra = sea_err;
rmsea = rmse;
ann_erra = ann_err;
ann_err_absa = ann_err_abs;
sea_err_absa = sea_err_abs;
clear rmse sea_err sea_err_abs tot_err ann_err ann_err_abs r2val pval
for i = 1:length(comp_out)
    type(i,1) = get_type(comp_out(i).valid);
    ind = comp_out(i).grand_id;
    idgauge = find(ids == ind);
    gauge_ts = all_data(idgauge).data;
    months = all_data(idgauge).months;
    cts = comp_out(i).cts;
    [rmse(i,:), sea_err(i,:), sea_err_abs(i,:), tot_err(i,:), ann_err(i,:), ann_err_abs(i,:),r2val(i,:),pval(i,:)] = get_gauge_error(gauge_ts, cts, months,1);
    vol(i,1) = comp_out(i).rv_mcm;
    latout = comp_out(i).lat;
    lonout = comp_out(i).lon;
    tout(i,1) = all_data(idgauge).type;
end
rmse_norm = 100*rmse./vol;
rmse_norm_abs = cat(2,nan(size(rmse_norm_abs(:,1))),rmse_norm_abs);




rm = rmse_norm_abs(type == 1,:);
%rm = 100*rm./vol(type == 1);
o(1,:) = median(rm);
se = sea_err(type == 1,:);
se = 100*se./vol(type == 1);
o(2,:) = median(se);
o(3,:) = median(abs(se));
te = tot_err(type == 1,:);
te = 100*te./vol(type == 1);
o(4,:) = median(te);
o(5,:) = median(abs(te));
an = ann_err(type == 1,:);
an = 100*an./vol(type == 1);
o(6,:) = median(an);
o(7,:) = median(abs(an));
o(8,:) = median(r2val(type ==1,:));

clear o
rm = rmse_norm_abs(:,[2 3 5]);
%rm = rmse_norm(:,[1 2 4]);

ns = sum(isnan(rm),2);
rm = rmse_norm_abs(ns == 0,:);
%rm = 100*rm./vol(ns == 0,:);
o(1,:) = median(rm);
se = sea_err(ns == 0,:);
se = 100*se./vol(ns == 0,:);
o(2,:) = median(se);
o(3,:) = median(abs(se));
te = tot_err(ns == 0,:);
te = 100*te./vol(ns == 0,:);
o(4,:) = median(te);
o(5,:) = median(abs(te));
an = ann_err(ns == 0,:);
an = 100*an./vol(ns == 0,:);
o(6,:) = median(an);
o(7,:) = median(abs(an));
o(8,:) = median(r2val(ns ==0,:));

%% figure 1
figure(1)

se = sea_err(type == 1,:);
se = 100*se./vol(type == 1);
subplot(2,2,1)
hold off
%yline(0)
%hold on
boxchart(se,'MarkerStyle','.');
ylim([-50 50])
ylabel('Seasonal Error (% of Capacity)')
xticklabels({'GLWS','GRS','GloLakes','GRDL-Y','GRDL-L'});
set(gca,'FontSize',14)

se = sea_err_abs(type == 1,:);
se = 100*se./vol(type == 1);
subplot(2,2,2)
boxchart(se,'MarkerStyle','.');
ylim([ 0 50])
ylabel('Seasonal Error (% of Capacity)')
xticklabels({'GLWS','GRS','GloLakes','GRDL-Y','GRDL-L'});
set(gca,'FontSize',14)


an = ann_err(type == 1,:);
an = 100*an./vol(type == 1);
subplot(2,2,3)
hold off
%yline(0)
hold on
boxchart(an,'MarkerStyle','.');
ylim([-50 50])
ylabel('Annual Mean Error (% of Capacity)')
xticklabels({'GLWS','GRS','GloLakes','GRDL-Y','GRDL-L'});
set(gca,'FontSize',14)

an = ann_err_abs(type == 1,:);
an = 100*an./vol(type == 1);
subplot(2,2,4)
boxchart(an,'MarkerStyle','.');
ylim([ 0 50])
ylabel('Annual Mean Error (% of Capacity)')
xticklabels({'GLWS','GRS','GloLakes','GRDL-Y','GRDL-L'});
set(gca,'FontSize',14)



%% figure 3
rm = rmse_norm(type == 1,:);
rma = rmse_norm_abs(type == 1,:);

grp1 = repmat(1:5,size(rm,1),1);
grp2 = repmat(1:5,size(rma,1),1);
clr1 = repmat(1,size(rm));
clr2 = repmat(2,size(rma));

y = [rm;rma];
x = [grp1;grp2];
c = [clr1;clr2];

x = x(:);
y = y(:);
c = c(:);
x=x*1.25;
figure(8)
subplot(1,3,1)
boxchart(x,y,'GroupByColor',c,'MarkerStyle','.');
ylim([ 0 150])
ylabel('RMSE (% of Capacity)');
xticks(1.25*[1 2 3 4 5])
xticklabels({'GLWS','GRS','GloLakes','GRDL-Y','GRDL-L'});
title([num2str(sum(type == 1)) ' Reservoirs'])
set(gca,'FontSize',14)
box on


rm = rmse_norm(:,[2 3 5]);
rma = rmse_norm_abs(:,[2 3 5]);
ns = sum(isnan(rm),2);
rm = rm(ns == 0,:);
rma = rma(ns == 0,:);


grp1 = repmat(1:3,size(rm,1),1);
grp2 = repmat(1:3,size(rma,1),1);
clr1 = repmat(1,size(rm));
clr2 = repmat(2,size(rma));

y = [rm;rma];
x = [grp1;grp2];
c = [clr1;clr2];
x = x(:);
y = y(:);
c = c(:);
x=x*1.25;


figure(8)
subplot(1,3,2)
boxchart(x,y,'GroupByColor',c,'MarkerStyle','.');
ylim([ 0 150])
ylabel('RMSE (% of Capacity)');
xticks(1.25*[1 2 3])
xticklabels({'GRS','GloLakes','GRDL-L'});
title([num2str(sum(sum(ns == 0))) ' Reservoirs'])
set(gca,'FontSize',14)
box on




rm = rmse_norm(:,[2  5]);
rma = rmse_norm_abs(:,[2  5]);
ns = sum(isnan(rm),2);
rm = rm(ns == 0,:);
rma = rma(ns == 0,:);

grp1 = repmat(1:2,size(rm,1),1);
grp2 = repmat(1:2,size(rma,1),1);
clr1 = repmat(1,size(rm));
clr2 = repmat(2,size(rma));

y = [rm;rma];
x = [grp1;grp2];
c = [clr1;clr2];
x = x(:);
y = y(:);
c = c(:);
x=x*1.25;

figure(8)
subplot(1,3,3)
boxchart(x,y,'GroupByColor',c,'MarkerStyle','.');
ylim([ 0 150])
ylabel('RMSE (% of Capacity)');
xticks(1.25*[1 2 ])
xticklabels({'GRS','GRDL-L'});
title([num2str(sum(sum(ns == 0))) ' Reservoirs'])
set(gca,'FontSize',14)
box on
legend('Relative','Absolute')



%% 

figure(3)
subplot(2,3,1)
rm = rmse_norm;
boxchart(rm(type == 1,:),'MarkerStyle','.');
ylim([ 0 50])
ylabel('RMSE (% of Capacity)');
xticklabels({'GLWS','GRS','GloLakes','GRDL-Y','GRDL-L'});
title([num2str(sum(type == 1)) ' Reservoirs'])
set(gca,'FontSize',14)
box on


subplot(2,3,2)
rm = rmse_norm(:,[2 3 5]);
ns = sum(isnan(rm),2);
boxchart(rm(ns == 0,:),'MarkerStyle','.');
ylim([ 0 50])
ylabel('RMSE (% of Capacity)');
xticklabels({'GRS','GloLakes','GRDL-L'});
title([num2str(sum(ns == 0)) ' Reservoirs'])
set(gca,'FontSize',14)
box on


subplot(2,3,3)
rm = rmse_norm(:,[2  5]);
ns = sum(isnan(rm),2);
boxchart(rm(ns == 0,:),'MarkerStyle','.');
ylim([ 0 50])
ylabel('RMSE (% of Capacity)');
xticklabels({'GRS','GRDL-L'});
title([num2str(sum(ns == 0)) ' Reservoirs'])
set(gca,'FontSize',14)
box on

subplot(2,3,4)
rm = rmse_norm_abs;
boxchart(rm(type == 1,:),'MarkerStyle','.');
ylim([ 0 100])
ylabel('RMSE (% of Capacity)');
xticklabels({'GRS','GloLakes','GRDL-Y','GRDL-L'});
title([num2str(sum(type == 1)) ' Reservoirs'])
set(gca,'FontSize',14)
box on


subplot(2,3,5)
rm = rmse_norm_abs(:,[1 2 4]);
ns = sum(isnan(rm),2);
boxchart(rm(ns == 0,:),'MarkerStyle','.');
ylim([ 0 100])
ylabel('RMSE (% of Capacity)');
xticklabels({'GRS','GloLakes','GRDL-L'});
title([num2str(sum(ns == 0)) ' Reservoirs'])
set(gca,'FontSize',14)
box on


subplot(2,3,6)
rm = rmse_norm_abs(:,[1  4]);
ns = sum(isnan(rm),2);
boxchart(rm(ns == 0,:),'MarkerStyle','.');
ylim([ 0 100])
ylabel('RMSE (% of Capacity)');
xticklabels({'GRS','GRDL-L'});
title([num2str(sum(ns == 0)) ' Reservoirs'])
set(gca,'FontSize',14)
box on

%% analysis of causes of error

%get predictor variables
cd('/Users/scooley2/Dropbox (University of Oregon)/Reservoir_Review/datasets');
%load('complete_dataset_mar22.mat');
gids = [out_data.grand_id]';

for i = 1:length(comp_out)
    ind = find(gids == comp_out(i).grand_id);
    tri_elv(i,1) = out_data(ind).tri_elv;
    sdibc(i,1) = 0.001*out_data(ind).sdibc;
    cloud(i,1) =out_data(ind).cloud;
    area(i,1) = out_data(ind).area;
    ind = comp_out(i).grand_id;
    idgauge = find(ids == ind);
    gauge_ts = all_data(idgauge).data;
    m = months;
    m(isnan(gauge_ts)) = [];
    gauge_ts(isnan(gauge_ts)) = [];
    yy = year(m);
 
        
    G = findgroups(yy);
    annual_max = splitapply(@max,gauge_ts',G);
    annual_min = splitapply(@min,gauge_ts',G);
    var = annual_max - annual_min;            
    vol(i,1) = comp_out(i).rv_mcm;
    gvar(i,1) = nanmean(var);
    tvar(i,1) = max(gauge_ts) - min(gauge_ts);
end
sdif = 100*gvar./tvar;
sea = nanmean(sea_err,2);
rm = nanmean(rmse_norm,2);
gvar = 100*gvar./vol;
ann = 100*nanmean(ann_err,2)./vol;
tot = 100*nanmean(abs(tot_err),2)./vol;


sea = 100*sea./vol;
%rm = 100*rm./vol;
%rm = sea;

[r,p] = corrcoef(tri_elv,rm);
op.trielvr = r(2);
op.trielvp = p(2);

[r,p] = corrcoef(sdibc,rm);
op.stibcr = r(2);
op.stibcp = p(2);
[r,p] = corrcoef(gvar,rm);
op.gvarr = r(2);
op.gvarp = p(2);
[r,p] = corrcoef(log(area),rm);
op.arear = r(2);
op.areavp = p(2);
[r,p] = corrcoef(log(vol),rm);
op.volr = r(2);
op.volp = p(2);
[r,p] = corrcoef(cloud,rm);
op.cloudr = r(2);
op.cloudp = p(2);
[r,p] = corrcoef(sdif,rm);
op.sdifr = r(2);
op.sdifp = p(2);

figure(5)
subplot(2,3,1)
cmap = lines(5);
c = cmap(1,:);
semilogx(vol,rm,'o','MarkerFaceColor',c,'MarkerEdgeColor',c);
box on
xlabel('Capacity (MCM)');
ylabel('RMSE (% Capacity)')
set(gca,'FontSize',14)
subplot(2,3,2)
semilogx(area,rm,'o','MarkerFaceColor',c,'MarkerEdgeColor',c);
box on
xlabel('Area (km^2)');
ylabel('RMSE (% Capacity)')
set(gca,'FontSize',14)

subplot(2,3,3)
scatter(gvar,rm,'o','MarkerFaceColor',c,'MarkerEdgeColor',c);
box on
xlabel('Mean Seasonal Variability (% Capacity)');
ylabel('RMSE (% Capacity)')
set(gca,'FontSize',14)

subplot(2,3,4)
scatter(sdibc,rm,'MarkerFaceColor',c,'MarkerEdgeColor',c);
box on
xlabel('Shoreline Development Index');
ylabel('RMSE (% Capacity)');
set(gca,'FontSize',14)

subplot(2,3,5)
scatter(tri_elv,rm,'MarkerFaceColor',c,'MarkerEdgeColor',c);
box on
xlabel('Terrain Ruggedness Index');
ylabel('RMSE (% Capacity)');
set(gca,'FontSize',14)

subplot(2,3,6)
scatter(cloud,rm,'MarkerFaceColor',c,'MarkerEdgeColor',c);
box on
xlabel('Cloudiness (%)');
ylabel('RMSE (% Capacity)');
set(gca,'FontSize',14)


rm = tot;
rr(1,1) = mean(rm);
rr(2,1) = median(rm);
rr(3,1) = prctile(rm,25);
rr(4,1) = prctile(rm,75);
rr(5,1) = 100*sum(abs(rm) < 5)./length(rm);
rr(6,1) = 100*sum(abs(rm) < 10)./length(rm);
rr(7,1) = 100*sum(abs(rm) < 25)./length(rm);
rr(8,1) = 100*sum(abs(rm) < 50)./length(rm);


vol = [comp_out.rv_mcm]';
area = [comp_out.area]';
gvar = [comp_out.sea_var]';
gvar = nanmean(gvar,2);

for i = 1:4;
    if i == 1; var = area; end
    if i == 2; var = vol; end
    if i == 3; var = gvar; end
    if i == 4; var = 100*gvar./vol; end
rp(1,i) = mean(var);
rp(2,i) = median(var);
rp(3,i) = max(var);
rp(4,i) = min(var);
end



