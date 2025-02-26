cd('/Users/sc961/University of Oregon Dropbox/Sarah Cooley/Reservoir_Review/datasets');
load('validation_dataset_feb5');
ids = [all_data.final_id]';


cd('/Users/sc961/University of Oregon Dropbox/Sarah Cooley/Reservoir_Review/datasets');
load('complete_dataset_feb12.mat');

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
    oo(i).type = type(i,1);
    oo(i).cts = cts;
    oo(i).gts = gauge_ts;
    oo(i).months = months;
    oo(i).ats = comp_out(i).ats;
    oo(i).name = all_data(idgauge).fname;
    oo(i).vol = comp_out(i).rv_mcm;
    oo(i).lat = latout;
    oo(i).lon = lonout;
    oo(i).grand_id = comp_out(i).grand_id;
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

rm = rmse_norm(type == 1,:);
rma = rmse_norm_abs(type == 1,:);

grp1 = repmat(1:5,size(rm,1),1);
grp2 = repmat(1:5,size(rma,1),1);
clr1 = repmat(1,size(rm));
clr2 = repmat(2,size(rma));

y = [rm;rma];
x = [grp1;grp2];
c = [clr1;clr2];
%%
%x = x(:);
%y = y(:);
%c = c(:);
%x=x*1.25;
figure(8)
map = {'#e3b505','#95190c','#610345','#107e7d','#044b7f'};
cmap = validatecolor(map, 'multiple');
colororder(cmap)


subplot(6,2,[1 3])
hold off
for i = 1:5
    if i == 1; hold off; else; hold on; end
b1 = boxchart(grp1(:,i),rm(:,i),'MarkerStyle','.','BoxWidth',0.5);
b1.BoxFaceColor = cmap(i,:);
end
hold on
ylim([ 0 80])
ylabel('RMSE (% of Capacity)');
xticks([1 2 3 4 5])
xticklabels({'GLWS','GRS','GloLakes','GRDL-Y','GRDL-L'});
%title([num2str(sum(type == 1)) ' Reservoirs'])
set(gca,'FontSize',14)
box on
xlim([0.2 5.8])

subplot(6,2,[2 4])
for i = 1:5
    if i == 1; hold off; else; hold on; end
b1 = boxchart(grp1(:,i)-0.2,rm(:,i),'MarkerStyle','.','BoxWidth',0.15);
b1.BoxFaceColor = cmap(i,:);
end

for i = 1:5
    b1 = boxchart(grp2(:,i)+0.2,rma(:,i),'MarkerStyle','.','BoxWidth',0.3);
    b1.BoxFaceColor = cmap(i,:);
end
hold on
ylim([ 0 140])
ylabel('RMSE (% of Capacity)');
xticks([1 2 3 4 5])
xticklabels({'GLWS','GRS','GloLakes','GRDL-Y','GRDL-L'});
%title([num2str(sum(type == 1)) ' Reservoirs'])
set(gca,'FontSize',14)
box on
xlim([0.2 5.8])




tt = [oo.type]';
op = oo;
op(tt ~= 1) = [];

ids = [132
451
2725
4946
5014
4776
6629
616];
%indids = [2 3 5 6 8 9 11 12];
indids = [5 7 9 11];
gids = [op.grand_id]';
for j =1:4
    figure(8)
    subplot(6,2,indids(j));
    gg = find(gids == ids(j));
    months = op(gg).months;
    gauge = op(gg).gts;
    gauge = gauge - gauge(1);
    gauge = gauge; %./op(gg).vol;
    cts = op(gg).cts;
    for i = 1:5
        if i == 1; hold off; else; hold on; end
        plot(months,cts(:,i),'LineWidth',1,'Color',cmap(i,:));
    end
    plot(months,gauge,'LineWidth',1,'Color','k');
    %title(op(ids(j)).name)
    %ylabel('Storage Change (% Capacity)')
    if j == 4; ylabel('Storage Change (% Capacity)'); end
    if j == 8; xlabel('Year'); end
set(gca,'FontSize',14)
end
indids = [6 8 10 12];
for j =1:4
    figure(8)
    subplot(6,2,indids(j));
    gg = find(gids == ids(j+4));
    months = op(gg).months;
    gauge = op(gg).gts;
    %gauge = gauge - gauge(1);
    gauge = gauge; %./op(gg).vol;
    cts = op(gg).ats;
    for i = 1:4
        if i == 1; hold off; else; hold on; end
        plot(months,cts(:,i),'LineWidth',1,'Color',cmap(i+1,:));
    end
    plot(months,gauge,'LineWidth',1,'Color','k');
    %title(op(ids(j)).name)
    %ylabel('Storage Change (% Capacity)')
    if j == 4; ylabel('Storage Change (% Capacity)'); end
    if j == 8; xlabel('Year'); end
set(gca,'FontSize',14)
end



