%cd('/Users/scooley2/Dropbox (University of Oregon)/Reservoir_Review/datasets');
cd('/Users/sc961/University of Oregon Dropbox/Sarah Cooley/Reservoir_Review/datasets');
load('complete_dataset_feb12.mat');

%1. Relative statistics, for ALL reservoirs, 473 reservoirs, and by dataset
%(RMSE, RMSE %, Seasonal Var, Seasonal Var %, R2, MBE,
%2. Absolute statistics, for ALL reservoirs, 473 reservoirs, and by dataset
%3. Correlation/predictor analyses, for ALL reservoirs, 473 reservoirs, and
%by dataset
%(both corr and partial corr controlling for XX)


%statistical analyses testing whether year is significant


%% construct dataset
c = 1;
for i = 1:length(out_data)
    
    cts = out_data(i).cts;
    ats = out_data(i).ats;
    valid = out_data(i).valid;
    vol = out_data(i).rv_mcm;
    if sum(valid) > 1
        
        %RMSE, R2, MBE
        [rmse,r2,mbe] = get_rmse_comparison(cts,valid);
        if sum(valid(2:end))>1
            [rmsea,r2a,mbea] = get_rmse_comparison(ats,valid(2:end));
        else
            rmsea = nan(5,1);
            r2a = nan(5,1);
            mbea = nan(5,1);
        end
        %Seasonal Var
        sea_var = out_data(i).sea_var;
        [sea_var_abs,sea_var] = get_comparison_values_sed_r(sea_var);

        %type
        type = get_type(valid);

        o(c).grand_id = out_data(i).grand_id;
        o(c).type = type;
        o(c).lat = out_data(i).lat;
        o(c).lon = out_data(i).lon;
        o(c).vol = out_data(i).rv_mcm;
        o(c).year = out_data(i).year;
        o(c).area = out_data(i).area;
        o(c).valid = valid;
        o(c).mean_elv = out_data(i).mean_elv;
        o(c).std_elv = out_data(i).std_elv;
        o(c).tri_elv = out_data(i).tri_elv;
        o(c).cloud = out_data(i).cloud;
        o(c).sdibc = out_data(i).sdibc;
        o(c).meansvar_norm = 100*nanmean(out_data(i).sea_var)./out_data(i).rv_mcm;
        o(c).trend = out_data(i).trend;
        o(c).ptrend = out_data(i).ptrend;
        o(c).continent = out_data(i).continent;
        o(c).basin = out_data(i).basin;
        o(c).rmse = rmse;
        o(c).rmse_norm = 100*rmse./vol;
        o(c).rmsea = rmsea;
        o(c).rmsea_norm = 100*rmsea./vol;
        o(c).r2 = r2;
        o(c).mbea = mbea;
        o(c).mbea_norm = 100*mbea./vol;
        o(c).sea_var = sea_var;
        o(c).sea_var_norm = 100*sea_var./vol;
        c = c+1;
    end
end

clear out_data

%% Statistical Analyses for Relative, 473 reservoirs

type = [o.type]';
co = o(type == 1);
table1(1,:) = median([co.rmse]');
table1(2,:) = median([co.rmse_norm]');
table1(3,:) = median([co.sea_var]');
table1(4,:) = median([co.sea_var_norm]');
table1(5,:) = median([co.r2]');
table1(6,:) = prctile([co.rmse_norm]',25);
table1(7,:) = prctile([co.rmse_norm]',75);


table2(1,:) = median([co.rmsea]');
table2(2,:) = median([co.rmsea_norm]');
table2(3,:) = median([co.mbea]');
table2(4,:) = median([co.mbea_norm]');
table2(5,:) = prctile([co.rmsea_norm]',25);
table2(6,:) = prctile([co.rmsea_norm]',75);

%% Bar Chart for 


figure(1)
map = {'#e3b505','#95190c','#610345','#107e7d','#044b7f'};
cmap = validatecolor(map, 'multiple');
colororder(cmap)
grp1 = repmat(1:5,473,1,1);

subplot(2,2,1)
hold off
var = [co.rmse_norm]';
for i = 1:5
    if i == 1; hold off; else; hold on; end
    b1 = boxchart(grp1(:,i),var(:,i),'MarkerStyle','.','BoxWidth',0.5);
    b1.BoxFaceColor = cmap(i,:);
end
hold on
ylim([ 0 100])
ylabel('Relative Change RMSE (% of Capacity)');
xticks([1 2 3 4 5])
xticklabels({'GLWS','GRS','GloLakes','GRDL-Y','GRDL-L'});
set(gca,'FontSize',14)
box on


subplot(2,2,2)
hold off

var = [co.sea_var_norm]';
for i = 1:5
    if i == 1; hold off; else; hold on; end
    b1 = boxchart(grp1(:,i),var(:,i),'MarkerStyle','.','BoxWidth',0.5);
    b1.BoxFaceColor = cmap(i,:);
end
hold on
ylim([ -50 50])
ylabel('Seasonal Variability Error (%)');
xticks([1 2 3 4 5])
xticklabels({'GLWS','GRS','GloLakes','GRDL-Y','GRDL-L'});
set(gca,'FontSize',14)
box on
    

subplot(2,2,3)
hold off
var = [co.rmsea_norm]';
var = cat(2,nan(473,1),var);
for i = 1:5
    if i == 1; hold off; else; hold on; end
    b1 = boxchart(grp1(:,i),var(:,i),'MarkerStyle','.','BoxWidth',0.5);
    b1.BoxFaceColor = cmap(i,:);
end
hold on
ylim([ 0 150])
ylabel('Absolute RMSE (% of Capacity)');
xticks([1 2 3 4 5])
xticklabels({'GLWS','GRS','GloLakes','GRDL-Y','GRDL-L'});
set(gca,'FontSize',14)
box on


subplot(2,2,4)
hold off
var = [co.mbea_norm]';
var = cat(2,nan(473,1),var);
for i = 1:5
    if i == 1; hold off; else; hold on; end
    b1 = boxchart(grp1(:,i),var(:,i),'MarkerStyle','.','BoxWidth',0.5);
    b1.BoxFaceColor = cmap(i,:);
end
hold on
ylim([-150 150])
ylabel('Absolute MBE (% of Capacity)');
xticks([1 2 3 4 5])
xticklabels({'GLWS','GRS','GloLakes','GRDL-Y','GRDL-L'});
set(gca,'FontSize',14)
box on


%% Statistics across 473 reservoirs and ALL reservoirs
rmse_norm = [co.rmse_norm]';
rmse_norm = nanmean(rmse_norm,2);
rmsea_norm = [co.rmsea_norm]';
rmsea_norm = nanmean(rmsea_norm,2);
table3(1,1) = median(rmse_norm);
table3(1,2) = prctile(rmse_norm,25)';
table3(1,3) = prctile(rmse_norm,75)';
table3(1,4) = mean(rmse_norm);
table3(2,1) = nanmedian(rmsea_norm);
table3(2,2) = prctile(rmsea_norm,25)';
table3(2,3) = prctile(rmsea_norm,75)';
table3(2,4) = mean(rmsea_norm);
rmse_norm = [o.rmse_norm]';
rmse_norm = nanmean(rmse_norm,2);
rmsea_norm = [o.rmsea_norm]';
rmsea_norm = nanmean(rmsea_norm,2);
table3(3,1) = nanmedian(rmse_norm);
table3(3,2) = prctile(rmse_norm,25)';
table3(3,3) = prctile(rmse_norm,75)';
table3(3,4) = mean(rmse_norm);
table3(4,1) = nanmedian(rmsea_norm);
table3(4,2) = prctile(rmsea_norm,25)';
table3(4,3) = prctile(rmsea_norm,75)';
table3(4,4) = nanmean(rmsea_norm);

%% Predictor Analyses
%starting with 473, individual to each dataset
%predictor vars = sea_var, cloud, elev, sdibc, vol, area, year

%first, individual correlations
rmse = [co.rmsea_norm]';
svar = [co.meansvar_norm]';
cloud = [co.cloud]';
elev = [co.tri_elv]';
sdi = [co.sdibc]';
year = [co.year]';
vol = log([co.vol]');
area = log([co.area]');
lat = abs([co.lat]');
testvars = [svar,cloud,elev,sdi,year,vol,area,lat];
for i = 1:5
    var = rmse(:,i);
    for j = 1:8
        tvar = testvars(:,j);
        [r,p] = corrcoef(var(~isnan(tvar)),tvar(~isnan(tvar)));
        table4(j,i*2-1) = r(2);
        table4(j,i*2) = p(2);
    end
end

%next, cross correlations
for i = 1:5
    var = rmse(:,i);
    tvar = [var,testvars];
    [rho,pval] = partialcorr(tvar);
    table5(:,i*2-1) = rho(:,1);
    table5(:,i*2) = pval(:,1);
end

%first, individual correlations
rmse = [co.rmsea_norm]';
for i = 1:4
    var = rmse(:,i);
    for j = 1:8
        tvar = testvars(:,j);
        [r,p] = corrcoef(var(~isnan(tvar)),tvar(~isnan(tvar)));
        table6(j,i*2-1) = r(2);
        table6(j,i*2) = p(2);
    end
end

%next, cross correlations
for i = 1:4
    var = rmse(:,i);
    tvar = [var,testvars];
    [rho,pval] = partialcorr(tvar);
    table7(:,i*2-1) = rho(:,1);
    table7(:,i*2) = pval(:,1);
end

%% correlations, for ALL data
rmse = [o.rmse_norm]';
svar = [o.meansvar_norm]';
cloud = [o.cloud]';
elev = [o.tri_elv]';
sdi = [o.sdibc]';
year = [o.year]';
vol = log([o.vol]');
area = log([o.area]');
lat = abs([o.lat]');
testvars = [svar,cloud,elev,sdi,year,vol,area,lat];
var = nanmean(rmse,2);
testvars(isnan(var),:) = [];
var(isnan(var),:) = [];
for j = 1:8
   tvar = testvars(:,j);
   [r,p] = corrcoef(var(~isnan(tvar)),tvar(~isnan(tvar)));
   table8(j,1) = r(2);
   table8(j,2) = p(2);
end

%next, cross correlations
for i = 1:7
    ind = find(isnan(testvars(:,i)));
    if i == 1;
        indout = ind;
    else
        indout = cat(1,indout,ind);
    end
end
var(indout,:) = [];
testvars(indout,:) = [];
tvar = [var,testvars];
[rho,pval] = partialcorr(tvar);
table8(:,3) = rho(2:end,1);
table8(:,4) = pval(2:end,1);


rm1 = rmse(year < 1999);
rm2 = rmse(year >= 1999);
[h,p] = ttest2(rm1,rm2);
table9(1,1) = nanmedian(rm1);
table9(1,2) = prctile(rm1,25);
table9(1,3) = prctile(rm1,75);
table9(1,4) = nanmean(rm1);
table9(2,1) = nanmedian(rm2);
table9(2,2) = prctile(rm2,25);
table9(2,3) = prctile(rm2,75);
table9(2,4) = nanmean(rm2);

