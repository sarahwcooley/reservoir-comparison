cd('/Users/scooley2/Dropbox (University of Oregon)/Reservoir_Review/datasets');
%cd('/Users/sarahcooley/Dropbox (University of Oregon)/Reservoir_Review/datasets');
load('complete_dataset_mar22.mat','compare_data');



for i = 1:length(compare_data)
    type(i,1) = get_type(compare_data(i).valid);
    rmse(i,:) = compare_data(i).rmse;
    vol(i,1) = compare_data(i).rv_mcm;
end


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



cd('/Users/scooley2/Library/CloudStorage/GoogleDrive-cooleysarahw@gmail.com/My Drive/000_Proposals/2023/NASA_SWOT/figures/World_Countries_Generalized');
wc = shaperead('World_Countries_Generalized.shp');
lat = [compare_data.lat]';
lon = [compare_data.lon]';
vert = [lon, lat];

for i = 1:length(wc)
    edge = [wc(i).X', wc(i).Y'];
    stat = inpoly2(vert,edge);
    if sum(stat) == 0
        stat = zeros(length(vert),1);
    end
    
    if i == 1 
        country = i*stat;
    else
        country = country + i*stat;
    end
end
count = 1;
for i = 1:length(wc)
    ind = find(country == i);
    if length(ind) >= 3
        cc(count).id = i;
        cc(count).name = wc(i).COUNTRY;
        cc(count).num_res = length(ind);
        cc(count).med_rmse = median(rm(ind));
        cc(count).mean_rmse = mean(rm(ind));
        cc(count).per_10r = 100*sum(rm(ind) < 10)./length(ind);
        cc(count).per_50r = 100*sum(rm(ind) > 50)./length(ind);
        cc(count).med_zsvar = median(zsvar(ind));
        cc(count).mean_zsvar = mean(zsvar(ind));
        cc(count).per_0z = 100*sum(zsvar(ind) < 0)./length(ind);
        cc(count).X = wc(i).X;
        cc(count).Y = wc(i).Y;
        cc(count).Geometry = 'Polygon';
        count = count + 1;
    end
end

nm = [cc.num_res]';
cp = cc;
cp(nm < 10) = [];
co = cc;
co(nm < 5) = [];

cd('/Users/scooley2/Dropbox (University of Oregon)/Reservoir_Review/datasets/basin_shapefiles/mar25');
shapewrite(cc,'country_data_v1_mar28.shp');
shapewrite(cp,'country_data_v1_mar28_large.shp');
shapewrite(co,'country_data_v1_mar28_5min.shp');

cpp = cc(nm >= 20);

for i = 1:length(cpp)
    cpp(i).hdi = hdi(i);
end
cpp(22) = [];
hdi(22) = [];
[hdis,B] = sort(hdi,'ascend');


ids = [cpp.id]';
lia = ismember(country,ids);
zvs = zsvar(lia);
cids = country(lia);

ids = ids(B);
sorted_ids = zeros(length(cids),1);
for i = 1:length(ids)
    ind = find(cids == ids(i));
    sorted_ids(ind) = i;
end

figure(2)
boxchart(sorted_ids,zvs)
ylim([-2 2])



per0 = [cpp.per_0z]';
per10 = [cpp.per_10r]';
zs = [cpp.med_zsvar]';
mp = [cpp.meanlat]';
%[hdis,B] = sort(hdi,'descend');
rms = [cpp.med_rmse]';


dd = zeros(101,1);
dd(data(:,4) < 60) = 1;
figure(1)
subplot(1,2,1)
hold off
scatter(data(dd == 1,1),data(dd == 1,3),70,'Filled');
hold on
scatter(data(dd == 0,1),data(dd == 0,3),70,'Filled');
box on
set(gca,'FontSize',18)
xlabel('Human Development Index')
ylabel('Median NRMSE (%)')
subplot(1,2,2)
hold off
scatter(data(dd == 1,1),data(dd == 1,2),70,'Filled');
hold on
scatter(data(dd == 0,1),data(dd == 0,2),70,'Filled');
box on
set(gca,'FontSize',18)
xlabel('Human Development Index')
ylabel('Median NRMSE Z-Score (var-normalized)')
legend({'Below 60°N','Above 60°N'})



hdi = data(:,1);
G = zeros(size(hdi));
thresh = [0 0.55 0.7 0.8 1];
for j = 1:length(thresh)-1
    ind = find(hdi > thresh(j) & hdi <= thresh(j+1));
    G(ind) = j;
end

figure(2)
subplot(2,1,1)
boxchart(G,data(:,3),'MarkerStyle','.')
subplot(2,1,2)
boxchart(G,(data(:,2)),'MarkerStyle','.')
set(gca,'FontSize',18)


G(data(:,4) >= 60) = 5;
figure(3)
subplot(2,1,1)
boxchart(G,data(:,3),'MarkerStyle','.','MarkerSize',15)
xticks([ 1 2 3 4 5])

xticklabels({'Low','Medium','High','Very High','Very High'});
xlabel('Human Development Index')
ylabel('Median NRMSE (% Capacity)')
ylim([0 40])
set(gca,'FontSize',18)
box on
xlim([0.35 5.65])
subplot(2,1,2)
boxchart(G,(data(:,2)),'MarkerStyle','+');
set(gca,'FontSize',18)
xticks([ 1 2 3 4 5])
xticklabels({'Low','Medium','High','Very High','Very High'});
xlabel('Human Development Index')
ylabel('Median NRMSE Z-Score (var-normalized)')
set(gca,'FontSize',18)
ylim([-0.7 1])
box on
xlim([0.35 5.65])







