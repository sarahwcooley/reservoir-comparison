cd('/Users/scooley2/Dropbox (University of Oregon)/Reservoir_Review/datasets');
%cd('/Users/sarahcooley/Dropbox (University of Oregon)/Reservoir_Review/datasets');
load('complete_dataset_mar22.mat');

for i = 1:length(out_data)
    type(i,1) = get_type(out_data(i).valid);
end

c = out_data(type == 1);
for i = 1:length(c)
    cts(:,:,i) = c(i).cts;
    basin(:,i) = c(i).basin;
end


cd('/Users/scooley2/Dropbox (University of Oregon)/Reservoir_Review/datasets/ancillary_data/untitled folder');
basins = shaperead('hydrobasins_dec21.shp');
basins(152) = [];

%% figure 1a seasonal variability 
clear shapeplot
for i = 1:length(basins)
    ind = find(basin == basins(i).id);
    shapeplot(i).Y = basins(i).Y;
    shapeplot(i).X = basins(i).X;
    if length(ind) >= 2
        ctsout = cts(:,:,ind);
        ctsout = 0.001*nansum(ctsout,3);
        valid = [1 1 1 1 1];
        remove = [0 0 0 0 0];
        [out] = get_complete_comparison_analysis(ctsout,months,valid,remove);
        trend = out.trend;
        trendp = out.trendp;
        trend_type = get_trend_type(trend,trendp,[1 1 1 1 1]);
        shapeplot(i).trend_type = trend_type;     
    else
        shapeplot(i).trend_type = -1;   
    end
    shapeplot(i).Geometry = 'Polygon';
    shapeplot(i).id = i;
end
cd('/Users/scooley2/Dropbox (University of Oregon)/Reservoir_Review/datasets/basin_shapefiles/mar25')
shapewrite(shapeplot,'trend_type_basins_overlap_mar25.shp');





%% all basins not intersecting

c = out_data;
clear valid cts
for i = 1:length(c)
   cts(:,:,i) = c(i).cts;
   basin(:,i) = c(i).basin;
end
ids = [c.grand_id]';
ii = length(basin)+1;
all = 1;
if all == 1
    cd('/Users/scooley2/Dropbox (University of Oregon)/Reservoir_Review/datasets');

    load('li_dataset_v2.mat');
    load('li_spatial_data.mat');
    li_ids = [li_data.grand_id]';
    lia = ismember(li_ids,ids);
    li_data(lia) = [];
    cont_li = spatial.basin;
    cont_li(lia) = [];
    ts = NaN(240,5);
    for j = 1:length(li_data)
        basin(ii,1) = cont_li(j);
        ts(:,2) = 1000*(li_data(j).out_ts - li_data(j).out_ts(1));
        cts(:,:,ii) = ts;
        valid(ii,:) = [0 1 0 0 0];
        ii = ii+1;
    end
end


cmap = lines(5);
for i = 1:length(c)
    type(i,1) = get_type(c(i).valid);
end
clear basin cts
for i = 1:length(c)
    cts(:,:,i) = c(i).cts;
    basin(:,i) = c(i).basin;
end
cd('/Users/scooley2/Dropbox (University of Oregon)/Reservoir_Review/datasets/ancillary_data/untitled folder');
basins = shaperead('hydrobasins_dec21.shp');
basins(152) = [];
clear shapeplot
for i = 1:length(basins)
    ind = find(basin == basins(i).id);
    shapeplot(i).Y = basins(i).Y;
    shapeplot(i).X = basins(i).X;
    if length(ind) >= 2
        ctsout = cts(:,:,ind);
        
        ctsout = 0.001*nansum(ctsout,3);
        vv = sum(ctsout);
        valid = [1 1 1 1 1];
        valid(vv == 0) = 0;
        remove = [0 0 0 0 0];
        [out] = get_complete_comparison_analysis(ctsout,months,valid,remove);
        trend = out.trend;
        trendp = out.trendp;
        valid(2) = [];
        trend_type = get_trend_type(trend,trendp,valid);
        shapeplot(i).trend_type = trend_type;     
    else
        shapeplot(i).trend_type = -1;   
    end
    shapeplot(i).Geometry = 'Polygon';
    shapeplot(i).id = i;
end
cd('/Users/scooley2/Dropbox (University of Oregon)/Reservoir_Review/datasets/basin_shapefiles')
shapewrite(shapeplot,'trend_type_basins_all_mar12.shp');