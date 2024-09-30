cd('/Users/scooley2/Dropbox (University of Oregon)/Reservoir_Review/datasets');
%cd('/Users/sarahcooley/Dropbox (University of Oregon)/Reservoir_Review/datasets');
load('complete_dataset_mar22.mat');
%load('complete_dataset_feb13.mat','compare_data');


for i = 1:length(out_data)
    type(i,1) = get_type(out_data(i).valid);
end
c = out_data;
clear valid
for i = 1:length(c)
    cts = c(i).cts;
    ctsout(:,:,i) = cts;
    glws(:,i) = cts(:,1);
    grs(:,i) = cts(:,2);
    glo(:,i) = cts(:,3);
    grdly(:,i) = cts(:,4);
    grdll(:,i) = cts(:,5);
    vol(i,1) = c(i).rv_mcm;
    valid(:,i) = c(i).valid;
    continent(i,1) = c(i).continent;
    basin(i,1) = c(i).basin;
    trend(i,:) = c(i).trend;
    ptrend(i,:) = c(i).ptrend;
    type(i,1) = get_type(c(i).valid);
end
ii = length(continent)+1;
valid = valid';
allv = 1;
if allv == 1
    load('li_dataset_v2.mat');
    load('li_spatial_data.mat');
    li_ids = [li_data.grand_id]';
    lia = ismember(li_ids,ids);
    li_data(lia) = [];
    cont_li = spatial.continent;
    bas_li = spatial.basin;
    cont_li(lia) = [];
    for j = 1:length(li_data)
        continent(ii,1) = cont_li(j);
        glws(:,ii) = NaN;
        grs(:,ii) = 1000*(li_data(j).out_ts - li_data(j).out_ts(1));
        glo(:,ii) = NaN;
        grdly(:,ii) = NaN;
        grdll(:,ii) = NaN;
        vol(ii,1) = li_data(j).rv_mcm;
        valid(ii,:) = [0 1 0 0 0];
        basin(ii,1) = bas_li(j);
        ts = 1000*(li_data(j).out_ts - li_data(j).out_ts(1));
         [m,w] = size(ts);
         if w > m; ts = ts'; end
        [tr,ptr] = get_trend(ts,months);
        trend(ii,:) = [NaN tr NaN NaN NaN];
        ptrend(ii,:) = [NaN ptr NaN NaN NaN];
        type(ii,1) = 5;
        ctsout(:,:,ii) = [glws(:,ii) grs(:,ii) glo(:,ii) grdly(:,ii) grdll(:,ii)];
        ii = ii+1;
    end
end


%% all

cd('/Users/scooley2/Dropbox (University of Oregon)/Reservoir_Review/datasets/ancillary_data/untitled folder');
basins = shaperead('hydrobasins_dec21.shp');
basins(152) = [];
clear shapeplot
c1 = 1;
c2 = 1;
c3 = 1;
c4 = 1;
for i = 1:length(basins)
    ind = find(basin == basins(i).id);
    if length(ind) >= 3
        ctsout1 = ctsout(:,:,ind);
        voli = vol(ind);
        for j = 1:5
            ids = ~isnan(ctsout1(10,j,:));
            sumvol(1,j) = 0.001*sum(voli(ids(1,1,:)));
        end

        ctsout1 = 0.001*nansum(ctsout1,3);
        %ctsout1 = 100*ctsout1./sumvol;
        valid = [1 1 1 1 1 ];
        remove = [0 0 0 0 0];
        [out] = get_complete_comparison_analysis(ctsout1,months,valid,remove);
        trend = out.trend;
        trendp = out.trendp;
        trend_type = get_trend_type(trend,trendp,[1 1 1 1 1]);
        shapeplot(i).trend_type = trend_type;     
        if trend_type == 1 || trend_type == 3
            trout = mean(trend(trendp < 0.05));
            shapehigh(c1).Y = basins(i).Y;
            shapehigh(c1).X = basins(i).X;
            shapehigh(c1).Geometry = 'Polygon';
            shapehigh(c1).trout_a = trout;
            shapehigh(c1).id = i;
            shapehigh(c1).num_res = length(ind);
            c1 = c1 + 1;
        end
        if trend_type == 2 || trend_type == 4
            if trend_type == 2; trout = mean(trend(trend < 0 & trendp < 0.05)); end
            if trend_type == 4; trout = mean(trend(trend > 0 & trendp < 0.05)); end
            shapelow(c2).Y = basins(i).Y;
            shapelow(c2).X = basins(i).X;
            shapelow(c2).Geometry = 'Polygon';
            shapelow(c2).trout = trout;
            c2 = c2 + 1;
        end
        if trend_type > 4
            shapenone(c3).Y = basins(i).Y;
            shapenone(c3).X = basins(i).X;
            shapenone(c3).Geometry = 'Polygon';
            shapenone(c3).type = trend_type;
            c3 = c3+1;
        end
    else
    shapenodata(c4).Y = basins(i).Y;
    shapenodata(c4).X = basins(i).X;
    shapenodata(c4).Geometry = 'Polygon';
    shapenodata(c4).id = i;
    c4 = c4+1;
    end
end
cd('/Users/scooley2/Dropbox (University of Oregon)/Reservoir_Review/datasets/basin_shapefiles/mar19')
shapewrite(shapehigh,'trend_type_basins_all_high_val_mar18.shp');
%shapewrite(shapelow,'trend_type_basins_all_low_mar18.shp');
%shapewrite(shapenodata,'trend_type_basins_all_nodata_mar18.shp');
%shapewrite(shapenone,'trend_type_basins_all_none_mar18.shp');

%% plot figure
