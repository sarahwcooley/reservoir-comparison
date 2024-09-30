cd('/Users/sc961/University of Oregon Dropbox/Sarah Cooley/Reservoir_Review/datasets');
%cd('/Users/sarahcooley/Dropbox (University of Oregon)/Reservoir_Review/datasets');
load('complete_dataset_mar22.mat');

c = out_data; %(type == 1);
year = [out_data.year]';
c(year < 1999) = [];
clear valid
map = {'#e3b505','#95190c','#610345','#107e7d','#044b7f'};
cmap = validatecolor(map, 'multiple');
colororder(cmap)

for i = 1:length(c)
    type(i,1) = get_type(c(i).valid);
end
%c = c(type == 1);
for i = 1:length(c)
    cts = c(i).trend_ts;
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
validv = valid;
clear o


% for j = 1:5
%     tr = trend(:,j); %(type == 1,j);
%     ptr = ptrend(:,j); % (type == 1,j);
%    o(1,j) = 100*nansum(ptr < 0.05)./sum(~isnan(ptr));
%    o(2,j) = 100*nansum(ptr < 0.05 & tr > 0)./sum(~isnan(ptr));
%    o(3,j) = 100*nansum(ptr < 0.05 & tr < 0)./sum(~isnan(ptr));
%    o(4,j) = 100*nansum(ptr > 0.05)./sum(~isnan(ptr));
%    
% 
% end


count = 6;
nn = unique(continent);
%nn(1) = [];
%nn(3) = [];
clear o

titles = {'Africa','Asia','North America','Oceania','South America','Europe'};
for p = 1:6
    t = find(continent == nn(p));

glwsts = 0.001*nansum(glws(:,t),2);
grsts = 0.001*nansum(grs(:,t),2);
glots = 0.001*nansum(glo(:,t),2);
grdlyts = 0.001*nansum(grdly(:,t),2);
grdllts = 0.001*nansum(grdll(:,t),2);
validout = valid(t,:);
sumvalid = sum(validout);
[tr,ptr] = get_trend(glwsts,months,1);
if ptr > 0.05; tr = 0; end
o(p).glwst = tr;
o(p).glwsp = ptr;
[tr,ptr] = get_trend(grsts,months,1);
if ptr > 0.05; tr = 0; end
o(p).grst = tr;
o(p).grsp = ptr;
[tr,ptr] = get_trend(glots,months,1);
if ptr > 0.05; tr = 0; end
o(p).glot = tr;
o(p).glop = ptr;
[tr,ptr] = get_trend(grdlyts,months,1);
if ptr > 0.05; tr = 0; end
o(p).grdlyt = tr;
o(p).grdlyp = ptr;
[tr,ptr] = get_trend(grdllts,months,1);
if ptr > 0.05; tr = 0; end
o(p).grdllt = tr;
o(p).grdllp = ptr;
o(p).num_res = length(t);

figure(3)
subplot(7,3,count);
hold off
plot(months,glwsts,'Color',cmap(1,:),'LineWidth',1.5);
hold on
plot(months,grsts,'Color',cmap(2,:),'LineWidth',1.5);
plot(months,glots,'Color',cmap(3,:),'LineWidth',1.5);
plot(months,grdlyts,'Color',cmap(4,:),'LineWidth',1.5);
plot(months,grdllts,'Color',cmap(5,:),'LineWidth',1.5);
%legend(['GLWS (n=' num2str(sumvalid(1)) ')'],['GRS (n=' num2str(sumvalid(2)) ')'],['GloLakes (n=' num2str(sumvalid(3)) ')'],...
   % ['GRDL-Y (n=' num2str(sumvalid(4)) ')'],['GRDL-L (n=' num2str(sumvalid(5)) ')']);
set(gca,'FontSize',12);
%ylabel('Volume (km^3)');
count = count + 3;
%title(titles{p});

end

sumvalid = sum(valid);
glwsts = 0.001*nansum(glws,2);
grsts = 0.001*nansum(grs,2);
glots = 0.001*nansum(glo,2);
grdlyts = 0.001*nansum(grdly,2);
grdllts = 0.001*nansum(grdll,2);
figure(3)
subplot(7,3,3)
hold off
plot(months,glwsts,'Color',cmap(1,:),'LineWidth',1.5);
hold on
plot(months,grsts,'Color',cmap(2,:),'LineWidth',1.5);
plot(months,glots,'Color',cmap(3,:),'LineWidth',1.5);
plot(months,grdlyts,'Color',cmap(4,:),'LineWidth',1.5);
plot(months,grdllts,'Color',cmap(5,:),'LineWidth',1.5);
%legend(['GLWS (n=' num2str(sumvalid(1)) ')'],['GRS (n=' num2str(sumvalid(2)) ')'],['GloLakes (n=' num2str(sumvalid(3)) ')'],...
    %['GRDL-Y (n=' num2str(sumvalid(4)) ')'],['GRDL-L (n=' num2str(sumvalid(5)) ')']);
%ylabel('Volume (km^3)');
%title('Global')
set(gca,'FontSize',12);


%% part 2
clear all
cd('/Users/scooley2/Dropbox (University of Oregon)/Reservoir_Review/datasets');
%cd('/Users/sarahcooley/Dropbox (University of Oregon)/Reservoir_Review/datasets');
load('complete_dataset_mar22.mat');

c = out_data; %(type == 1);
year = [out_data.year]';
c(year >= 1999) = [];

clear valid
map = {'#e3b505','#95190c','#610345','#107e7d','#044b7f'};
cmap = validatecolor(map, 'multiple');
colororder(cmap)

for i = 1:length(c)
    type(i,1) = get_type(c(i).valid);
end
for i = 1:length(c)
    cts = c(i).trend_ts;
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
validv = valid;
clear o




count = 5;
nn = unique(continent);
nn(1) = [];
%nn(3) = [];
clear o

titles = {'Africa','Asia','North America','Oceania','South America','Europe'};
for p = 1:6
    t = find(continent == nn(p));

glwsts = 0.001*nansum(glws(:,t),2);
grsts = 0.001*nansum(grs(:,t),2);
glots = 0.001*nansum(glo(:,t),2);
grdlyts = 0.001*nansum(grdly(:,t),2);
grdllts = 0.001*nansum(grdll(:,t),2);
validout = valid(t,:);
sumvalid = sum(validout);
[tr,ptr] = get_trend(glwsts,months,1);
if ptr > 0.05; tr = 0; end
o(p).glwst = tr;
o(p).glwsp = ptr;
[tr,ptr] = get_trend(grsts,months,1);
if ptr > 0.05; tr = 0; end
o(p).grst = tr;
o(p).grsp = ptr;
[tr,ptr] = get_trend(glots,months,1);
if ptr > 0.05; tr = 0; end
o(p).glot = tr;
o(p).glop = ptr;
[tr,ptr] = get_trend(grdlyts,months,1);
if ptr > 0.05; tr = 0; end
o(p).grdlyt = tr;
o(p).grdlyp = ptr;
[tr,ptr] = get_trend(grdllts,months,1);
if ptr > 0.05; tr = 0; end
o(p).grdllt = tr;
o(p).grdllp = ptr;
o(p).num_res = length(t);

figure(3)
subplot(7,3,count);
hold off
plot(months,glwsts,'Color',cmap(1,:),'LineWidth',1.5);
hold on
plot(months,grsts,'Color',cmap(2,:),'LineWidth',1.5);
plot(months,glots,'Color',cmap(3,:),'LineWidth',1.5);
plot(months,grdlyts,'Color',cmap(4,:),'LineWidth',1.5);
plot(months,grdllts,'Color',cmap(5,:),'LineWidth',1.5);
%legend(['GLWS (n=' num2str(sumvalid(1)) ')'],['GRS (n=' num2str(sumvalid(2)) ')'],['GloLakes (n=' num2str(sumvalid(3)) ')'],...
%    ['GRDL-Y (n=' num2str(sumvalid(4)) ')'],['GRDL-L (n=' num2str(sumvalid(5)) ')']);
set(gca,'FontSize',12);
%ylabel('Volume (km^3)');
count = count + 3;
%title(titles{p});

end

sumvalid = sum(valid);
glwsts = 0.001*nansum(glws,2);
grsts = 0.001*nansum(grs,2);
glots = 0.001*nansum(glo,2);
grdlyts = 0.001*nansum(grdly,2);
grdllts = 0.001*nansum(grdll,2);
figure(3)
subplot(7,3,2)
hold off
plot(months,glwsts,'Color',cmap(1,:),'LineWidth',1.5);
hold on
plot(months,grsts,'Color',cmap(2,:),'LineWidth',1.5);
plot(months,glots,'Color',cmap(3,:),'LineWidth',1.5);
plot(months,grdlyts,'Color',cmap(4,:),'LineWidth',1.5);
plot(months,grdllts,'Color',cmap(5,:),'LineWidth',1.5);
%legend(['GLWS (n=' num2str(sumvalid(1)) ')'],['GRS (n=' num2str(sumvalid(2)) ')'],['GloLakes (n=' num2str(sumvalid(3)) ')'],...
%    ['GRDL-Y (n=' num2str(sumvalid(4)) ')'],['GRDL-L (n=' num2str(sumvalid(5)) ')']);
%ylabel('Volume (km^3)');
%title('Global')
set(gca,'FontSize',12);

%% part 3
clear all
cd('/Users/scooley2/Dropbox (University of Oregon)/Reservoir_Review/datasets');
%cd('/Users/sarahcooley/Dropbox (University of Oregon)/Reservoir_Review/datasets');
load('complete_dataset_mar22.mat');

c = out_data; %(type == 1);
year = [out_data.year]';

clear valid
map = {'#e3b505','#95190c','#610345','#107e7d','#044b7f'};
cmap = validatecolor(map, 'multiple');
colororder(cmap)

for i = 1:length(c)
    type(i,1) = get_type(c(i).valid);
end
for i = 1:length(c)
    cts = c(i).trend_ts;
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
validv = valid;
clear o




count = 4;
nn = unique(continent);
nn(1) = [];
%nn(3) = [];
clear o

titles = {'Africa','Asia','North America','Oceania','South America','Europe'};
for p = 1:6
    t = find(continent == nn(p));

glwsts = 0.001*nansum(glws(:,t),2);
grsts = 0.001*nansum(grs(:,t),2);
glots = 0.001*nansum(glo(:,t),2);
grdlyts = 0.001*nansum(grdly(:,t),2);
grdllts = 0.001*nansum(grdll(:,t),2);
validout = valid(t,:);
sumvalid = sum(validout);
[tr,ptr] = get_trend(glwsts,months,1);
if ptr > 0.05; tr = 0; end
o(p).glwst = tr;
o(p).glwsp = ptr;
[tr,ptr] = get_trend(grsts,months,1);
if ptr > 0.05; tr = 0; end
o(p).grst = tr;
o(p).grsp = ptr;
[tr,ptr] = get_trend(glots,months,1);
if ptr > 0.05; tr = 0; end
o(p).glot = tr;
o(p).glop = ptr;
[tr,ptr] = get_trend(grdlyts,months,1);
if ptr > 0.05; tr = 0; end
o(p).grdlyt = tr;
o(p).grdlyp = ptr;
[tr,ptr] = get_trend(grdllts,months,1);
if ptr > 0.05; tr = 0; end
o(p).grdllt = tr;
o(p).grdllp = ptr;
o(p).num_res = length(t);

figure(3)
subplot(7,3,count);
hold off
plot(months,glwsts,'Color',cmap(1,:),'LineWidth',1.5);
hold on
plot(months,grsts,'Color',cmap(2,:),'LineWidth',1.5);
plot(months,glots,'Color',cmap(3,:),'LineWidth',1.5);
plot(months,grdlyts,'Color',cmap(4,:),'LineWidth',1.5);
plot(months,grdllts,'Color',cmap(5,:),'LineWidth',1.5);
%legend(['GLWS (n=' num2str(sumvalid(1)) ')'],['GRS (n=' num2str(sumvalid(2)) ')'],['GloLakes (n=' num2str(sumvalid(3)) ')'],...
%    ['GRDL-Y (n=' num2str(sumvalid(4)) ')'],['GRDL-L (n=' num2str(sumvalid(5)) ')']);
set(gca,'FontSize',12);
%ylabel('Volume (km^3)');
count = count + 3;
%title(titles{p});

end

sumvalid = sum(valid);
glwsts = 0.001*nansum(glws,2);
grsts = 0.001*nansum(grs,2);
glots = 0.001*nansum(glo,2);
grdlyts = 0.001*nansum(grdly,2);
grdllts = 0.001*nansum(grdll,2);
figure(3)
subplot(7,3,1)
hold off
plot(months,glwsts,'Color',cmap(1,:),'LineWidth',1.5);
hold on
plot(months,grsts,'Color',cmap(2,:),'LineWidth',1.5);
plot(months,glots,'Color',cmap(3,:),'LineWidth',1.5);
plot(months,grdlyts,'Color',cmap(4,:),'LineWidth',1.5);
plot(months,grdllts,'Color',cmap(5,:),'LineWidth',1.5);
%legend(['GLWS (n=' num2str(sumvalid(1)) ')'],['GRS (n=' num2str(sumvalid(2)) ')'],['GloLakes (n=' num2str(sumvalid(3)) ')'],...
%    ['GRDL-Y (n=' num2str(sumvalid(4)) ')'],['GRDL-L (n=' num2str(sumvalid(5)) ')']);
%ylabel('Volume (km^3)');
%title('Global')
set(gca,'FontSize',12);

