clear all; close all; clc;
imt = 3; ngm = 7;
im = 144; jm = 90; 
Tlat = [90:-2:-90]; %  2deg lat cell
Tlon = [-180:2.5:180];%  2.5deg long cell
[Tlg,Tlt]=meshgrid(Tlon,Tlat);
N_hemisph = zeros(im,jm); N_hemisph(1:im,jm/2+1:jm) = 1;
S_hemisph = zeros(im,jm); S_hemisph(1:im,1:jm/2) = 1;
% i_site = 38; j_site = 65; %MMSF
% i_site=35; j_site=84; %bug at this cell
i_site=44; j_site=29; %irrigation bug at this cell
% i_site = 82; j_site = 76; %Hyytiala
% i_site = 24; j_site = 65; %Vaira
% i_site = 43; j_site = 66; %NYC
% i_site = 51; j_site = 44; %TNF
% i_site = 87; j_site = 46; %Mpala
%--------------------------------------------------------------------------
% Read in fraction of cells covered by land
ncload ('/NoBackup/modelE_out/JAN1902.ijEA2c20irrig.nc','frac_land')
fearth=rot90(flipud(frac_land),3);

% Read in area of grid cell
fid = fopen('/Users/mpuma/Documents/LSM_2X2_5/area_2x2_5.bi'...
    ,'r','ieee-be.l64');
if fid == -1
    disp(message)
end
skip = fread(fid,[4]);
grid_area = fread(fid,[im,jm],'float64','ieee-be.l64'); % in km2
skip = fread(fid,[4]);
fclose(fid);
%--------------------------------------------------------------------------

land_area=grid_area.*fearth;
earth_land=sum(sum(land_area));
%earth_land = 130,580,000 km ice-free land;%148,940,000 km land
%--------------------------------------------------------------------------
fid = fopen('/Users/mpuma/Documents/LSM_2X2_5/CD144X90.ext',...
    'r','ieee-be');
if fid == -1
    disp(message)
end


skip = fread(fid,[4]);
CD_title = fread(fid,[80],'char','ieee-be');
temp = fread(fid,[im jm],'float32','ieee-be'); 
z0m= 30./10.^(temp);
site_z0m=z0m(i_site,j_site);
skip = fread(fid,[4]);

fclose(fid);
figure(5)
z0m(i_site,j_site) = 10000;
var_plot=z0m;
var_plot(find(fearth==0))=NaN;
TsK = rot90(var_plot);
TsK1 = TsK;
TsK1(1:jm,im+1) = TsK(1:jm,im);
TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
m_proj('Equidistant Cylindrical','lon',0,'lat',90);
% set the figure background to white
set(gcf,'color','white')
m_pcolor(Tlg,Tlt,TsK1);
shading flat; % gets rid of the grid box edges
colormap(jet);%colormap(flipud(HSV));
m_coast('color',[0.3 0.3 0.3]);
m_grid('box','fancy','color',[0 0 0]);
h=colorbar('EastOutside');
title('Roughness length for momentum, z0m [m]')
%--------------------------------------------------------------------------
fid = fopen('/Users/mpuma/Documents/LSM_2X2_5/V144X90_no_crops.ext',...
    'r','ieee-be');
if fid == -1
    disp(message)
end

% Read in 10 Matthews vegetation types
for i = 1:10
    skip = fread(fid,[4]);
    veg_title = fread(fid,[80],'char','ieee-be');
    temp = fread(fid,[im jm],'float32','ieee-be'); 
    site_veg(i)=temp(i_site,j_site);
    temp(find(fearth==0))=NaN;
    VEG(:,:,i) = temp;
    skip = fread(fid,[4]);
end
fclose(fid);
%--------------------------------------------------------------------------
fid = fopen('/Users/mpuma/Documents/LSM_2X2_5/w_soil_max.bin',...
    'r','ieee-be.l64');
if fid == -1
    disp(message)
end

for k = 1:6
    for ibv =1:2 %1 - bare; 2 - vegetated
       skip = fread(fid,[4]);
       temp = fread(fid,[im jm],'float64','ieee-be.l64'); 
       w_soil_max(k,ibv,:,:) = temp;
       skip = fread(fid,[4]);
    end
end
fclose(fid);
%--------------------------------------------------------------------------

% fid = fopen('/Users/mpuma/Documents/LSM_2X2_5/CROPS_144X90N.nocasp.ext',...
%     'r','ieee-be');
% if fid == -1
%     disp(message)
% end
% % Read in CROP cover (24 time periods: 1700, 1750, 1800, 1850:10:1980;
% % 1986:1:1992
% for i = 1:1
%     skip = fread(fid,[4]);
%     CROPS_title = fread(fid,[80],'char','ieee-be');
%     temp = fread(fid,[im jm],'float32','ieee-be'); 
%     %temp(find(fearth==0))=NaN;
%     CROPS(:,:,i) = temp;
%     skip = fread(fid,[4]);
% end
% fclose(fid);
% puma=CROPS(:,:,1);
%--------------------------------------------------------------------------
% site data
veg_vector = zeros(10,1);
veg_vector(:,1) = VEG(i_site,j_site,:);
bugcell=fearth(i_site,j_site);
%--------------------------------------------------------------------------
figure(1)
subplot(2,1,1)
var_plot=top_index;
var_plot(find(fearth==0))=NaN;
TsK = rot90(var_plot);
TsK1 = TsK;
TsK1(1:jm,im+1) = TsK(1:jm,im);
TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
m_proj('Equidistant Cylindrical','lon',0,'lat',90);
% set the figure background to white
set(gcf,'color','white')
m_pcolor(Tlg,Tlt,TsK1);
shading flat; % gets rid of the grid box edges
colormap(jet);%colormap(flipud(HSV));
m_coast('color',[0.3 0.3 0.3]);
m_grid('box','fancy','color',[0 0 0]);
h=colorbar('EastOutside');
title('Topographic Index')

subplot(2,1,2)
var_plot=top_dev;
var_plot(find(fearth==0))=NaN;
TsK = rot90(var_plot);
TsK1 = TsK;
TsK1(1:jm,im+1) = TsK(1:jm,im);
TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
m_proj('Equidistant Cylindrical','lon',0,'lat',90);
% set the figure background to white
set(gcf,'color','white')
m_pcolor(Tlg,Tlt,TsK1);
shading flat; % gets rid of the grid box edges
colormap(jet);%colormap(flipud(HSV));
m_coast('color',[0.3 0.3 0.3]);
m_grid('box','fancy','color',[0 0 0]);
h=colorbar('EastOutside');
title('Standard Deviation of Topographic Index')
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
figure(3)
var_plot=sl;
var_plot(find(fearth==0))=NaN;
TsK = rot90(var_plot);
TsK1 = TsK;
TsK1(1:jm,im+1) = TsK(1:jm,im);
TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
m_proj('Equidistant Cylindrical','lon',0,'lat',90);
% set the figure background to white
set(gcf,'color','white')
m_pcolor(Tlg,Tlt,TsK1);
shading flat; % gets rid of the grid box edges
colormap(jet);%colormap(flipud(HSV));
m_coast('color',[0.3 0.3 0.3]);
m_grid('box','fancy','color',[0 0 0]);
h=colorbar('EastOutside');
title('sl: slope???')

%--------------------------------------------------------------------------
figure(4)
subplot(3,2,1)
var_plot=dz(:,:,1);
var_plot(find(fearth==0))=NaN;
TsK = rot90(var_plot);
TsK1 = TsK;
TsK1(1:jm,im+1) = TsK(1:jm,im);
TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
m_proj('Equidistant Cylindrical','lon',0,'lat',90);
% set the figure background to white
set(gcf,'color','white')
m_pcolor(Tlg,Tlt,TsK1);
shading flat; % gets rid of the grid box edges
colormap(jet);%colormap(flipud(HSV));
caxis([0 1.5])
m_coast('color',[0.3 0.3 0.3]);
m_grid('box','fancy','color',[0 0 0]);
h=colorbar('EastOutside');
title('dz(1)')

subplot(3,2,2)
var_plot=dz(:,:,2);
var_plot(find(fearth==0))=NaN;
TsK = rot90(var_plot);
TsK1 = TsK;
TsK1(1:jm,im+1) = TsK(1:jm,im);
TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
m_proj('Equidistant Cylindrical','lon',0,'lat',90);
% set the figure background to white
set(gcf,'color','white')
m_pcolor(Tlg,Tlt,TsK1);
shading flat; % gets rid of the grid box edges
colormap(jet);%colormap(flipud(HSV));
caxis([0 1.5])
m_coast('color',[0.3 0.3 0.3]);
m_grid('box','fancy','color',[0 0 0]);
h=colorbar('EastOutside');
title('dz(2)')

subplot(3,2,3)
var_plot=dz(:,:,3);
var_plot(find(fearth==0))=NaN;
TsK = rot90(var_plot);
TsK1 = TsK;
TsK1(1:jm,im+1) = TsK(1:jm,im);
TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
m_proj('Equidistant Cylindrical','lon',0,'lat',90);
% set the figure background to white
set(gcf,'color','white')
m_pcolor(Tlg,Tlt,TsK1);
shading flat; % gets rid of the grid box edges
colormap(jet);%colormap(flipud(HSV));
caxis([0 1.5])
m_coast('color',[0.3 0.3 0.3]);
m_grid('box','fancy','color',[0 0 0]);
h=colorbar('EastOutside');
title('dz(3)')

subplot(3,2,4)
var_plot=dz(:,:,4);
var_plot(find(fearth==0))=NaN;
TsK = rot90(var_plot);
TsK1 = TsK;
TsK1(1:jm,im+1) = TsK(1:jm,im);
TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
m_proj('Equidistant Cylindrical','lon',0,'lat',90);
% set the figure background to white
set(gcf,'color','white')
m_pcolor(Tlg,Tlt,TsK1);
shading flat; % gets rid of the grid box edges
colormap(jet);%colormap(flipud(HSV));
caxis([0 1.5])
m_coast('color',[0.3 0.3 0.3]);
m_grid('box','fancy','color',[0 0 0]);
h=colorbar('EastOutside');
title('dz(4)')

subplot(3,2,5)
var_plot=dz(:,:,5);
var_plot(find(fearth==0))=NaN;
TsK = rot90(var_plot);
TsK1 = TsK;
TsK1(1:jm,im+1) = TsK(1:jm,im);
TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
m_proj('Equidistant Cylindrical','lon',0,'lat',90);
% set the figure background to white
set(gcf,'color','white')
m_pcolor(Tlg,Tlt,TsK1);
shading flat; % gets rid of the grid box edges
colormap(jet);%colormap(flipud(HSV));
caxis([0 1.5])
m_coast('color',[0.3 0.3 0.3]);
m_grid('box','fancy','color',[0 0 0]);
h=colorbar('EastOutside');
title('dz(5)')

subplot(3,2,6)
var_plot=dz(:,:,6);
var_plot(find(fearth==0))=NaN;
TsK = rot90(var_plot);
TsK1 = TsK;
TsK1(1:jm,im+1) = TsK(1:jm,im);
TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
m_proj('Equidistant Cylindrical','lon',0,'lat',90);
% set the figure background to white
set(gcf,'color','white')
m_pcolor(Tlg,Tlt,TsK1);
shading flat; % gets rid of the grid box edges
colormap(jet);%colormap(flipud(HSV));
caxis([0 1.5])
m_coast('color',[0.3 0.3 0.3]);
m_grid('box','fancy','color',[0 0 0]);
h=colorbar('EastOutside');
title('dz(6)')
%--------------------------------------------------------------------------



% Print global maps of vegetation fractions
print_vegfrac=1;
if (print_vegfrac == 1)
    var_avg=VEG(:,:,1);
    figure(21)
    TsK = rot90(var_avg);
    TsK1 = TsK;
    TsK1(1:jm,im+1) = TsK(1:jm,im);
    TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
    m_proj('Equidistant Cylindrical','lon',0,'lat',90);
    % set the figure background to white
    set(gcf,'color','white')
    m_pcolor(Tlg,Tlt,TsK1);
    shading flat; % gets rid of the grid box edges
    colormap(jet);%colormap(flipud(HSV));
    m_coast('color',[0.3 0.3 0.3]);
    m_grid('box','fancy','color',[0 0 0]);
    caxis([0 1])
    h=colorbar('EastOutside');
    title('Bright Bare Soil Fraction')

    var_avg=VEG(:,:,10);
    figure(22)
    TsK = rot90(var_avg);
    TsK1 = TsK;
    TsK1(1:jm,im+1) = TsK(1:jm,im);
    TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
    m_proj('Equidistant Cylindrical','lon',0,'lat',90);
    % set the figure background to white
    set(gcf,'color','white')
    m_pcolor(Tlg,Tlt,TsK1);
    shading flat; % gets rid of the grid box edges
    colormap(jet);%colormap(flipud(HSV));
    m_coast('color',[0.3 0.3 0.3]);
    m_grid('box','fancy','color',[0 0 0]);
    caxis([0 1])
    h=colorbar('EastOutside');
    title('Dark Bare Soil Fraction')

    var_avg=VEG(:,:,2);
    figure(23)
    TsK = rot90(var_avg);
    TsK1 = TsK;
    TsK1(1:jm,im+1) = TsK(1:jm,im);
    TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
    m_proj('Equidistant Cylindrical','lon',0,'lat',90);
    % set the figure background to white
    set(gcf,'color','white')
    m_pcolor(Tlg,Tlt,TsK1);
    shading flat; % gets rid of the grid box edges
    colormap(jet);%colormap(flipud(HSV));
    m_coast('color',[0.3 0.3 0.3]);
    m_grid('box','fancy','color',[0 0 0]);
    caxis([0 1])
    h=colorbar('EastOutside');
    title('Tundra Fraction')

    var_avg=VEG(:,:,3);
    figure(24)
    TsK = rot90(var_avg);
    TsK1 = TsK;
    TsK1(1:jm,im+1) = TsK(1:jm,im);
    TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
    m_proj('Equidistant Cylindrical','lon',0,'lat',90);
    % set the figure background to white
    set(gcf,'color','white')
    m_pcolor(Tlg,Tlt,TsK1);
    shading flat; % gets rid of the grid box edges
    colormap(jet);%colormap(flipud(HSV));
    m_coast('color',[0.3 0.3 0.3]);
    m_grid('box','fancy','color',[0 0 0]);
    caxis([0 1])
    h=colorbar('EastOutside');
    title('Grassland Fraction')

    var_avg=VEG(:,:,4);
    figure(25)
    TsK = rot90(var_avg);
    TsK1 = TsK;
    TsK1(1:jm,im+1) = TsK(1:jm,im);
    TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
    m_proj('Equidistant Cylindrical','lon',0,'lat',90);
    % set the figure background to white
    set(gcf,'color','white')
    m_pcolor(Tlg,Tlt,TsK1);
    shading flat; % gets rid of the grid box edges
    colormap(jet);%colormap(flipud(HSV));
    m_coast('color',[0.3 0.3 0.3]);
    m_grid('box','fancy','color',[0 0 0]);
    caxis([0 1])
    h=colorbar('EastOutside');
    title('Shrub/Grassland Fraction')

    var_avg=VEG(:,:,5);
    figure(26)
    TsK = rot90(var_avg);
    TsK1 = TsK;
    TsK1(1:jm,im+1) = TsK(1:jm,im);
    TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
    m_proj('Equidistant Cylindrical','lon',0,'lat',90);
    % set the figure background to white
    set(gcf,'color','white')
    m_pcolor(Tlg,Tlt,TsK1);
    shading flat; % gets rid of the grid box edges
    colormap(jet);%colormap(flipud(HSV));
    m_coast('color',[0.3 0.3 0.3]);
    m_grid('box','fancy','color',[0 0 0]);
    caxis([0 1])
    h=colorbar('EastOutside');
    title('Tree/Grassland Fraction')

    var_avg=VEG(:,:,6);
    figure(27)
    TsK = rot90(var_avg);
    TsK1 = TsK;
    TsK1(1:jm,im+1) = TsK(1:jm,im);
    TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
    m_proj('Equidistant Cylindrical','lon',0,'lat',90);
    % set the figure background to white
    set(gcf,'color','white')
    m_pcolor(Tlg,Tlt,TsK1);
    shading flat; % gets rid of the grid box edges
    colormap(jet);%colormap(flipud(HSV));
    m_coast('color',[0.3 0.3 0.3]);
    m_grid('box','fancy','color',[0 0 0]);
    caxis([0 1])
    h=colorbar('EastOutside');
    title('Deciduous Forest Fraction')

    var_avg=VEG(:,:,7);
    figure(28)
    TsK = rot90(var_avg);
    TsK1 = TsK;
    TsK1(1:jm,im+1) = TsK(1:jm,im);
    TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
    m_proj('Equidistant Cylindrical','lon',0,'lat',90);
    % set the figure background to white
    set(gcf,'color','white')
    m_pcolor(Tlg,Tlt,TsK1);
    shading flat; % gets rid of the grid box edges
    colormap(jet);%colormap(flipud(HSV));
    m_coast('color',[0.3 0.3 0.3]);
    m_grid('box','fancy','color',[0 0 0]);
    caxis([0 1])
    h=colorbar('EastOutside');
    title('Evergreen Forest Fraction')

    var_avg=VEG(:,:,8);
    figure(29)
    TsK = rot90(var_avg);
    TsK1 = TsK;
    TsK1(1:jm,im+1) = TsK(1:jm,im);
    TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
    m_proj('Equidistant Cylindrical','lon',0,'lat',90);
    % set the figure background to white
    set(gcf,'color','white')
    m_pcolor(Tlg,Tlt,TsK1);
    shading flat; % gets rid of the grid box edges
    colormap(jet);%colormap(flipud(HSV));
    m_coast('color',[0.3 0.3 0.3]);
    m_grid('box','fancy','color',[0 0 0]);
    caxis([0 1])
    h=colorbar('EastOutside');
    title('Rainforest Fraction')

    var_avg=VEG(:,:,9);
    figure(30)
    TsK = rot90(var_avg);
    TsK1 = TsK;
    TsK1(1:jm,im+1) = TsK(1:jm,im);
    TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
    m_proj('Equidistant Cylindrical','lon',0,'lat',90);
    % set the figure background to white
    set(gcf,'color','white')
    m_pcolor(Tlg,Tlt,TsK1);
    shading flat; % gets rid of the grid box edges
    colormap(jet);%colormap(flipud(HSV));
    m_coast('color',[0.3 0.3 0.3]);
    m_grid('box','fancy','color',[0 0 0]);
    caxis([0 1])
    h=colorbar('EastOutside');
    title('Cultivation Fraction')
end

% fearth = earth fraction
ngm = 6;%+1   % # of soil layers:% 6 soil(+1 is canopy: layer 0)
LS_NFRAC = 3; % # of land surface fracs:1)bare,2)vegetated,3)lake covered
imt = 5;      % # of soil textures
% shc_soil_texture specific heat capacity of soil texture (J/K/M^3)
% shc_soil_texture(imt)= [2d6,2d6,2d6,2.5d6,2.4d6];
% topmodel input data and standard deviation of the elevation
%--------------------------------------------------------------------------
% fid = fopen('/Users/mpuma/Documents/LSM_2X2_5/lsm_data_dump_f_fort.953'...
%     ,'r','ieee-be.l64');
% if fid == -1
%     disp(message)
% end
% %Initialize matrices
% fearth = zeros(im,jm);
% top_index = zeros(im,jm);
% top_dev = zeros(im,jm);
% dz = zeros(im,jm,ngm); 
% q = zeros(im,jm,imt,ngm);
% qk = zeros(im,jm,imt,ngm);
% sl= zeros(im,jm);
% 
% skip = fread(fid,[4]);
% fearth = fread(fid,[im,jm],'float64','ieee-be.l64');
% top_index = fread(fid,[im,jm],'float64','ieee-be.l64');
% top_dev = fread(fid,[im,jm],'float64','ieee-be.l64');
% 
% for i=1:ngm
%     dz(:,:,i) = fread(fid,[im,jm],'float64','ieee-be.l64');
% end
% for i=1:imt
%     for j=ngm
%         temp = fread(fid,[im,jm],'float64','ieee-be.l64');
%         q(:,:,i,j)=temp;
%     end
% end
% for i=1:imt
%     for j=ngm
%         temp = fread(fid,[im,jm],'float64','ieee-be.l64');
%         qk(:,:,i,j) =temp;
%     end
% end
% sl = fread(fid,[im,jm],'float64','ieee-be.l64');
% fclose(fid);
% 
% q_site  = zeros(imt,ngm);
% q_site(:,:)=q(i_site,j_site,:,:);
