% EVAPOTRANSPIRATION PLOTS (2 deg lat x 2.5 deg long)
% FIGURE A (figures 1-4):
% - annual total ET
% - annual transpiration/ET
% - annual soil evaporation/ET
% - annual canopy evaporation/ET
% 
% FIGURE B (figures 5-7)
% - annual potential evaporation
% - annual precipitation
% - annual dryness index

clear all; close all; clc;
im = 144; jm = 90; 
Tlat = [90:-2:-90]; %  2deg lat cell
Tlon = [-180:2.5:180];%  2.5deg long cell
[Tlg,Tlt]=meshgrid(Tlon,Tlat);
N_hemisph = zeros(im,jm); N_hemisph(1:im,jm/2+1:jm) = 1;
S_hemisph = zeros(im,jm); S_hemisph(1:im,1:jm/2) = 1;

%%%num_yrs_data = ??;
num_yrs = 1;
num_spinup_yr = 0;
yr_start = 1986-num_spinup_yr;
yr_end = yr_start + num_yrs - 1;
num_mnth = 12*num_yrs;
tlsm_mnth = yr_start:1/12:(yr_end+11/12);
tlsm_yr = (yr_start+0.5):(yr_end+0.5);

rhow = 1000; % water density kg/m3
mm_fromkg_per_m2 = 1/rhow*1000; %kg/m2=> m => mm
m_fromkg_per_m2 =  1/rhow;      %kg/m2 => m
dz=[0.09999;0.1725;0.2977;0.5137;0.8864;1.529];
%--------------------------------------------------------------------------
% fearth = earth fraction
ngm = 6;%+1   % # of soil layers:% 6 soil(+1 is canopy: layer 0)
LS_NFRAC = 3; % # of land surface fractions:1)bare, 2)vegetated, 3)covered by lake
imt = 5;      % # of soil textures
% shc_soil_texture specific heat capacity of soil texture (J/K/M^3)
%shc_soil_texture(imt)= [2d6,2d6,2d6,2.5d6,2.4d6];
% topmodel input data and standard deviation of the elevation

%--------------------------------------------------------------------------
fid = fopen('/Users/mpuma/Documents/LSM_2X2_5/lsm_data_dump_f_fort.953'...
    ,'r','ieee-be.l64');
if fid == -1
    disp(message)
end
%Initialize matrices
fearth = zeros(im,jm);
top_index = zeros(im,jm);
top_dev = zeros(im,jm);
dz = zeros(im,jm,ngm); 
q = zeros(im,jm,imt,ngm);
qk = zeros(im,jm,imt,ngm);
sl= zeros(im,jm);

skip = fread(fid,[4]);
fearth = fread(fid,[im,jm],'float64','ieee-be.l64');
top_index = fread(fid,[im,jm],'float64','ieee-be.l64');
top_dev = fread(fid,[im,jm],'float64','ieee-be.l64');
for i=1:ngm
    dz(:,:,i) = fread(fid,[im,jm],'float64','ieee-be.l64');
end
for i=1:imt
    for j=ngm
        q(:,:,i,j) = fread(fid,[im,jm],'float64','ieee-be.l64');
    end
end
for i=1:imt
    for j=ngm
        qk(:,:,i,j) = fread(fid,[im,jm],'float64','ieee-be.l64');
    end
end
sl = fread(fid,[im,jm],'float64','ieee-be.l64');
fclose(fid);

land_mask=zeros(im,jm); land_mask(fearth>0)=1;
i_site=35; j_site=84; %GSWP2 does not have land here
land_mask(i_site,j_site)=0;
%--------------------------------------------------------------------------
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
land_area=grid_area.*fearth;
earth_land=sum(sum(land_area));
%earth_land = 130,580,000 km ice-free land;%148,940,000 km land

%--------------------------------------------------------------------------
% Total evaporation [kg H2O/m2]
% aevap = aevap + (evap_tot(1)*fb +evap_tot(2)*fv)*dts
fid = fopen('/NoBackup/LSM_output/GHY_old/fort.980','r','ieee-be');
if fid == -1
    disp(message)
end

%Initialize matrices
temp = zeros(im,jm);
temp_all = zeros(im,jm,num_mnth);

for i = 1:num_mnth
    skip = fread(fid,[4]);
    title_read = fread(fid,[80],'*char');
    temp = fread(fid,[im,jm],'float32','ieee-be'); %float64 for real*8; float32 for real*4
    temp_all(:,:,i)=temp;
    skip = fread(fid,[4]);
end

E_tot_mean_jan = zeros(im,jm); E_tot_mean_feb = zeros(im,jm);
E_tot_mean_mar = zeros(im,jm); E_tot_mean_apr = zeros(im,jm);
E_tot_mean_may = zeros(im,jm); E_tot_mean_jun = zeros(im,jm);
E_tot_mean_jul = zeros(im,jm); E_tot_mean_aug = zeros(im,jm);
E_tot_mean_sep = zeros(im,jm); E_tot_mean_oct = zeros(im,jm);
E_tot_mean_nov = zeros(im,jm); E_tot_mean_dec = zeros(im,jm);

% Sum individual months over all years of simulation
for i = (num_spinup_yr+1):num_yrs
    E_tot_mean_jan = E_tot_mean_jan+temp_all(:,:,1+(i-1)*12);
    E_tot_mean_feb = E_tot_mean_feb+temp_all(:,:,2+(i-1)*12);
    E_tot_mean_mar = E_tot_mean_mar+temp_all(:,:,3+(i-1)*12);
    E_tot_mean_apr = E_tot_mean_apr+temp_all(:,:,4+(i-1)*12);
    E_tot_mean_may = E_tot_mean_may+temp_all(:,:,5+(i-1)*12);
    E_tot_mean_jun = E_tot_mean_jun+temp_all(:,:,6+(i-1)*12);
    E_tot_mean_jul = E_tot_mean_jul+temp_all(:,:,7+(i-1)*12);
    E_tot_mean_aug = E_tot_mean_aug+temp_all(:,:,8+(i-1)*12);
    E_tot_mean_sep = E_tot_mean_sep+temp_all(:,:,9+(i-1)*12);
    E_tot_mean_oct = E_tot_mean_oct+temp_all(:,:,10+(i-1)*12);
    E_tot_mean_nov = E_tot_mean_nov+temp_all(:,:,11+(i-1)*12);
    E_tot_mean_dec = E_tot_mean_dec+temp_all(:,:,12+(i-1)*12);
end

% Compute average monthly values [kg/m2/month]
num_yrs_avg = num_yrs-num_spinup_yr;
E_tot_mean_jan = E_tot_mean_jan/num_yrs_avg;
E_tot_mean_feb = E_tot_mean_feb/num_yrs_avg;
E_tot_mean_mar = E_tot_mean_mar/num_yrs_avg;
E_tot_mean_apr = E_tot_mean_apr/num_yrs_avg;
E_tot_mean_may = E_tot_mean_may/num_yrs_avg;
E_tot_mean_jun = E_tot_mean_jun/num_yrs_avg;
E_tot_mean_jul = E_tot_mean_jul/num_yrs_avg;
E_tot_mean_aug = E_tot_mean_aug/num_yrs_avg;
E_tot_mean_sep = E_tot_mean_sep/num_yrs_avg;
E_tot_mean_oct = E_tot_mean_oct/num_yrs_avg;
E_tot_mean_nov = E_tot_mean_nov/num_yrs_avg;
E_tot_mean_dec = E_tot_mean_dec/num_yrs_avg;

% Compute annual average values [kg/m2/yr]
E_tot_mean_annual = (E_tot_mean_jan + E_tot_mean_feb+ E_tot_mean_mar + E_tot_mean_apr+ ...
                E_tot_mean_may + E_tot_mean_jun+ E_tot_mean_jul + E_tot_mean_aug+ ...
                E_tot_mean_sep + E_tot_mean_oct+ E_tot_mean_nov + E_tot_mean_dec);

% Compute global total
Etot = E_tot_mean_annual; Etot(find(land_mask<=0))=0; 
Etot=Etot.*(land_area*1000000)/rhow; %kg/m^2/yr * m^2 / (kg/m3)
total_Etot = sum(sum(Etot(find(land_mask>0))))/(1000000*earth_land);%m/yr;
total_Etot=1000*total_Etot; %mm/yr

% Compute time series        
% Monthly average
for i = 1:num_mnth 
   E_tot_time(i) = mean(mean(temp_all(:,:,i).*(temp_all(:,:,i)>=0)));
end
% Yearly average
for i = 1:num_yrs 
   E_tot_tyr(i) = mean(E_tot_time(12*(i-1)+1:i*12));   
end

%--------------------------------------------------------------------------
% Canopy evaporation [kg H2O/m2]
% aevapw = aevapw + evapvw*fw*(1-fr_snow(2)*fm)*fv*dts
fid = fopen('/NoBackup/LSM_output/GHY_old/fort.981','r','ieee-be');
if fid == -1
    disp(message)
end

%Initialize matrices
temp = zeros(im,jm);
temp_all = zeros(im,jm,num_mnth);

for i = 1:num_mnth
    skip = fread(fid,[4]);
    title_read = fread(fid,[80],'*char');
    temp = fread(fid,[im,jm],'float32','ieee-be'); %float64 for real*8; float32 for real*4
    temp_all(:,:,i)=temp;
    skip = fread(fid,[4]); 
end

E_can_mean_jan = zeros(im,jm); E_can_mean_feb = zeros(im,jm);
E_can_mean_mar = zeros(im,jm); E_can_mean_apr = zeros(im,jm);
E_can_mean_may = zeros(im,jm); E_can_mean_jun = zeros(im,jm);
E_can_mean_jul = zeros(im,jm); E_can_mean_aug = zeros(im,jm);
E_can_mean_sep = zeros(im,jm); E_can_mean_oct = zeros(im,jm);
E_can_mean_nov = zeros(im,jm); E_can_mean_dec = zeros(im,jm);

% Sum individual months over all years of simulation
for i = (num_spinup_yr+1):num_yrs
    E_can_mean_jan = E_can_mean_jan+temp_all(:,:,1+(i-1)*12);
    E_can_mean_feb = E_can_mean_feb+temp_all(:,:,2+(i-1)*12);
    E_can_mean_mar = E_can_mean_mar+temp_all(:,:,3+(i-1)*12);
    E_can_mean_apr = E_can_mean_apr+temp_all(:,:,4+(i-1)*12);
    E_can_mean_may = E_can_mean_may+temp_all(:,:,5+(i-1)*12);
    E_can_mean_jun = E_can_mean_jun+temp_all(:,:,6+(i-1)*12);
    E_can_mean_jul = E_can_mean_jul+temp_all(:,:,7+(i-1)*12);
    E_can_mean_aug = E_can_mean_aug+temp_all(:,:,8+(i-1)*12);
    E_can_mean_sep = E_can_mean_sep+temp_all(:,:,9+(i-1)*12);
    E_can_mean_oct = E_can_mean_oct+temp_all(:,:,10+(i-1)*12);
    E_can_mean_nov = E_can_mean_nov+temp_all(:,:,11+(i-1)*12);
    E_can_mean_dec = E_can_mean_dec+temp_all(:,:,12+(i-1)*12);
end

% Compute average monthly values [kg/m2/month]
num_yrs_avg = num_yrs-num_spinup_yr;
E_can_mean_jan = E_can_mean_jan/num_yrs_avg;
E_can_mean_feb = E_can_mean_feb/num_yrs_avg;
E_can_mean_mar = E_can_mean_mar/num_yrs_avg;
E_can_mean_apr = E_can_mean_apr/num_yrs_avg;
E_can_mean_may = E_can_mean_may/num_yrs_avg;
E_can_mean_jun = E_can_mean_jun/num_yrs_avg;
E_can_mean_jul = E_can_mean_jul/num_yrs_avg;
E_can_mean_aug = E_can_mean_aug/num_yrs_avg;
E_can_mean_sep = E_can_mean_sep/num_yrs_avg;
E_can_mean_oct = E_can_mean_oct/num_yrs_avg;
E_can_mean_nov = E_can_mean_nov/num_yrs_avg;
E_can_mean_dec = E_can_mean_dec/num_yrs_avg;

% Compute annual average values [kg/m2/yr]
E_can_mean_annual = (E_can_mean_jan + E_can_mean_feb+ E_can_mean_mar + E_can_mean_apr+ ...
                E_can_mean_may + E_can_mean_jun+ E_can_mean_jul + E_can_mean_aug+ ...
                E_can_mean_sep + E_can_mean_oct+ E_can_mean_nov + E_can_mean_dec);

% Compute global total
Ecan = E_can_mean_annual; Ecan(find(land_mask<=0))=0; 
Ecan=Ecan.*(land_area*1000000)/rhow; %kg/m^2/yr * m^2 / (kg/m3)
total_Ecan = sum(sum(Ecan(find(land_mask>0))))/(1000000*earth_land);%m/yr;
total_Ecan=1000*total_Ecan; %mm/yr


% Compute time series        
% Monthly average
for i = 1:num_mnth 
   E_can_time(i) = mean(mean(temp_all(:,:,i).*(temp_all(:,:,i)>=0)));
end
% Yearly average
for i = 1:num_yrs 
   E_can_tyr(i) = mean(E_can_time(12*(i-1)+1:i*12));   
end


%--------------------------------------------------------------------------
% Transpiration [kg H2O/m2]
% aevapd = aevapd + evapvd*fd*(1-fr_snow(2)*fm)*fv*dts
fid = fopen('/NoBackup/LSM_output/GHY_old/fort.982','r','ieee-be');
if fid == -1
    disp(message)
end

%Initialize matrices
temp = zeros(im,jm);
temp_all = zeros(im,jm,num_mnth);

for i = 1:num_mnth
    skip = fread(fid,[4]);
    title_read = fread(fid,[80],'*char');
    temp = fread(fid,[im,jm],'float32','ieee-be'); %float64 for real*8; float32 for real*4
    temp_all(:,:,i)=temp;
    skip = fread(fid,[4]);   
end

E_tran_mean_jan = zeros(im,jm); E_tran_mean_feb = zeros(im,jm);
E_tran_mean_mar = zeros(im,jm); E_tran_mean_apr = zeros(im,jm);
E_tran_mean_may = zeros(im,jm); E_tran_mean_jun = zeros(im,jm);
E_tran_mean_jul = zeros(im,jm); E_tran_mean_aug = zeros(im,jm);
E_tran_mean_sep = zeros(im,jm); E_tran_mean_oct = zeros(im,jm);
E_tran_mean_nov = zeros(im,jm); E_tran_mean_dec = zeros(im,jm);

% Sum individual months over all years of simulation
for i = (num_spinup_yr+1):num_yrs
    E_tran_mean_jan = E_tran_mean_jan+temp_all(:,:,1+(i-1)*12);
    E_tran_mean_feb = E_tran_mean_feb+temp_all(:,:,2+(i-1)*12);
    E_tran_mean_mar = E_tran_mean_mar+temp_all(:,:,3+(i-1)*12);
    E_tran_mean_apr = E_tran_mean_apr+temp_all(:,:,4+(i-1)*12);
    E_tran_mean_may = E_tran_mean_may+temp_all(:,:,5+(i-1)*12);
    E_tran_mean_jun = E_tran_mean_jun+temp_all(:,:,6+(i-1)*12);
    E_tran_mean_jul = E_tran_mean_jul+temp_all(:,:,7+(i-1)*12);
    E_tran_mean_aug = E_tran_mean_aug+temp_all(:,:,8+(i-1)*12);
    E_tran_mean_sep = E_tran_mean_sep+temp_all(:,:,9+(i-1)*12);
    E_tran_mean_oct = E_tran_mean_oct+temp_all(:,:,10+(i-1)*12);
    E_tran_mean_nov = E_tran_mean_nov+temp_all(:,:,11+(i-1)*12);
    E_tran_mean_dec = E_tran_mean_dec+temp_all(:,:,12+(i-1)*12);  
end

% Compute average monthly values [kg/m2/month]
num_yrs_avg = num_yrs-num_spinup_yr;
E_tran_mean_jan = E_tran_mean_jan/num_yrs_avg;
E_tran_mean_feb = E_tran_mean_feb/num_yrs_avg;
E_tran_mean_mar = E_tran_mean_mar/num_yrs_avg;
E_tran_mean_apr = E_tran_mean_apr/num_yrs_avg;
E_tran_mean_may = E_tran_mean_may/num_yrs_avg;
E_tran_mean_jun = E_tran_mean_jun/num_yrs_avg;
E_tran_mean_jul = E_tran_mean_jul/num_yrs_avg;
E_tran_mean_aug = E_tran_mean_aug/num_yrs_avg;
E_tran_mean_sep = E_tran_mean_sep/num_yrs_avg;
E_tran_mean_oct = E_tran_mean_oct/num_yrs_avg;
E_tran_mean_nov = E_tran_mean_nov/num_yrs_avg;
E_tran_mean_dec = E_tran_mean_dec/num_yrs_avg;

% Compute annual average values [kg/m2/yr]
E_tran_mean_annual = (E_tran_mean_jan + E_tran_mean_feb+ E_tran_mean_mar + E_tran_mean_apr+ ...
                E_tran_mean_may + E_tran_mean_jun+ E_tran_mean_jul + E_tran_mean_aug+ ...
                E_tran_mean_sep + E_tran_mean_oct+ E_tran_mean_nov + E_tran_mean_dec);

% Compute global value
Etran = E_tran_mean_annual; Etran(find(land_mask<=0))=0; 
Etran=Etran.*(land_area*1000000)/rhow; %kg/m^2/yr * m^2 / (kg/m3)
total_Etran = sum(sum(Etran(find(land_mask>0))))/(1000000*earth_land);%m/yr;
total_Etran=1000*total_Etran; %mm/yr

% Compute time series        
% Monthly average
for i = 1:num_mnth 
   E_tran_time(i) = mean(mean(temp_all(:,:,i).*(temp_all(:,:,i)>=0)));
end
% Yearly average
for i = 1:num_yrs 
   E_tran_tyr(i) = mean(E_tran_time(12*(i-1)+1:i*12));   
end

%--------------------------------------------------------------------------
% Bare soil evaporation [kg H2O/m2]
% aevapb = aevapb + evapb*(1.d0-fr_snow(1))*fb*dts
fid = fopen('/NoBackup/LSM_output/GHY_old/fort.983','r','ieee-be');
if fid == -1
    disp(message)
end

%Initialize matrices
temp = zeros(im,jm);
temp_all = zeros(im,jm,num_mnth);

for i = 1:num_mnth
    skip = fread(fid,[4]);
    title_read = fread(fid,[80],'*char');
    temp = fread(fid,[im,jm],'float32','ieee-be'); %float64 for real*8; float32 for real*4
    temp_all(:,:,i)=temp;
    skip = fread(fid,[4]);
end

E_bare_mean_jan = zeros(im,jm); E_bare_mean_feb = zeros(im,jm);
E_bare_mean_mar = zeros(im,jm); E_bare_mean_apr = zeros(im,jm);
E_bare_mean_may = zeros(im,jm); E_bare_mean_jun = zeros(im,jm);
E_bare_mean_jul = zeros(im,jm); E_bare_mean_aug = zeros(im,jm);
E_bare_mean_sep = zeros(im,jm); E_bare_mean_oct = zeros(im,jm);
E_bare_mean_nov = zeros(im,jm); E_bare_mean_dec = zeros(im,jm);

% Sum individual months over all years of simulation
for i = (num_spinup_yr+1):num_yrs
    E_bare_mean_jan = E_bare_mean_jan+temp_all(:,:,1+(i-1)*12);
    E_bare_mean_feb = E_bare_mean_feb+temp_all(:,:,2+(i-1)*12);
    E_bare_mean_mar = E_bare_mean_mar+temp_all(:,:,3+(i-1)*12);
    E_bare_mean_apr = E_bare_mean_apr+temp_all(:,:,4+(i-1)*12);
    E_bare_mean_may = E_bare_mean_may+temp_all(:,:,5+(i-1)*12);
    E_bare_mean_jun = E_bare_mean_jun+temp_all(:,:,6+(i-1)*12);
    E_bare_mean_jul = E_bare_mean_jul+temp_all(:,:,7+(i-1)*12);
    E_bare_mean_aug = E_bare_mean_aug+temp_all(:,:,8+(i-1)*12);
    E_bare_mean_sep = E_bare_mean_sep+temp_all(:,:,9+(i-1)*12);
    E_bare_mean_oct = E_bare_mean_oct+temp_all(:,:,10+(i-1)*12);
    E_bare_mean_nov = E_bare_mean_nov+temp_all(:,:,11+(i-1)*12);
    E_bare_mean_dec = E_bare_mean_dec+temp_all(:,:,12+(i-1)*12);
end

% Compute average monthly values [kg/m2/month]
num_yrs_avg = num_yrs-num_spinup_yr;
E_bare_mean_jan = E_bare_mean_jan/num_yrs_avg;
E_bare_mean_feb = E_bare_mean_feb/num_yrs_avg;
E_bare_mean_mar = E_bare_mean_mar/num_yrs_avg;
E_bare_mean_apr = E_bare_mean_apr/num_yrs_avg;
E_bare_mean_may = E_bare_mean_may/num_yrs_avg;
E_bare_mean_jun = E_bare_mean_jun/num_yrs_avg;
E_bare_mean_jul = E_bare_mean_jul/num_yrs_avg;
E_bare_mean_aug = E_bare_mean_aug/num_yrs_avg;
E_bare_mean_sep = E_bare_mean_sep/num_yrs_avg;
E_bare_mean_oct = E_bare_mean_oct/num_yrs_avg;
E_bare_mean_nov = E_bare_mean_nov/num_yrs_avg;
E_bare_mean_dec = E_bare_mean_dec/num_yrs_avg;

% Compute annual average values [kg/m2/yr]
E_bare_mean_annual = (E_bare_mean_jan + E_bare_mean_feb+ E_bare_mean_mar + E_bare_mean_apr+ ...
                E_bare_mean_may + E_bare_mean_jun+ E_bare_mean_jul + E_bare_mean_aug+ ...
                E_bare_mean_sep + E_bare_mean_oct+ E_bare_mean_nov + E_bare_mean_dec);

% Compute global total
Ebare = E_bare_mean_annual; Ebare(find(land_mask<=0))=0; 
Ebare=Ebare.*(land_area*1000000)/rhow; %kg/m^2/yr * m^2 / (kg/m3)
total_Ebare = sum(sum(Ebare(find(land_mask>0))))/(1000000*earth_land);%m/yr;
total_Ebare=1000*total_Ebare; %mm/yr

% Compute time series        
% Monthly average
for i = 1:num_mnth 
   E_bare_time(i) = mean(mean(temp_all(:,:,i).*(temp_all(:,:,i)>=0)));
end
% Yearly average
for i = 1:num_yrs 
   E_bare_tyr(i) = mean(E_bare_time(12*(i-1)+1:i*12));   
end

%--------------------------------------------------------------------------
% Potential evaporation [kg H2O/m2]
fid = fopen('/NoBackup/LSM_output/GHY_old/fort.998','r','ieee-be');
if fid == -1
    disp(message)
end

%Initialize matrices
temp = zeros(im,jm);
temp_all = zeros(im,jm,num_mnth);

for i = 1:num_mnth
    skip = fread(fid,[4]);
    title_read = fread(fid,[80],'*char');
    temp = fread(fid,[im,jm],'float32','ieee-be'); %float64 for real*8; float32 for real*4
    temp_all(:,:,i)=temp;
    skip = fread(fid,[4]);
end

E_pot_mean_jan = zeros(im,jm); E_pot_mean_feb = zeros(im,jm);
E_pot_mean_mar = zeros(im,jm); E_pot_mean_apr = zeros(im,jm);
E_pot_mean_may = zeros(im,jm); E_pot_mean_jun = zeros(im,jm);
E_pot_mean_jul = zeros(im,jm); E_pot_mean_aug = zeros(im,jm);
E_pot_mean_sep = zeros(im,jm); E_pot_mean_oct = zeros(im,jm);
E_pot_mean_nov = zeros(im,jm); E_pot_mean_dec = zeros(im,jm);

% Sum individual months over all years of simulation
for i = (num_spinup_yr+1):num_yrs
    E_pot_mean_jan = E_pot_mean_jan+temp_all(:,:,1+(i-1)*12);
    E_pot_mean_feb = E_pot_mean_feb+temp_all(:,:,2+(i-1)*12);
    E_pot_mean_mar = E_pot_mean_mar+temp_all(:,:,3+(i-1)*12);
    E_pot_mean_apr = E_pot_mean_apr+temp_all(:,:,4+(i-1)*12);
    E_pot_mean_may = E_pot_mean_may+temp_all(:,:,5+(i-1)*12);
    E_pot_mean_jun = E_pot_mean_jun+temp_all(:,:,6+(i-1)*12);
    E_pot_mean_jul = E_pot_mean_jul+temp_all(:,:,7+(i-1)*12);
    E_pot_mean_aug = E_pot_mean_aug+temp_all(:,:,8+(i-1)*12);
    E_pot_mean_sep = E_pot_mean_sep+temp_all(:,:,9+(i-1)*12);
    E_pot_mean_oct = E_pot_mean_oct+temp_all(:,:,10+(i-1)*12);
    E_pot_mean_nov = E_pot_mean_nov+temp_all(:,:,11+(i-1)*12);
    E_pot_mean_dec = E_pot_mean_dec+temp_all(:,:,12+(i-1)*12); 
end

% Compute average monthly values [kg/m2/month]
num_yrs_avg = num_yrs-num_spinup_yr;
E_pot_mean_jan = E_pot_mean_jan/num_yrs_avg;
E_pot_mean_feb = E_pot_mean_feb/num_yrs_avg;
E_pot_mean_mar = E_pot_mean_mar/num_yrs_avg;
E_pot_mean_apr = E_pot_mean_apr/num_yrs_avg;
E_pot_mean_may = E_pot_mean_may/num_yrs_avg;
E_pot_mean_jun = E_pot_mean_jun/num_yrs_avg;
E_pot_mean_jul = E_pot_mean_jul/num_yrs_avg;
E_pot_mean_aug = E_pot_mean_aug/num_yrs_avg;
E_pot_mean_sep = E_pot_mean_sep/num_yrs_avg;
E_pot_mean_oct = E_pot_mean_oct/num_yrs_avg;
E_pot_mean_nov = E_pot_mean_nov/num_yrs_avg;
E_pot_mean_dec = E_pot_mean_dec/num_yrs_avg;

% Assign monthly values to single matrix [kg/m2/month]
monthly_Epot_avg(:,:,1)=E_pot_mean_jan;
monthly_Epot_avg(:,:,2)=E_pot_mean_feb;
monthly_Epot_avg(:,:,3)=E_pot_mean_mar;
monthly_Epot_avg(:,:,4)=E_pot_mean_apr;
monthly_Epot_avg(:,:,5)=E_pot_mean_may;
monthly_Epot_avg(:,:,6)=E_pot_mean_jun;
monthly_Epot_avg(:,:,7)=E_pot_mean_jul;
monthly_Epot_avg(:,:,8)=E_pot_mean_aug;
monthly_Epot_avg(:,:,9)=E_pot_mean_sep;
monthly_Epot_avg(:,:,10)=E_pot_mean_oct;
monthly_Epot_avg(:,:,11)=E_pot_mean_nov;
monthly_Epot_avg(:,:,12)=E_pot_mean_dec;

% Compute annual average values [kg/m2/yr]
E_pot_mean_annual = (E_pot_mean_jan + E_pot_mean_feb+ E_pot_mean_mar + E_pot_mean_apr+ ...
                E_pot_mean_may + E_pot_mean_jun+ E_pot_mean_jul + E_pot_mean_aug+ ...
                E_pot_mean_sep + E_pot_mean_oct+ E_pot_mean_nov + E_pot_mean_dec);

% Compute global total
Epot = E_pot_mean_annual; Epot(find(land_mask<=0))=0; 
Epot=Epot.*(land_area*1000000)/rhow; %kg/m^2/yr * m^2 / (kg/m3)
total_Epot = sum(sum(Epot(find(land_mask>0))))/(1000000*earth_land);%m/yr;
total_Epot=1000*total_Epot; %mm/yr

% Compute time series        
% Monthly average
for i = 1:num_mnth 
   E_pot_time(i) = mean(mean(temp_all(:,:,i).*(temp_all(:,:,i)>=0)));
end
% Yearly average
for i = 1:num_yrs 
   E_pot_tyr(i) = mean(E_pot_time(12*(i-1)+1:i*12));   
end

%--------------------------------------------------------------------------
% Bare soil moisture layer 1 [m]
fid = fopen('/NoBackup/LSM_output/GHY_old/fort.1003','r','ieee-be');
if fid == -1
    disp(message)
end

layer_thickness = 0.09999;
%Initialize matrices
temp = zeros(im,jm);
temp_all = zeros(im,jm,num_mnth);

for i = 1:num_mnth
    skip = fread(fid,[4]);
    title_read = fread(fid,[80],'*char');
    temp = fread(fid,[im,jm],'float32','ieee-be'); %float64 for real*8; float32 for real*4
    temp_all(:,:,i)=temp/layer_thickness;
    skip = fread(fid,[4]);
end

SW_bare_mean_jan = zeros(im,jm); SW_bare_mean_feb = zeros(im,jm);
SW_bare_mean_mar = zeros(im,jm); SW_bare_mean_apr = zeros(im,jm);
SW_bare_mean_may = zeros(im,jm); SW_bare_mean_jun = zeros(im,jm);
SW_bare_mean_jul = zeros(im,jm); SW_bare_mean_aug = zeros(im,jm);
SW_bare_mean_sep = zeros(im,jm); SW_bare_mean_oct = zeros(im,jm);
SW_bare_mean_nov = zeros(im,jm); SW_bare_mean_dec = zeros(im,jm);

% Sum individual months over all years of simulation
for i = (num_spinup_yr+1):num_yrs
    SW_bare_mean_jan = SW_bare_mean_jan+temp_all(:,:,1+(i-1)*12);
    SW_bare_mean_feb = SW_bare_mean_feb+temp_all(:,:,2+(i-1)*12);
    SW_bare_mean_mar = SW_bare_mean_mar+temp_all(:,:,3+(i-1)*12);
    SW_bare_mean_apr = SW_bare_mean_apr+temp_all(:,:,4+(i-1)*12);
    SW_bare_mean_may = SW_bare_mean_may+temp_all(:,:,5+(i-1)*12);
    SW_bare_mean_jun = SW_bare_mean_jun+temp_all(:,:,6+(i-1)*12);
    SW_bare_mean_jul = SW_bare_mean_jul+temp_all(:,:,7+(i-1)*12);
    SW_bare_mean_aug = SW_bare_mean_aug+temp_all(:,:,8+(i-1)*12);
    SW_bare_mean_sep = SW_bare_mean_sep+temp_all(:,:,9+(i-1)*12);
    SW_bare_mean_oct = SW_bare_mean_oct+temp_all(:,:,10+(i-1)*12);
    SW_bare_mean_nov = SW_bare_mean_nov+temp_all(:,:,11+(i-1)*12);
    SW_bare_mean_dec = SW_bare_mean_dec+temp_all(:,:,12+(i-1)*12);
end

% Compute average monthly values [-]
num_yrs_avg = num_yrs-num_spinup_yr;
SW_bare_mean_jan = SW_bare_mean_jan/num_yrs_avg;
SW_bare_mean_feb = SW_bare_mean_feb/num_yrs_avg;
SW_bare_mean_mar = SW_bare_mean_mar/num_yrs_avg;
SW_bare_mean_apr = SW_bare_mean_apr/num_yrs_avg;
SW_bare_mean_may = SW_bare_mean_may/num_yrs_avg;
SW_bare_mean_jun = SW_bare_mean_jun/num_yrs_avg;
SW_bare_mean_jul = SW_bare_mean_jul/num_yrs_avg;
SW_bare_mean_aug = SW_bare_mean_aug/num_yrs_avg;
SW_bare_mean_sep = SW_bare_mean_sep/num_yrs_avg;
SW_bare_mean_oct = SW_bare_mean_oct/num_yrs_avg;
SW_bare_mean_nov = SW_bare_mean_nov/num_yrs_avg;
SW_bare_mean_dec = SW_bare_mean_dec/num_yrs_avg;

% Compute annual average values [-]
SW_bare_mean_annual = (SW_bare_mean_jan + SW_bare_mean_feb+ SW_bare_mean_mar + SW_bare_mean_apr+ ...
                SW_bare_mean_may + SW_bare_mean_jun+ SW_bare_mean_jul + SW_bare_mean_aug+ ...
                SW_bare_mean_sep + SW_bare_mean_oct+ SW_bare_mean_nov + SW_bare_mean_dec)/12;

% Compute global total
SW_bare = SW_bare_mean_annual; SW_bare(find(land_mask<=0))=0; 
SW_bare=SW_bare.*(land_area*1000000); %[-] * m^2
total_SW_bare = sum(sum(SW_bare(find(land_mask>0))))/(1000000*earth_land);%m^2/m^2
total_SW_bare=total_SW_bare; %[-]

% Compute time series        
% Monthly average
for i = 1:num_mnth 
   SW_bare_time(i) = mean(mean(temp_all(:,:,i).*(temp_all(:,:,i)>=0)));
end
% Yearly average
for i = 1:num_yrs 
   SW_bare_tyr(i) = mean(SW_bare_time(12*(i-1)+1:i*12));   
end

%--------------------------------------------------------------------------
% Bare soil moisture layer 2 [m]
fid = fopen('/NoBackup/LSM_output/GHY_old/fort.1004','r','ieee-be');
if fid == -1
    disp(message)
end

layer_thickness = 0.1725;
%Initialize matrices
temp = zeros(im,jm);
temp_all = zeros(im,jm,num_mnth);

for i = 1:num_mnth
    skip = fread(fid,[4]);
    title_read = fread(fid,[80],'*char');
    temp = fread(fid,[im,jm],'float32','ieee-be'); %float64 for real*8; float32 for real*4
    temp_all(:,:,i)=temp/layer_thickness;
    skip = fread(fid,[4]);
end

SW2_bare_mean_jan = zeros(im,jm); SW2_bare_mean_feb = zeros(im,jm);
SW2_bare_mean_mar = zeros(im,jm); SW2_bare_mean_apr = zeros(im,jm);
SW2_bare_mean_may = zeros(im,jm); SW2_bare_mean_jun = zeros(im,jm);
SW2_bare_mean_jul = zeros(im,jm); SW2_bare_mean_aug = zeros(im,jm);
SW2_bare_mean_sep = zeros(im,jm); SW2_bare_mean_oct = zeros(im,jm);
SW2_bare_mean_nov = zeros(im,jm); SW2_bare_mean_dec = zeros(im,jm);

% Sum individual months over all years of simulation
for i = (num_spinup_yr+1):num_yrs
    SW2_bare_mean_jan = SW2_bare_mean_jan+temp_all(:,:,1+(i-1)*12);
    SW2_bare_mean_feb = SW2_bare_mean_feb+temp_all(:,:,2+(i-1)*12);
    SW2_bare_mean_mar = SW2_bare_mean_mar+temp_all(:,:,3+(i-1)*12);
    SW2_bare_mean_apr = SW2_bare_mean_apr+temp_all(:,:,4+(i-1)*12);
    SW2_bare_mean_may = SW2_bare_mean_may+temp_all(:,:,5+(i-1)*12);
    SW2_bare_mean_jun = SW2_bare_mean_jun+temp_all(:,:,6+(i-1)*12);
    SW2_bare_mean_jul = SW2_bare_mean_jul+temp_all(:,:,7+(i-1)*12);
    SW2_bare_mean_aug = SW2_bare_mean_aug+temp_all(:,:,8+(i-1)*12);
    SW2_bare_mean_sep = SW2_bare_mean_sep+temp_all(:,:,9+(i-1)*12);
    SW2_bare_mean_oct = SW2_bare_mean_oct+temp_all(:,:,10+(i-1)*12);
    SW2_bare_mean_nov = SW2_bare_mean_nov+temp_all(:,:,11+(i-1)*12);
    SW2_bare_mean_dec = SW2_bare_mean_dec+temp_all(:,:,12+(i-1)*12);
end

% Compute average monthly values [-]
num_yrs_avg = num_yrs-num_spinup_yr;
SW2_bare_mean_jan = SW2_bare_mean_jan/num_yrs_avg;
SW2_bare_mean_feb = SW2_bare_mean_feb/num_yrs_avg;
SW2_bare_mean_mar = SW2_bare_mean_mar/num_yrs_avg;
SW2_bare_mean_apr = SW2_bare_mean_apr/num_yrs_avg;
SW2_bare_mean_may = SW2_bare_mean_may/num_yrs_avg;
SW2_bare_mean_jun = SW2_bare_mean_jun/num_yrs_avg;
SW2_bare_mean_jul = SW2_bare_mean_jul/num_yrs_avg;
SW2_bare_mean_aug = SW2_bare_mean_aug/num_yrs_avg;
SW2_bare_mean_sep = SW2_bare_mean_sep/num_yrs_avg;
SW2_bare_mean_oct = SW2_bare_mean_oct/num_yrs_avg;
SW2_bare_mean_nov = SW2_bare_mean_nov/num_yrs_avg;
SW2_bare_mean_dec = SW2_bare_mean_dec/num_yrs_avg;

% Compute annual average values [-]
SW2_bare_mean_annual = (SW2_bare_mean_jan + SW2_bare_mean_feb+ SW2_bare_mean_mar + SW2_bare_mean_apr+ ...
                SW2_bare_mean_may + SW2_bare_mean_jun+ SW2_bare_mean_jul + SW2_bare_mean_aug+ ...
                SW2_bare_mean_sep + SW2_bare_mean_oct+ SW2_bare_mean_nov + SW2_bare_mean_dec)/12;

% Compute global total
SW2_bare = SW2_bare_mean_annual; SW2_bare(find(land_mask<=0))=0; 
SW2_bare=SW2_bare.*(land_area*1000000); %[-] * m^2
total_SW2_bare = sum(sum(SW2_bare(find(land_mask>0))))/(1000000*earth_land);%m^2/m^2
total_SW2_bare=total_SW2_bare; %[-]

% Compute time series        
% Monthly average
for i = 1:num_mnth 
   SW2_bare_time(i) = mean(mean(temp_all(:,:,i).*(temp_all(:,:,i)>=0)));
end
% Yearly average
for i = 1:num_yrs 
   SW2_bare_tyr(i) = mean(SW2_bare_time(12*(i-1)+1:i*12));   
end
%--------------------------------------------------------------------------
% Bare soil moisture layer 3 [m]
fid = fopen('/NoBackup/LSM_output/GHY_old/fort.1005','r','ieee-be');
if fid == -1
    disp(message)
end

layer_thickness = 0.2977;
%Initialize matrices
temp = zeros(im,jm);
temp_all = zeros(im,jm,num_mnth);

for i = 1:num_mnth
    skip = fread(fid,[4]);
    title_read = fread(fid,[80],'*char');
    temp = fread(fid,[im,jm],'float32','ieee-be'); %float64 for real*8; float32 for real*4
    temp_all(:,:,i)=temp/layer_thickness;
    skip = fread(fid,[4]);
end

SW3_bare_mean_jan = zeros(im,jm); SW3_bare_mean_feb = zeros(im,jm);
SW3_bare_mean_mar = zeros(im,jm); SW3_bare_mean_apr = zeros(im,jm);
SW3_bare_mean_may = zeros(im,jm); SW3_bare_mean_jun = zeros(im,jm);
SW3_bare_mean_jul = zeros(im,jm); SW3_bare_mean_aug = zeros(im,jm);
SW3_bare_mean_sep = zeros(im,jm); SW3_bare_mean_oct = zeros(im,jm);
SW3_bare_mean_nov = zeros(im,jm); SW3_bare_mean_dec = zeros(im,jm);

% Sum individual months over all years of simulation
for i = (num_spinup_yr+1):num_yrs
    SW3_bare_mean_jan = SW3_bare_mean_jan+temp_all(:,:,1+(i-1)*12);
    SW3_bare_mean_feb = SW3_bare_mean_feb+temp_all(:,:,2+(i-1)*12);
    SW3_bare_mean_mar = SW3_bare_mean_mar+temp_all(:,:,3+(i-1)*12);
    SW3_bare_mean_apr = SW3_bare_mean_apr+temp_all(:,:,4+(i-1)*12);
    SW3_bare_mean_may = SW3_bare_mean_may+temp_all(:,:,5+(i-1)*12);
    SW3_bare_mean_jun = SW3_bare_mean_jun+temp_all(:,:,6+(i-1)*12);
    SW3_bare_mean_jul = SW3_bare_mean_jul+temp_all(:,:,7+(i-1)*12);
    SW3_bare_mean_aug = SW3_bare_mean_aug+temp_all(:,:,8+(i-1)*12);
    SW3_bare_mean_sep = SW3_bare_mean_sep+temp_all(:,:,9+(i-1)*12);
    SW3_bare_mean_oct = SW3_bare_mean_oct+temp_all(:,:,10+(i-1)*12);
    SW3_bare_mean_nov = SW3_bare_mean_nov+temp_all(:,:,11+(i-1)*12);
    SW3_bare_mean_dec = SW3_bare_mean_dec+temp_all(:,:,12+(i-1)*12);
end

% Compute average monthly values [-]
num_yrs_avg = num_yrs-num_spinup_yr;
SW3_bare_mean_jan = SW3_bare_mean_jan/num_yrs_avg;
SW3_bare_mean_feb = SW3_bare_mean_feb/num_yrs_avg;
SW3_bare_mean_mar = SW3_bare_mean_mar/num_yrs_avg;
SW3_bare_mean_apr = SW3_bare_mean_apr/num_yrs_avg;
SW3_bare_mean_may = SW3_bare_mean_may/num_yrs_avg;
SW3_bare_mean_jun = SW3_bare_mean_jun/num_yrs_avg;
SW3_bare_mean_jul = SW3_bare_mean_jul/num_yrs_avg;
SW3_bare_mean_aug = SW3_bare_mean_aug/num_yrs_avg;
SW3_bare_mean_sep = SW3_bare_mean_sep/num_yrs_avg;
SW3_bare_mean_oct = SW3_bare_mean_oct/num_yrs_avg;
SW3_bare_mean_nov = SW3_bare_mean_nov/num_yrs_avg;
SW3_bare_mean_dec = SW3_bare_mean_dec/num_yrs_avg;

% Compute annual average values [-]
SW3_bare_mean_annual = (SW3_bare_mean_jan + SW3_bare_mean_feb+ SW3_bare_mean_mar + SW3_bare_mean_apr+ ...
                SW3_bare_mean_may + SW3_bare_mean_jun+ SW3_bare_mean_jul + SW3_bare_mean_aug+ ...
                SW3_bare_mean_sep + SW3_bare_mean_oct+ SW3_bare_mean_nov + SW3_bare_mean_dec)/12;

% Compute global total
SW3_bare = SW3_bare_mean_annual; SW3_bare(find(land_mask<=0))=0; 
SW3_bare=SW3_bare.*(land_area*1000000); %[-] * m^2
total_SW3_bare = sum(sum(SW3_bare(find(land_mask>0))))/(1000000*earth_land);%m^2/m^2
total_SW3_bare=total_SW3_bare; %[-]

% Compute time series        
% Monthly average
for i = 1:num_mnth 
   SW3_bare_time(i) = mean(mean(temp_all(:,:,i).*(temp_all(:,:,i)>=0)));
end
% Yearly average
for i = 1:num_yrs 
   SW3_bare_tyr(i) = mean(SW3_bare_time(12*(i-1)+1:i*12));   
end
%--------------------------------------------------------------------------
% Bare soil moisture layer 4 [m]
fid = fopen('/NoBackup/LSM_output/GHY_old/fort.1006','r','ieee-be');
if fid == -1
    disp(message)
end

layer_thickness = 0.5137;
%Initialize matrices
temp = zeros(im,jm);
temp_all = zeros(im,jm,num_mnth);

for i = 1:num_mnth
    skip = fread(fid,[4]);
    title_read = fread(fid,[80],'*char');
    temp = fread(fid,[im,jm],'float32','ieee-be'); %float64 for real*8; float32 for real*4
    temp_all(:,:,i)=temp/layer_thickness;
    skip = fread(fid,[4]);
end

SW4_bare_mean_jan = zeros(im,jm); SW4_bare_mean_feb = zeros(im,jm);
SW4_bare_mean_mar = zeros(im,jm); SW4_bare_mean_apr = zeros(im,jm);
SW4_bare_mean_may = zeros(im,jm); SW4_bare_mean_jun = zeros(im,jm);
SW4_bare_mean_jul = zeros(im,jm); SW4_bare_mean_aug = zeros(im,jm);
SW4_bare_mean_sep = zeros(im,jm); SW4_bare_mean_oct = zeros(im,jm);
SW4_bare_mean_nov = zeros(im,jm); SW4_bare_mean_dec = zeros(im,jm);

% Sum individual months over all years of simulation
for i = (num_spinup_yr+1):num_yrs
    SW4_bare_mean_jan = SW4_bare_mean_jan+temp_all(:,:,1+(i-1)*12);
    SW4_bare_mean_feb = SW4_bare_mean_feb+temp_all(:,:,2+(i-1)*12);
    SW4_bare_mean_mar = SW4_bare_mean_mar+temp_all(:,:,3+(i-1)*12);
    SW4_bare_mean_apr = SW4_bare_mean_apr+temp_all(:,:,4+(i-1)*12);
    SW4_bare_mean_may = SW4_bare_mean_may+temp_all(:,:,5+(i-1)*12);
    SW4_bare_mean_jun = SW4_bare_mean_jun+temp_all(:,:,6+(i-1)*12);
    SW4_bare_mean_jul = SW4_bare_mean_jul+temp_all(:,:,7+(i-1)*12);
    SW4_bare_mean_aug = SW4_bare_mean_aug+temp_all(:,:,8+(i-1)*12);
    SW4_bare_mean_sep = SW4_bare_mean_sep+temp_all(:,:,9+(i-1)*12);
    SW4_bare_mean_oct = SW4_bare_mean_oct+temp_all(:,:,10+(i-1)*12);
    SW4_bare_mean_nov = SW4_bare_mean_nov+temp_all(:,:,11+(i-1)*12);
    SW4_bare_mean_dec = SW4_bare_mean_dec+temp_all(:,:,12+(i-1)*12);
end

% Compute average monthly values [-]
num_yrs_avg = num_yrs-num_spinup_yr;
SW4_bare_mean_jan = SW4_bare_mean_jan/num_yrs_avg;
SW4_bare_mean_feb = SW4_bare_mean_feb/num_yrs_avg;
SW4_bare_mean_mar = SW4_bare_mean_mar/num_yrs_avg;
SW4_bare_mean_apr = SW4_bare_mean_apr/num_yrs_avg;
SW4_bare_mean_may = SW4_bare_mean_may/num_yrs_avg;
SW4_bare_mean_jun = SW4_bare_mean_jun/num_yrs_avg;
SW4_bare_mean_jul = SW4_bare_mean_jul/num_yrs_avg;
SW4_bare_mean_aug = SW4_bare_mean_aug/num_yrs_avg;
SW4_bare_mean_sep = SW4_bare_mean_sep/num_yrs_avg;
SW4_bare_mean_oct = SW4_bare_mean_oct/num_yrs_avg;
SW4_bare_mean_nov = SW4_bare_mean_nov/num_yrs_avg;
SW4_bare_mean_dec = SW4_bare_mean_dec/num_yrs_avg;

% Compute annual average values [-]
SW4_bare_mean_annual = (SW4_bare_mean_jan + SW4_bare_mean_feb+ SW4_bare_mean_mar + SW4_bare_mean_apr+ ...
                SW4_bare_mean_may + SW4_bare_mean_jun+ SW4_bare_mean_jul + SW4_bare_mean_aug+ ...
                SW4_bare_mean_sep + SW4_bare_mean_oct+ SW4_bare_mean_nov + SW4_bare_mean_dec)/12;

% Compute global total
SW4_bare = SW4_bare_mean_annual; SW4_bare(find(land_mask<=0))=0; 
SW4_bare=SW4_bare.*(land_area*1000000); %[-] * m^2
total_SW4_bare = sum(sum(SW4_bare(find(land_mask>0))))/(1000000*earth_land);%m^2/m^2
total_SW4_bare=total_SW4_bare; %[-]

% Compute time series        
% Monthly average
for i = 1:num_mnth 
   SW4_bare_time(i) = mean(mean(temp_all(:,:,i).*(temp_all(:,:,i)>=0)));
end
% Yearly average
for i = 1:num_yrs 
   SW4_bare_tyr(i) = mean(SW4_bare_time(12*(i-1)+1:i*12));   
end
%--------------------------------------------------------------------------
% Bare soil moisture layer 5 [m]
fid = fopen('/NoBackup/LSM_output/GHY_old/fort.1007','r','ieee-be');
if fid == -1
    disp(message)
end

layer_thickness = 0.8864;
%Initialize matrices
temp = zeros(im,jm);
temp_all = zeros(im,jm,num_mnth);

for i = 1:num_mnth
    skip = fread(fid,[4]);
    title_read = fread(fid,[80],'*char');
    temp = fread(fid,[im,jm],'float32','ieee-be'); %float64 for real*8; float32 for real*4
    temp_all(:,:,i)=temp/layer_thickness;
    skip = fread(fid,[4]);
end

SW5_bare_mean_jan = zeros(im,jm); SW5_bare_mean_feb = zeros(im,jm);
SW5_bare_mean_mar = zeros(im,jm); SW5_bare_mean_apr = zeros(im,jm);
SW5_bare_mean_may = zeros(im,jm); SW5_bare_mean_jun = zeros(im,jm);
SW5_bare_mean_jul = zeros(im,jm); SW5_bare_mean_aug = zeros(im,jm);
SW5_bare_mean_sep = zeros(im,jm); SW5_bare_mean_oct = zeros(im,jm);
SW5_bare_mean_nov = zeros(im,jm); SW5_bare_mean_dec = zeros(im,jm);

% Sum individual months over all years of simulation
for i = (num_spinup_yr+1):num_yrs
    SW5_bare_mean_jan = SW5_bare_mean_jan+temp_all(:,:,1+(i-1)*12);
    SW5_bare_mean_feb = SW5_bare_mean_feb+temp_all(:,:,2+(i-1)*12);
    SW5_bare_mean_mar = SW5_bare_mean_mar+temp_all(:,:,3+(i-1)*12);
    SW5_bare_mean_apr = SW5_bare_mean_apr+temp_all(:,:,4+(i-1)*12);
    SW5_bare_mean_may = SW5_bare_mean_may+temp_all(:,:,5+(i-1)*12);
    SW5_bare_mean_jun = SW5_bare_mean_jun+temp_all(:,:,6+(i-1)*12);
    SW5_bare_mean_jul = SW5_bare_mean_jul+temp_all(:,:,7+(i-1)*12);
    SW5_bare_mean_aug = SW5_bare_mean_aug+temp_all(:,:,8+(i-1)*12);
    SW5_bare_mean_sep = SW5_bare_mean_sep+temp_all(:,:,9+(i-1)*12);
    SW5_bare_mean_oct = SW5_bare_mean_oct+temp_all(:,:,10+(i-1)*12);
    SW5_bare_mean_nov = SW5_bare_mean_nov+temp_all(:,:,11+(i-1)*12);
    SW5_bare_mean_dec = SW5_bare_mean_dec+temp_all(:,:,12+(i-1)*12);
end

% Compute average monthly values [-]
num_yrs_avg = num_yrs-num_spinup_yr;
SW5_bare_mean_jan = SW5_bare_mean_jan/num_yrs_avg;
SW5_bare_mean_feb = SW5_bare_mean_feb/num_yrs_avg;
SW5_bare_mean_mar = SW5_bare_mean_mar/num_yrs_avg;
SW5_bare_mean_apr = SW5_bare_mean_apr/num_yrs_avg;
SW5_bare_mean_may = SW5_bare_mean_may/num_yrs_avg;
SW5_bare_mean_jun = SW5_bare_mean_jun/num_yrs_avg;
SW5_bare_mean_jul = SW5_bare_mean_jul/num_yrs_avg;
SW5_bare_mean_aug = SW5_bare_mean_aug/num_yrs_avg;
SW5_bare_mean_sep = SW5_bare_mean_sep/num_yrs_avg;
SW5_bare_mean_oct = SW5_bare_mean_oct/num_yrs_avg;
SW5_bare_mean_nov = SW5_bare_mean_nov/num_yrs_avg;
SW5_bare_mean_dec = SW5_bare_mean_dec/num_yrs_avg;

% Compute annual average values [-]
SW5_bare_mean_annual = (SW5_bare_mean_jan + SW5_bare_mean_feb+ SW5_bare_mean_mar + SW5_bare_mean_apr+ ...
                SW5_bare_mean_may + SW5_bare_mean_jun+ SW5_bare_mean_jul + SW5_bare_mean_aug+ ...
                SW5_bare_mean_sep + SW5_bare_mean_oct+ SW5_bare_mean_nov + SW5_bare_mean_dec)/12;

% Compute global total
SW5_bare = SW5_bare_mean_annual; SW5_bare(find(land_mask<=0))=0; 
SW5_bare=SW5_bare.*(land_area*1000000); %[-] * m^2
total_SW5_bare = sum(sum(SW5_bare(find(land_mask>0))))/(1000000*earth_land);%m^2/m^2
total_SW5_bare=total_SW5_bare; %[-]

% Compute time series        
% Monthly average
for i = 1:num_mnth 
   SW5_bare_time(i) = mean(mean(temp_all(:,:,i).*(temp_all(:,:,i)>=0)));
end
% Yearly average
for i = 1:num_yrs 
   SW5_bare_tyr(i) = mean(SW5_bare_time(12*(i-1)+1:i*12));   
end
%--------------------------------------------------------------------------
% Bare soil moisture layer 6 [m]
fid = fopen('/NoBackup/LSM_output/GHY_old/fort.1008','r','ieee-be');
if fid == -1
    disp(message)
end

layer_thickness = 1.529;
%Initialize matrices
temp = zeros(im,jm);
temp_all = zeros(im,jm,num_mnth);

for i = 1:num_mnth
    skip = fread(fid,[4]);
    title_read = fread(fid,[80],'*char');
    temp = fread(fid,[im,jm],'float32','ieee-be'); %float64 for real*8; float32 for real*4
    temp_all(:,:,i)=temp/layer_thickness;
    skip = fread(fid,[4]);
end

SW6_bare_mean_jan = zeros(im,jm); SW6_bare_mean_feb = zeros(im,jm);
SW6_bare_mean_mar = zeros(im,jm); SW6_bare_mean_apr = zeros(im,jm);
SW6_bare_mean_may = zeros(im,jm); SW6_bare_mean_jun = zeros(im,jm);
SW6_bare_mean_jul = zeros(im,jm); SW6_bare_mean_aug = zeros(im,jm);
SW6_bare_mean_sep = zeros(im,jm); SW6_bare_mean_oct = zeros(im,jm);
SW6_bare_mean_nov = zeros(im,jm); SW6_bare_mean_dec = zeros(im,jm);

% Sum individual months over all years of simulation
for i = (num_spinup_yr+1):num_yrs
    SW6_bare_mean_jan = SW6_bare_mean_jan+temp_all(:,:,1+(i-1)*12);
    SW6_bare_mean_feb = SW6_bare_mean_feb+temp_all(:,:,2+(i-1)*12);
    SW6_bare_mean_mar = SW6_bare_mean_mar+temp_all(:,:,3+(i-1)*12);
    SW6_bare_mean_apr = SW6_bare_mean_apr+temp_all(:,:,4+(i-1)*12);
    SW6_bare_mean_may = SW6_bare_mean_may+temp_all(:,:,5+(i-1)*12);
    SW6_bare_mean_jun = SW6_bare_mean_jun+temp_all(:,:,6+(i-1)*12);
    SW6_bare_mean_jul = SW6_bare_mean_jul+temp_all(:,:,7+(i-1)*12);
    SW6_bare_mean_aug = SW6_bare_mean_aug+temp_all(:,:,8+(i-1)*12);
    SW6_bare_mean_sep = SW6_bare_mean_sep+temp_all(:,:,9+(i-1)*12);
    SW6_bare_mean_oct = SW6_bare_mean_oct+temp_all(:,:,10+(i-1)*12);
    SW6_bare_mean_nov = SW6_bare_mean_nov+temp_all(:,:,11+(i-1)*12);
    SW6_bare_mean_dec = SW6_bare_mean_dec+temp_all(:,:,12+(i-1)*12);
end

% Compute average monthly values [-]
num_yrs_avg = num_yrs-num_spinup_yr;
SW6_bare_mean_jan = SW6_bare_mean_jan/num_yrs_avg;
SW6_bare_mean_feb = SW6_bare_mean_feb/num_yrs_avg;
SW6_bare_mean_mar = SW6_bare_mean_mar/num_yrs_avg;
SW6_bare_mean_apr = SW6_bare_mean_apr/num_yrs_avg;
SW6_bare_mean_may = SW6_bare_mean_may/num_yrs_avg;
SW6_bare_mean_jun = SW6_bare_mean_jun/num_yrs_avg;
SW6_bare_mean_jul = SW6_bare_mean_jul/num_yrs_avg;
SW6_bare_mean_aug = SW6_bare_mean_aug/num_yrs_avg;
SW6_bare_mean_sep = SW6_bare_mean_sep/num_yrs_avg;
SW6_bare_mean_oct = SW6_bare_mean_oct/num_yrs_avg;
SW6_bare_mean_nov = SW6_bare_mean_nov/num_yrs_avg;
SW6_bare_mean_dec = SW6_bare_mean_dec/num_yrs_avg;

% Compute annual average values [-]
SW6_bare_mean_annual = (SW6_bare_mean_jan + SW6_bare_mean_feb+ SW6_bare_mean_mar + SW6_bare_mean_apr+ ...
                SW6_bare_mean_may + SW6_bare_mean_jun+ SW6_bare_mean_jul + SW6_bare_mean_aug+ ...
                SW6_bare_mean_sep + SW6_bare_mean_oct+ SW6_bare_mean_nov + SW6_bare_mean_dec)/12;

% Compute global total
SW6_bare = SW6_bare_mean_annual; SW6_bare(find(land_mask<=0))=0; 
SW6_bare=SW6_bare.*(land_area*1000000); %[-] * m^2
total_SW6_bare = sum(sum(SW6_bare(find(land_mask>0))))/(1000000*earth_land);%m^2/m^2
total_SW6_bare=total_SW6_bare; %[-]

% Compute time series        
% Monthly average
for i = 1:num_mnth 
   SW6_bare_time(i) = mean(mean(temp_all(:,:,i).*(temp_all(:,:,i)>=0)));
end
% Yearly average
for i = 1:num_yrs 
   SW6_bare_tyr(i) = mean(SW6_bare_time(12*(i-1)+1:i*12));   
end
%--------------------------------------------------------------------------
% LSM runoff
fid = fopen('/NoBackup/LSM_output/GHY_old/fort.985','r','ieee-be');
if fid == -1
    disp(message)
end
%Initialize matrices
temp = zeros(im,jm);
temp_all = zeros(im,jm,num_mnth);

for i = 1:num_mnth
    skip = fread(fid,[4]);
    title_read = fread(fid,[80],'*char');
    temp = fread(fid,[im,jm],'float32','ieee-be');
    temp_all(:,:,i)=temp;
    skip = fread(fid,[4]);  
end
% underground runoff
fid = fopen('/NoBackup/LSM_output/GHY_old/fort.986','r','ieee-be');
if fid == -1
    disp(message)
end

%Initialize matrices
temp2 = zeros(im,jm);
temp2_all = zeros(im,jm,num_mnth);

for i = 1:num_mnth
    skip = fread(fid,[4]);
    title_read = fread(fid,[80],'*char');
    temp2 = fread(fid,[im,jm],'float32','ieee-be');
    temp2_all(:,:,i)=temp2;
    skip = fread(fid,[4]);
end

% Sum individual months over all years of simulation
Rs_mean_jan = zeros(im,jm); Rs_mean_feb = zeros(im,jm);
Rs_mean_mar = zeros(im,jm); Rs_mean_apr = zeros(im,jm);
Rs_mean_may = zeros(im,jm); Rs_mean_jun = zeros(im,jm);
Rs_mean_jul = zeros(im,jm); Rs_mean_aug = zeros(im,jm);
Rs_mean_sep = zeros(im,jm); Rs_mean_oct = zeros(im,jm);
Rs_mean_nov = zeros(im,jm); Rs_mean_dec = zeros(im,jm);

Ru_mean_jan = zeros(im,jm); Ru_mean_feb = zeros(im,jm);
Ru_mean_mar = zeros(im,jm); Ru_mean_apr = zeros(im,jm);
Ru_mean_may = zeros(im,jm); Ru_mean_jun = zeros(im,jm);
Ru_mean_jul = zeros(im,jm); Ru_mean_aug = zeros(im,jm);
Ru_mean_sep = zeros(im,jm); Ru_mean_oct = zeros(im,jm);
Ru_mean_nov = zeros(im,jm); Ru_mean_dec = zeros(im,jm);

for i = (num_spinup_yr+1):num_yrs
    Rs_mean_jan = Rs_mean_jan+temp_all(:,:,1+(i-1)*12);
    Rs_mean_feb = Rs_mean_feb+temp_all(:,:,2+(i-1)*12);
    Rs_mean_mar = Rs_mean_mar+temp_all(:,:,3+(i-1)*12);
    Rs_mean_apr = Rs_mean_apr+temp_all(:,:,4+(i-1)*12);
    Rs_mean_may = Rs_mean_may+temp_all(:,:,5+(i-1)*12);
    Rs_mean_jun = Rs_mean_jun+temp_all(:,:,6+(i-1)*12);
    Rs_mean_jul = Rs_mean_jul+temp_all(:,:,7+(i-1)*12);
    Rs_mean_aug = Rs_mean_aug+temp_all(:,:,8+(i-1)*12);
    Rs_mean_sep = Rs_mean_sep+temp_all(:,:,9+(i-1)*12);
    Rs_mean_oct = Rs_mean_oct+temp_all(:,:,10+(i-1)*12);
    Rs_mean_nov = Rs_mean_nov+temp_all(:,:,11+(i-1)*12);
    Rs_mean_dec = Rs_mean_dec+temp_all(:,:,12+(i-1)*12);
    
    Ru_mean_jan = Ru_mean_jan+temp2_all(:,:,1+(i-1)*12);
    Ru_mean_feb = Ru_mean_feb+temp2_all(:,:,2+(i-1)*12);
    Ru_mean_mar = Ru_mean_mar+temp2_all(:,:,3+(i-1)*12);
    Ru_mean_apr = Ru_mean_apr+temp2_all(:,:,4+(i-1)*12);
    Ru_mean_may = Ru_mean_may+temp2_all(:,:,5+(i-1)*12);
    Ru_mean_jun = Ru_mean_jun+temp2_all(:,:,6+(i-1)*12);
    Ru_mean_jul = Ru_mean_jul+temp2_all(:,:,7+(i-1)*12);
    Ru_mean_aug = Ru_mean_aug+temp2_all(:,:,8+(i-1)*12);
    Ru_mean_sep = Ru_mean_sep+temp2_all(:,:,9+(i-1)*12);
    Ru_mean_oct = Ru_mean_oct+temp2_all(:,:,10+(i-1)*12);
    Ru_mean_nov = Ru_mean_nov+temp2_all(:,:,11+(i-1)*12);
    Ru_mean_dec = Ru_mean_dec+temp2_all(:,:,12+(i-1)*12);
end

% Compute average monthly values [kg/m2/month]
num_yrs_avg = num_yrs-num_spinup_yr;
Rs_mean_jan = Rs_mean_jan/num_yrs_avg;
Rs_mean_feb = Rs_mean_feb/num_yrs_avg;
Rs_mean_mar = Rs_mean_mar/num_yrs_avg;
Rs_mean_apr = Rs_mean_apr/num_yrs_avg;
Rs_mean_may = Rs_mean_may/num_yrs_avg;
Rs_mean_jun = Rs_mean_jun/num_yrs_avg;
Rs_mean_jul = Rs_mean_jul/num_yrs_avg;
Rs_mean_aug = Rs_mean_aug/num_yrs_avg;
Rs_mean_sep = Rs_mean_sep/num_yrs_avg;
Rs_mean_oct = Rs_mean_oct/num_yrs_avg;
Rs_mean_nov = Rs_mean_nov/num_yrs_avg;
Rs_mean_dec = Rs_mean_dec/num_yrs_avg;

Ru_mean_jan = Ru_mean_jan/num_yrs_avg;
Ru_mean_feb = Ru_mean_feb/num_yrs_avg;
Ru_mean_mar = Ru_mean_mar/num_yrs_avg;
Ru_mean_apr = Ru_mean_apr/num_yrs_avg;
Ru_mean_may = Ru_mean_may/num_yrs_avg;
Ru_mean_jun = Ru_mean_jun/num_yrs_avg;
Ru_mean_jul = Ru_mean_jul/num_yrs_avg;
Ru_mean_aug = Ru_mean_aug/num_yrs_avg;
Ru_mean_sep = Ru_mean_sep/num_yrs_avg;
Ru_mean_oct = Ru_mean_oct/num_yrs_avg;
Ru_mean_nov = Ru_mean_nov/num_yrs_avg;
Ru_mean_dec = Ru_mean_dec/num_yrs_avg;
          
% Compute annual average values [kg/m2/yr]
Rs_mean_annual = (Rs_mean_jan + Rs_mean_feb+ Rs_mean_mar + Rs_mean_apr+ ...
                Rs_mean_may + Rs_mean_jun+ Rs_mean_jul + Rs_mean_aug+ ...
                Rs_mean_sep + Rs_mean_oct+ Rs_mean_nov + Rs_mean_dec)/12;
Ru_mean_annual = (Ru_mean_jan + Ru_mean_feb+ Ru_mean_mar + Ru_mean_apr+ ...
                Ru_mean_may + Ru_mean_jun+ Ru_mean_jul + Ru_mean_aug+ ...
                Ru_mean_sep + Ru_mean_oct+ Ru_mean_nov + Ru_mean_dec)/12;
Rs_mean10yr  = mean(mean(Rs_mean_annual));
Ru_mean10yr  = mean(mean(Ru_mean_annual));
Rsu_mean10yr = Rs_mean10yr+Ru_mean10yr;

% Compute global total
Rs = Rs_mean_annual; Rs(find(land_mask<=0))=0; 
Rs=Rs.*(land_area*1000000)/rhow; %kg/m^2/yr * m^2 / (kg/m3)
total_Rs = sum(sum(Rs(find(land_mask>0))))/(1000000*earth_land);%m/yr;
total_Rs=1000*total_Rs; %mm/yr

Ru = Ru_mean_annual; Ru(find(land_mask<=0))=0; 
Ru=Ru.*(land_area*1000000)/rhow; %kg/m^2/yr * m^2 / (kg/m3)
total_Ru = sum(sum(Ru(find(land_mask>0))))/(1000000*earth_land);%m/yr;
total_Ru=1000*total_Ru; %mm/yr


% Compute time series        
% Monthly average
for i = 1:num_mnth
   Rs_time(i) = mean(mean(temp_all(:,:,i).*(temp_all(:,:,i)>=0)))*mm_fromkg_per_m2;
   Ru_time(i) = mean(mean(temp2_all(:,:,i).*(temp2_all(:,:,i)>=0)))*mm_fromkg_per_m2;
   Rsu_time(i)= Rs_time(i)+Ru_time(i);
end
% Yearly average
for i = 1:num_yrs
   Rs_tyr(i) = mean(Rs_time(12*(i-1)+1:i*12));
   Ru_tyr(i) = mean(Ru_time(12*(i-1)+1:i*12));
   Rsu_tyr(i)= Rs_tyr(i)+Ru_tyr(i);    
end
%--------------------------------------------------------------------------
% Gross Primary Productivity [kg C/m2]
fid = fopen('/NoBackup/LSM_output/GHY_old/fort.987','r','ieee-be');
if fid == -1
    disp(message)
end

%Initialize matrices
temp = zeros(im,jm);
temp_all = zeros(im,jm,num_mnth);

for i = 1:num_mnth
    skip = fread(fid,[4]);
    title_read = fread(fid,[80],'*char');
    temp = fread(fid,[im,jm],'float32','ieee-be'); %float64 for real*8; float32 for real*4
    temp_all(:,:,i)=temp;
    skip = fread(fid,[4]);
end
%--------------------------------------------------------------------------
GPP_mean_jan = zeros(im,jm); GPP_mean_feb = zeros(im,jm);
GPP_mean_mar = zeros(im,jm); GPP_mean_apr = zeros(im,jm);
GPP_mean_may = zeros(im,jm); GPP_mean_jun = zeros(im,jm);
GPP_mean_jul = zeros(im,jm); GPP_mean_aug = zeros(im,jm);
GPP_mean_sep = zeros(im,jm); GPP_mean_oct = zeros(im,jm);
GPP_mean_nov = zeros(im,jm); GPP_mean_dec = zeros(im,jm);

% Sum individual months over all years of simulation
for i = (num_spinup_yr+1):num_yrs
    GPP_mean_jan = GPP_mean_jan+temp_all(:,:,1+(i-1)*12);
    GPP_mean_feb = GPP_mean_feb+temp_all(:,:,2+(i-1)*12);
    GPP_mean_mar = GPP_mean_mar+temp_all(:,:,3+(i-1)*12);
    GPP_mean_apr = GPP_mean_apr+temp_all(:,:,4+(i-1)*12);
    GPP_mean_may = GPP_mean_may+temp_all(:,:,5+(i-1)*12);
    GPP_mean_jun = GPP_mean_jun+temp_all(:,:,6+(i-1)*12);
    GPP_mean_jul = GPP_mean_jul+temp_all(:,:,7+(i-1)*12);
    GPP_mean_aug = GPP_mean_aug+temp_all(:,:,8+(i-1)*12);
    GPP_mean_sep = GPP_mean_sep+temp_all(:,:,9+(i-1)*12);
    GPP_mean_oct = GPP_mean_oct+temp_all(:,:,10+(i-1)*12);
    GPP_mean_nov = GPP_mean_nov+temp_all(:,:,11+(i-1)*12);
    GPP_mean_dec = GPP_mean_dec+temp_all(:,:,12+(i-1)*12);
end

% Compute average monthly values [kg C /month]
num_yrs_avg = num_yrs-num_spinup_yr;
GPP_mean_jan = GPP_mean_jan/num_yrs_avg;
GPP_mean_feb = GPP_mean_feb/num_yrs_avg;
GPP_mean_mar = GPP_mean_mar/num_yrs_avg;
GPP_mean_apr = GPP_mean_apr/num_yrs_avg;
GPP_mean_may = GPP_mean_may/num_yrs_avg;
GPP_mean_jun = GPP_mean_jun/num_yrs_avg;
GPP_mean_jul = GPP_mean_jul/num_yrs_avg;
GPP_mean_aug = GPP_mean_aug/num_yrs_avg;
GPP_mean_sep = GPP_mean_sep/num_yrs_avg;
GPP_mean_oct = GPP_mean_oct/num_yrs_avg;
GPP_mean_nov = GPP_mean_nov/num_yrs_avg;
GPP_mean_dec = GPP_mean_dec/num_yrs_avg;

% Compute annual average values [kg C/yr]
GPP_mean_annual = (GPP_mean_jan + GPP_mean_feb+ GPP_mean_mar + GPP_mean_apr+ ...
                GPP_mean_may + GPP_mean_jun+ GPP_mean_jul + GPP_mean_aug+ ...
                GPP_mean_sep + GPP_mean_oct+ GPP_mean_nov + GPP_mean_dec);

% Compute global total carbon
gpp = GPP_mean_annual; gpp(find(land_mask<=0))=0; 
gpp=gpp.*land_area*1000000;%kg C/m^2/yr * m^2
total_gpp = sum(sum(gpp(find(land_mask>0)))); %kg C  /yr
total_gpp = total_gpp/1e12; %Pg C   (10^15 grams)

% Compute time series        
% Monthly average
for i = 1:num_mnth 
   GPP_time(i) = mean(mean(temp_all(:,:,i).*(temp_all(:,:,i)>=0)));
end
% Yearly average
for i = 1:num_yrs 
   GPP_tyr(i) = mean(GPP_time(12*(i-1)+1:i*12));   
end

%--------------------------------------------------------------------------
% Autotrophic respiration [kg C/m2]
fid = fopen('/NoBackup/LSM_output/GHY_old/fort.988','r','ieee-be');
if fid == -1
    disp(message)
end

%Initialize matrices
temp = zeros(im,jm);
temp_all = zeros(im,jm,num_mnth);

for i = 1:num_mnth
    skip = fread(fid,[4]);
    title_read = fread(fid,[80],'*char');
    temp = fread(fid,[im,jm],'float32','ieee-be'); %float64 for real*8; float32 for real*4
    temp_all(:,:,i)=temp;
    skip = fread(fid,[4]);    
end

R_auto_mean_jan = zeros(im,jm); R_auto_mean_feb = zeros(im,jm);
R_auto_mean_mar = zeros(im,jm); R_auto_mean_apr = zeros(im,jm);
R_auto_mean_may = zeros(im,jm); R_auto_mean_jun = zeros(im,jm);
R_auto_mean_jul = zeros(im,jm); R_auto_mean_aug = zeros(im,jm);
R_auto_mean_sep = zeros(im,jm); R_auto_mean_oct = zeros(im,jm);
R_auto_mean_nov = zeros(im,jm); R_auto_mean_dec = zeros(im,jm);

% Sum individual months over all years of simulation
for i = (num_spinup_yr+1):num_yrs
    R_auto_mean_jan = R_auto_mean_jan+temp_all(:,:,1+(i-1)*12);
    R_auto_mean_feb = R_auto_mean_feb+temp_all(:,:,2+(i-1)*12);
    R_auto_mean_mar = R_auto_mean_mar+temp_all(:,:,3+(i-1)*12);
    R_auto_mean_apr = R_auto_mean_apr+temp_all(:,:,4+(i-1)*12);
    R_auto_mean_may = R_auto_mean_may+temp_all(:,:,5+(i-1)*12);
    R_auto_mean_jun = R_auto_mean_jun+temp_all(:,:,6+(i-1)*12);
    R_auto_mean_jul = R_auto_mean_jul+temp_all(:,:,7+(i-1)*12);
    R_auto_mean_aug = R_auto_mean_aug+temp_all(:,:,8+(i-1)*12);
    R_auto_mean_sep = R_auto_mean_sep+temp_all(:,:,9+(i-1)*12);
    R_auto_mean_oct = R_auto_mean_oct+temp_all(:,:,10+(i-1)*12);
    R_auto_mean_nov = R_auto_mean_nov+temp_all(:,:,11+(i-1)*12);
    R_auto_mean_dec = R_auto_mean_dec+temp_all(:,:,12+(i-1)*12);    
end

% Compute average monthly values [kg C /month]
num_yrs_avg = num_yrs-num_spinup_yr;
R_auto_mean_jan = R_auto_mean_jan/num_yrs_avg;
R_auto_mean_feb = R_auto_mean_feb/num_yrs_avg;
R_auto_mean_mar = R_auto_mean_mar/num_yrs_avg;
R_auto_mean_apr = R_auto_mean_apr/num_yrs_avg;
R_auto_mean_may = R_auto_mean_may/num_yrs_avg;
R_auto_mean_jun = R_auto_mean_jun/num_yrs_avg;
R_auto_mean_jul = R_auto_mean_jul/num_yrs_avg;
R_auto_mean_aug = R_auto_mean_aug/num_yrs_avg;
R_auto_mean_sep = R_auto_mean_sep/num_yrs_avg;
R_auto_mean_oct = R_auto_mean_oct/num_yrs_avg;
R_auto_mean_nov = R_auto_mean_nov/num_yrs_avg;
R_auto_mean_dec = R_auto_mean_dec/num_yrs_avg;

% Compute annual average values [kg C/yr]        
R_auto_mean_annual = (R_auto_mean_jan + R_auto_mean_feb+ R_auto_mean_mar + R_auto_mean_apr+ ...
                R_auto_mean_may + R_auto_mean_jun+ R_auto_mean_jul + R_auto_mean_aug+ ...
                R_auto_mean_sep + R_auto_mean_oct+ R_auto_mean_nov + R_auto_mean_dec);

% Compute global total carbon
R_auto = R_auto_mean_annual; R_auto(find(land_mask<=0))=0; %kg C/m^2
R_auto=R_auto.*land_area*1000000;%kg C/m^2/yr * m^2
total_R_auto = sum(sum(R_auto(find(land_mask>0)))); %kg C  /yr
total_R_auto = total_R_auto/1e12; %Pg C   (10^15 grams)

% Compute time series 
% Monthly average
for i = 1:num_mnth
   R_auto_time(i) = mean(mean(temp_all(:,:,i).*(temp_all(:,:,i)>=0)));
end
% Yearly average
for i = 1:num_yrs
   R_auto_tyr(i) = mean(R_auto_time(12*(i-1)+1:i*12));   
end
%--------------------------------------------------------------------------
% Compute net primary productivity from GPP and R_auto
NPP_mean_annual=GPP_mean_annual-R_auto_mean_annual;
total_NPP =total_gpp-total_R_auto;
NPP_time=GPP_time-R_auto_time;
NPP_tyr=GPP_tyr-R_auto_tyr;

NPP_mean_jan = GPP_mean_jan - R_auto_mean_jan;
NPP_mean_feb = GPP_mean_feb - R_auto_mean_feb;
NPP_mean_mar = GPP_mean_mar - R_auto_mean_mar;
NPP_mean_apr = GPP_mean_apr - R_auto_mean_apr;
NPP_mean_may = GPP_mean_may - R_auto_mean_may;
NPP_mean_jun = GPP_mean_jun - R_auto_mean_jun;
NPP_mean_jul = GPP_mean_jul - R_auto_mean_jul;
NPP_mean_aug = GPP_mean_aug - R_auto_mean_aug;
NPP_mean_sep = GPP_mean_sep - R_auto_mean_sep;
NPP_mean_oct = GPP_mean_oct - R_auto_mean_oct;
NPP_mean_nov = GPP_mean_nov - R_auto_mean_nov;
NPP_mean_dec = GPP_mean_dec - R_auto_mean_dec;
%--------------------------------------------------------------------------
% soil respiration [kg C/m2]
fid = fopen('/NoBackup/LSM_output/GHY_old/fort.989','r','ieee-be');
if fid == -1
    disp(message)
end

%Initialize matrices
temp = zeros(im,jm);
temp_all = zeros(im,jm,num_mnth);

for i = 1:num_mnth
    skip = fread(fid,[4]);
    title_read = fread(fid,[80],'*char');
    temp = fread(fid,[im,jm],'float32','ieee-be'); %float64 for real*8; float32 for real*4
    temp_all(:,:,i)=temp;
    skip = fread(fid,[4]); 
end
R_soil_mean_jan = zeros(im,jm); R_soil_mean_feb = zeros(im,jm);
R_soil_mean_mar = zeros(im,jm); R_soil_mean_apr = zeros(im,jm);
R_soil_mean_may = zeros(im,jm); R_soil_mean_jun = zeros(im,jm);
R_soil_mean_jul = zeros(im,jm); R_soil_mean_aug = zeros(im,jm);
R_soil_mean_sep = zeros(im,jm); R_soil_mean_oct = zeros(im,jm);
R_soil_mean_nov = zeros(im,jm); R_soil_mean_dec = zeros(im,jm);

% Sum individual months over all years of simulation
for i = (num_spinup_yr+1):num_yrs
    R_soil_mean_jan = R_soil_mean_jan+temp_all(:,:,1+(i-1)*12);
    R_soil_mean_feb = R_soil_mean_feb+temp_all(:,:,2+(i-1)*12);
    R_soil_mean_mar = R_soil_mean_mar+temp_all(:,:,3+(i-1)*12);
    R_soil_mean_apr = R_soil_mean_apr+temp_all(:,:,4+(i-1)*12);
    R_soil_mean_may = R_soil_mean_may+temp_all(:,:,5+(i-1)*12);
    R_soil_mean_jun = R_soil_mean_jun+temp_all(:,:,6+(i-1)*12);
    R_soil_mean_jul = R_soil_mean_jul+temp_all(:,:,7+(i-1)*12);
    R_soil_mean_aug = R_soil_mean_aug+temp_all(:,:,8+(i-1)*12);
    R_soil_mean_sep = R_soil_mean_sep+temp_all(:,:,9+(i-1)*12);
    R_soil_mean_oct = R_soil_mean_oct+temp_all(:,:,10+(i-1)*12);
    R_soil_mean_nov = R_soil_mean_nov+temp_all(:,:,11+(i-1)*12);
    R_soil_mean_dec = R_soil_mean_dec+temp_all(:,:,12+(i-1)*12);
end

% Compute average monthly values [kg C /month]
num_yrs_avg = num_yrs-num_spinup_yr;
R_soil_mean_jan = R_soil_mean_jan/num_yrs_avg;
R_soil_mean_feb = R_soil_mean_feb/num_yrs_avg;
R_soil_mean_mar = R_soil_mean_mar/num_yrs_avg;
R_soil_mean_apr = R_soil_mean_apr/num_yrs_avg;
R_soil_mean_may = R_soil_mean_may/num_yrs_avg;
R_soil_mean_jun = R_soil_mean_jun/num_yrs_avg;
R_soil_mean_jul = R_soil_mean_jul/num_yrs_avg;
R_soil_mean_aug = R_soil_mean_aug/num_yrs_avg;
R_soil_mean_sep = R_soil_mean_sep/num_yrs_avg;
R_soil_mean_oct = R_soil_mean_oct/num_yrs_avg;
R_soil_mean_nov = R_soil_mean_nov/num_yrs_avg;
R_soil_mean_dec = R_soil_mean_dec/num_yrs_avg;

% Compute annual average values [kg C/yr]        
R_soil_mean_annual = (R_soil_mean_jan + R_soil_mean_feb+ R_soil_mean_mar + R_soil_mean_apr+ ...
                R_soil_mean_may + R_soil_mean_jun+ R_soil_mean_jul + R_soil_mean_aug+ ...
                R_soil_mean_sep + R_soil_mean_oct+ R_soil_mean_nov + R_soil_mean_dec);

% Compute global total carbon
R_soil = R_soil_mean_annual; R_soil(find(land_mask<=0))=0; %kg C/m^2
R_soil=R_soil.*land_area*1000000;%kg C/m^2/yr * m^2
total_R_soil = sum(sum(R_soil(find(land_mask>0)))); %kg C  /yr
total_R_soil = total_R_soil/1e12; %Pg C   (10^15 grams)

% Compute time series 
% Monthly average
for i = 1:num_mnth
   R_soil_time(i) = mean(mean(temp_all(:,:,i).*(temp_all(:,:,i)>=0)));
end
% Yearly average
for i = 1:num_yrs
   R_soil_tyr(i) = mean(R_soil_time(12*(i-1)+1:i*12));   
end

% Compute net ecosystem exchange from GPP, R_auto, and R_soil
NEE_mean_annual=-(GPP_mean_annual-R_auto_mean_annual-R_soil_mean_annual);
total_NEE =-(total_gpp-total_R_auto-total_R_soil);
NEE_time=-(GPP_time-R_auto_time-R_soil_time);
NEE_tyr=-(GPP_tyr-R_auto_tyr-R_soil_tyr);

% %--------------------------------------------------------------------------
% % 3-hourly data with resolution of 4 deg latitude x 5 deg longitude
% % functions month_rain and month_snow return average intensity for 
% % 12 months of a year with units kg/m^2/s
% 
% rain_min = 0; rain_max = 0.03; % Rainfall rate [kg/m2/s]  % GSWP
% snow_min = 0; snow_max = 0.007;% Snowfall rate [kg/m2/s] % GSWP
% monthly_rain_avg=zeros(im,jm,12);
% monthly_snow_avg=zeros(im,jm,12);
% 
% year = '1986'; hr_feb = 28*8;
% month_rain_time(:,:,1:12)=month_rain(year,hr_feb,land_mask,im,jm,rain_min,rain_max);
% month_snow_time(:,:,1:12)=month_snow(year,hr_feb,land_mask,im,jm,snow_min,snow_max);
% monthly_rain_avg=monthly_rain_avg + month_rain_time(:,:,1:12);
% monthly_snow_avg=monthly_snow_avg + month_snow_time(:,:,1:12);
% 
% year = '1987'; hr_feb = 28*8;
% month_rain_time(:,:,13:24)=month_rain(year,hr_feb,land_mask,im,jm,rain_min,rain_max);
% month_snow_time(:,:,13:24)=month_snow(year,hr_feb,land_mask,im,jm,snow_min,snow_max);
% monthly_rain_avg=monthly_rain_avg + month_rain_time(:,:,13:24);
% monthly_snow_avg=monthly_snow_avg + month_snow_time(:,:,13:24);
% 
% year = '1988'; hr_feb = 29*8;
% month_rain_time(:,:,25:36)=month_rain(year,hr_feb,land_mask,im,jm,rain_min,rain_max);
% month_snow_time(:,:,25:36)=month_snow(year,hr_feb,land_mask,im,jm,snow_min,snow_max);
% monthly_rain_avg=monthly_rain_avg + month_rain_time(:,:,25:36);
% monthly_snow_avg=monthly_snow_avg + month_snow_time(:,:,25:36);
% 
% year = '1989'; hr_feb = 28*8;
% month_rain_time(:,:,37:48)=month_rain(year,hr_feb,land_mask,im,jm,rain_min,rain_max);
% month_snow_time(:,:,37:48)=month_snow(year,hr_feb,land_mask,im,jm,snow_min,snow_max);
% monthly_rain_avg=monthly_rain_avg + month_rain_time(:,:,37:48);
% monthly_snow_avg=monthly_snow_avg + month_snow_time(:,:,37:48);
% 
% year = '1990'; hr_feb = 28*8;
% month_rain_time(:,:,49:60)=month_rain(year,hr_feb,land_mask,im,jm,rain_min,rain_max);
% month_snow_time(:,:,49:60)=month_snow(year,hr_feb,land_mask,im,jm,snow_min,snow_max);
% monthly_rain_avg=monthly_rain_avg + month_rain_time(:,:,49:60);
% monthly_snow_avg=monthly_snow_avg + month_snow_time(:,:,49:60);
% 
% year = '1991'; hr_feb = 28*8;
% month_rain_time(:,:,61:72)=month_rain(year,hr_feb,land_mask,im,jm,rain_min,rain_max);
% month_snow_time(:,:,61:72)=month_snow(year,hr_feb,land_mask,im,jm,snow_min,snow_max);
% monthly_rain_avg=monthly_rain_avg + month_rain_time(:,:,61:72);
% monthly_snow_avg=monthly_snow_avg + month_snow_time(:,:,61:72);
% 
% year = '1992'; hr_feb = 29*8;
% month_rain_time(:,:,73:84)=month_rain(year,hr_feb,land_mask,im,jm,rain_min,rain_max);
% month_snow_time(:,:,73:84)=month_snow(year,hr_feb,land_mask,im,jm,snow_min,snow_max);
% monthly_rain_avg=monthly_rain_avg + month_rain_time(:,:,73:84);
% monthly_snow_avg=monthly_snow_avg + month_snow_time(:,:,73:84);
% 
% year = '1993'; hr_feb = 28*8;
% month_rain_time(:,:,85:96)=month_rain(year,hr_feb,land_mask,im,jm,rain_min,rain_max);
% month_snow_time(:,:,85:96)=month_snow(year,hr_feb,land_mask,im,jm,snow_min,snow_max);
% monthly_rain_avg=monthly_rain_avg + month_rain_time(:,:,85:96);
% monthly_snow_avg=monthly_snow_avg + month_snow_time(:,:,85:96);
% 
% year = '1994'; hr_feb = 28*8;
% month_rain_time(:,:,97:108)=month_rain(year,hr_feb,land_mask,im,jm,rain_min,rain_max);
% month_snow_time(:,:,97:108)=month_snow(year,hr_feb,land_mask,im,jm,snow_min,snow_max);
% monthly_rain_avg=monthly_rain_avg + month_rain_time(:,:,97:108);
% monthly_snow_avg=monthly_snow_avg + month_snow_time(:,:,97:108);
% 
% year = '1995'; hr_feb = 28*8;
% month_rain_time(:,:,109:120)=month_rain(year,hr_feb,land_mask,im,jm,rain_min,rain_max);
% month_snow_time(:,:,109:120)=month_snow(year,hr_feb,land_mask,im,jm,snow_min,snow_max);
% monthly_rain_avg=monthly_rain_avg + month_rain_time(:,:,109:120);
% monthly_snow_avg=monthly_snow_avg + month_snow_time(:,:,109:120);
% 
% % Average over the num of years of data years [kg/m^2/s]
% monthly_rain_avg  = monthly_rain_avg /num_yrs_data;
% monthly_snow_avg  = monthly_snow_avg /num_yrs_data;
% 
% % Global monthly rainfall averages [kg/m^2/s]
% for i = 1:120
%     Rain_time(i)=mean(mean(month_rain_time(:,:,i)))/num_yrs_data;
%     Snow_time(i)=mean(mean(month_snow_time(:,:,i)))/num_yrs_data;
% end
% 
% % Add rain and snow together for total precip. [kg/m^2/month]
% avg_num_days = [31 28.2 31 30 31 30 31 31 30 31 30 31];
% for i=1:12
%     dt = 86400*avg_num_days(i); % seconds per month
%     monthly_precip_avg(:,:,i)=(monthly_rain_avg(:,:,i) +...
%                   monthly_snow_avg(:,:,i))*dt;
% end
% 
% % Compute monthly average dryness index
% monthly_dryness  = monthly_Epot_avg./monthly_precip_avg;
% 
% %--------------------------------------------------------------------------
% seasonal_Epot_avg= zeros(im,jm,4);
% seasonal_precip  = zeros(im,jm,4);
% seasonal_dryness = zeros(im,jm,4);
% 
% annual_dryness = zeros(im,jm);
% annual_precip  = zeros(im,jm);
% annual_Epot_avg= zeros(im,jm);
% 
% seasonal_R_soil_avg(:,:,1) =  (R_soil_mean_dec+ R_soil_mean_jan + R_soil_mean_feb) /3;
% seasonal_R_soil_avg(:,:,2) =  (R_soil_mean_mar+ R_soil_mean_apr + R_soil_mean_may) /3;
% seasonal_R_soil_avg(:,:,3) =  (R_soil_mean_jun+ R_soil_mean_jul + R_soil_mean_aug) /3;
% seasonal_R_soil_avg(:,:,4) =  (R_soil_mean_sep+ R_soil_mean_oct + R_soil_mean_nov) /3;
% 
% seasonal_GPP_avg(:,:,1) =  (GPP_mean_dec+ GPP_mean_jan + GPP_mean_feb) /3;
% seasonal_GPP_avg(:,:,2) =  (GPP_mean_mar+ GPP_mean_apr + GPP_mean_may) /3;
% seasonal_GPP_avg(:,:,3) =  (GPP_mean_jun+ GPP_mean_jul + GPP_mean_aug) /3;
% seasonal_GPP_avg(:,:,4) =  (GPP_mean_sep+ GPP_mean_oct + GPP_mean_nov) /3;
% 
% seasonal_NPP_avg(:,:,1) =  (NPP_mean_dec+ NPP_mean_jan + NPP_mean_feb) /3;
% seasonal_NPP_avg(:,:,2) =  (NPP_mean_mar+ NPP_mean_apr + NPP_mean_may) /3;
% seasonal_NPP_avg(:,:,3) =  (NPP_mean_jun+ NPP_mean_jul + NPP_mean_aug) /3;
% seasonal_NPP_avg(:,:,4) =  (NPP_mean_sep+ NPP_mean_oct + NPP_mean_nov) /3;
% 
% seasonal_Epot_avg(:,:,1) = (monthly_Epot_avg(:,:,12)+ monthly_Epot_avg(:,:,1) + monthly_Epot_avg(:,:,2)) /3;
% seasonal_Epot_avg(:,:,2) = (monthly_Epot_avg(:,:,3) + monthly_Epot_avg(:,:,4) + monthly_Epot_avg(:,:,5)) /3;
% seasonal_Epot_avg(:,:,3) = (monthly_Epot_avg(:,:,6) + monthly_Epot_avg(:,:,7) + monthly_Epot_avg(:,:,8)) /3;
% seasonal_Epot_avg(:,:,4) = (monthly_Epot_avg(:,:,9) + monthly_Epot_avg(:,:,10)+ monthly_Epot_avg(:,:,11))/3;
% 
% seasonal_precip(:,:,1) = (monthly_precip_avg(:,:,12) + monthly_precip_avg(:,:,1)+monthly_precip_avg(:,:,2))/3;
% seasonal_precip(:,:,2) = (monthly_precip_avg(:,:,3)  + monthly_precip_avg(:,:,4)+monthly_precip_avg(:,:,5))/3;
% seasonal_precip(:,:,3) = (monthly_precip_avg(:,:,6)  + monthly_precip_avg(:,:,7)+monthly_precip_avg(:,:,8))/3;
% seasonal_precip(:,:,4) = (monthly_precip_avg(:,:,9)  + monthly_precip_avg(:,:,10)+monthly_precip_avg(:,:,11))/3;
% 
% seasonal_dryness(:,:,1) = (monthly_Epot_avg(:,:,1)./monthly_precip_avg(:,:,1));
% seasonal_dryness(:,:,2) = (monthly_Epot_avg(:,:,2)./monthly_precip_avg(:,:,2));
% seasonal_dryness(:,:,3) = (monthly_Epot_avg(:,:,3)./monthly_precip_avg(:,:,3));
% seasonal_dryness(:,:,4) = (monthly_Epot_avg(:,:,4)./monthly_precip_avg(:,:,4));
% 
% annual_precip(:,:)  = (seasonal_precip(:,:,1)  + seasonal_precip(:,:,2)  + seasonal_precip(:,:,3)  + seasonal_precip(:,:,4))  /4;
% annual_Epot_avg(:,:)= (seasonal_Epot_avg(:,:,1)+ seasonal_Epot_avg(:,:,2)+ seasonal_Epot_avg(:,:,3)+ seasonal_Epot_avg(:,:,4))/4;
% annual_dryness(:,:) = annual_Epot_avg(:,:)./annual_precip;
% 
% annual_dryness(find(land_mask==0))=NaN;
% annual_precip(find(land_mask==0))=NaN;
% annual_Epot_avg(find(land_mask==0))=NaN;

%--------------------------------------------------------------------------
% Compute evapotranspiration fractions

E_frac_can=total_Ecan/total_Etot;
E_frac_tran=total_Etran/total_Etot;
E_frac_bare=total_Ebare/total_Etot;
E_frac_evapveg=1-E_frac_can-E_frac_tran-E_frac_bare;

E_evapveg_mean_annual=E_tot_mean_annual-E_can_mean_annual-E_tran_mean_annual-E_bare_mean_annual;
% Compute global total
Evg = E_evapveg_mean_annual; Evg(find(land_mask<=0))=0; 
total_Evg = mean(mean(Evg(find(land_mask>0))));%*mm_fromkg_per_m2;


% Global plot of ET
var_avg = E_tot_mean_annual*mm_fromkg_per_m2; %mm/yr;
myscale=[1e0 1e3];
figure(1)
var_avg(find(land_mask<=0))=NaN;
TsK = rot90(var_avg);
TsK1 = TsK;
TsK1(1:jm,im+1) = TsK(1:jm,im);
TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
m_proj('Equidistant Cylindrical','lon',0,'lat',90);
%m_proj('Miller Cylindrical','lon',0,'lat',90);
%m_proj('mollweide','lon',[-180 180],'lat',[-90 90]);
% set the figure background to white
set(gcf,'color','white')
m_pcolor(Tlg,Tlt,TsK1);
%m_pcolor(Tlg,Tlt,log10(TsK1));
shading flat; % gets rid of the grid box edges
colormap(jet);%colormap(flipud(HSV));
m_coast('color',[0.3 0.3 0.3]);
m_grid('box','fancy','color',[0 0 0]);
%caxis(log10(myscale))
%caxis([0 1000])
h=colorbar('EastOutside');
%ticks_wanted=unique([myscale(1),get(h,'YTick'),myscale(2)]);
%ticks_wanted=10.^(0:3);
%caxis(log10(myscale))
%set(h,'YTick',log10(ticks_wanted));
%set(h,'YTickLabel',ticks_wanted);
set(get(h,'title'),'string','mm/yr');
title('Mean evapotranspiration from bare & vegetated soil for 1983')
text(0,-1,['Land mean = ',num2str(total_Etot,'%10.2f'),' mm/yr'])

% Global plot of Ecan
var_avg = E_can_mean_annual;%./E_tot_mean_annual; 
myscale=[1e0 1e3];
figure(2)
var_avg(find(land_mask<=0))=NaN;
TsK = rot90(var_avg);
TsK1 = TsK;
TsK1(1:jm,im+1) = TsK(1:jm,im);
TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
m_proj('Equidistant Cylindrical','lon',0,'lat',90);
%m_proj('Miller Cylindrical','lon',0,'lat',90);
%m_proj('mollweide','lon',[-180 180],'lat',[-90 90]);
% set the figure background to white
set(gcf,'color','white')
m_pcolor(Tlg,Tlt,TsK1);
%m_pcolor(Tlg,Tlt,log10(TsK1));
shading flat; % gets rid of the grid box edges
colormap(jet);%colormap(flipud(HSV));
m_coast('color',[0.3 0.3 0.3]);
m_grid('box','fancy','color',[0 0 0]);
%caxis(log10(myscale))
%caxis([0 1])
h=colorbar('EastOutside');
%ticks_wanted=unique([myscale(1),get(h,'YTick'),myscale(2)]);
%ticks_wanted=10.^(0:3);
%caxis(log10(myscale))
%set(h,'YTick',log10(ticks_wanted));
%set(h,'YTickLabel',ticks_wanted);
%set(get(h,'title'),'string','mm/yr');
title('Canopy Evaporation for 1983');%as fraction of total ET for 1983')
text(0,-0.85,['Land mean = ',num2str(total_Ecan,'%10.2f'),' mm/yr'])
text(0,-1.1,['f_{Ecan} = ',num2str(E_frac_can,'%10.2f')])

% Global plot of Etran
var_avg = E_tran_mean_annual;%./E_tot_mean_annual; 
myscale=[1e0 1e3];
figure(3)
var_avg(find(land_mask<=0))=NaN;
TsK = rot90(var_avg);
TsK1 = TsK;
TsK1(1:jm,im+1) = TsK(1:jm,im);
TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
m_proj('Equidistant Cylindrical','lon',0,'lat',90);
%m_proj('Miller Cylindrical','lon',0,'lat',90);
%m_proj('mollweide','lon',[-180 180],'lat',[-90 90]);
% set the figure background to white
set(gcf,'color','white')
m_pcolor(Tlg,Tlt,TsK1);
%m_pcolor(Tlg,Tlt,log10(TsK1));
shading flat; % gets rid of the grid box edges
colormap(jet);%colormap(flipud(HSV));
m_coast('color',[0.3 0.3 0.3]);
m_grid('box','fancy','color',[0 0 0]);
%caxis(log10(myscale))
%caxis([0 1])
h=colorbar('EastOutside');
%ticks_wanted=unique([myscale(1),get(h,'YTick'),myscale(2)]);
%ticks_wanted=10.^(0:3);
%caxis(log10(myscale))
%set(h,'YTick',log10(ticks_wanted));
%set(h,'YTickLabel',ticks_wanted);
%set(get(h,'title'),'string','mm/yr');
title('Transpiration for 1983');%as fraction of total ET for 1983')
text(0,-0.85,[' Land mean  = ',num2str(total_Etran,'%10.2f'),' mm/yr'])
text(0,-1.1,['f_{Etran} = ',num2str(E_frac_tran,'%10.2f')])

% Global plot of Ebare
var_avg = E_bare_mean_annual;%./E_tot_mean_annual; 
myscale=[1e0 1e3];
figure(4)
var_avg(find(land_mask<=0))=NaN;
TsK = rot90(var_avg);
TsK1 = TsK;
TsK1(1:jm,im+1) = TsK(1:jm,im);
TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
m_proj('Equidistant Cylindrical','lon',0,'lat',90);
%m_proj('Miller Cylindrical','lon',0,'lat',90);
%m_proj('mollweide','lon',[-180 180],'lat',[-90 90]);
% set the figure background to white
set(gcf,'color','white')
m_pcolor(Tlg,Tlt,TsK1);
%m_pcolor(Tlg,Tlt,log10(TsK1));
shading flat; % gets rid of the grid box edges
colormap(jet);%colormap(flipud(HSV));
m_coast('color',[0.3 0.3 0.3]);
m_grid('box','fancy','color',[0 0 0]);
%caxis(log10(myscale))
%caxis([0 1])
h=colorbar('EastOutside');
%ticks_wanted=unique([myscale(1),get(h,'YTick'),myscale(2)]);
%ticks_wanted=10.^(0:3);
%caxis(log10(myscale))
%set(h,'YTick',log10(ticks_wanted));
%set(h,'YTickLabel',ticks_wanted);
%set(get(h,'title'),'string','mm/yr');
title('Bare soil evaporation for 1983');%as fraction of total ET for 1983')
text(0,-0.85,['Land mean = ',num2str(total_Ebare,'%10.2f'),' mm/yr'])
%text(0,-1.1,['f_{Ebare} = ',num2str(E_frac_bare,'%10.2f')])

% Global plot of Eevapveg
var_avg = E_evapveg_mean_annual;%./E_tot_mean_annual; 
myscale=[1e0 1e3];
figure(14)
var_avg(find(land_mask<=0))=NaN;
TsK = rot90(var_avg);
TsK1 = TsK;
TsK1(1:jm,im+1) = TsK(1:jm,im);
TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
m_proj('Equidistant Cylindrical','lon',0,'lat',90);
%m_proj('Miller Cylindrical','lon',0,'lat',90);
%m_proj('mollweide','lon',[-180 180],'lat',[-90 90]);
% set the figure background to white
set(gcf,'color','white')
m_pcolor(Tlg,Tlt,TsK1);
%m_pcolor(Tlg,Tlt,log10(TsK1));
shading flat; % gets rid of the grid box edges
colormap(jet);%colormap(flipud(HSV));
m_coast('color',[0.3 0.3 0.3]);
m_grid('box','fancy','color',[0 0 0]);
%caxis(log10(myscale))
%caxis([0 1])
h=colorbar('EastOutside');
%ticks_wanted=unique([myscale(1),get(h,'YTick'),myscale(2)]);
%ticks_wanted=10.^(0:3);
%caxis(log10(myscale))
%set(h,'YTick',log10(ticks_wanted));
%set(h,'YTickLabel',ticks_wanted);
%set(get(h,'title'),'string','mm/yr');
title('Evaporation from vegetated soil for 1983');% as fraction of total ET for 1983')
text(0,-0.85,['Land mean = ',num2str(total_Evg,'%10.2f'),' mm/yr'])
text(0,-1.1,['f_{Evg} = ',num2str(E_frac_evapveg,'%10.2f')])


% Global plot of Epot
var_avg = E_pot_mean_annual*mm_fromkg_per_m2; %mm/yr;
myscale=[1e0 1e3];
figure(5)
var_avg(find(land_mask<=0))=NaN;
TsK = rot90(var_avg);
TsK1 = TsK;
TsK1(1:jm,im+1) = TsK(1:jm,im);
TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
m_proj('Equidistant Cylindrical','lon',0,'lat',90);
%m_proj('Miller Cylindrical','lon',0,'lat',90);
%m_proj('mollweide','lon',[-180 180],'lat',[-90 90]);
% set the figure background to white
set(gcf,'color','white')
m_pcolor(Tlg,Tlt,TsK1);
%m_pcolor(Tlg,Tlt,log10(TsK1));
shading flat; % gets rid of the grid box edges
colormap(jet);%colormap(flipud(HSV));
m_coast('color',[0.3 0.3 0.3]);
m_grid('box','fancy','color',[0 0 0]);
%caxis(log10(myscale))
%caxis([0 1200])
h=colorbar('EastOutside');
%ticks_wanted=unique([myscale(1),get(h,'YTick'),myscale(2)]);
%ticks_wanted=10.^(0:3);
%caxis(log10(myscale))
%set(h,'YTick',log10(ticks_wanted));
%set(h,'YTickLabel',ticks_wanted);
set(get(h,'title'),'string','mm/yr');
title('Potential evaporation for 1983')
text(0,-1,['Land mean = ',num2str(total_Epot,'%10.2f'),' mm/yr'])


% Global plot of total runoff
var_avg = (Rs_mean_annual+Ru_mean_annual)*mm_fromkg_per_m2;
myscale=[1e0 1e3];
figure(7)
var_avg(find(land_mask<=0))=NaN;
TsK = rot90(var_avg);
TsK1 = TsK;
TsK1(1:jm,im+1) = TsK(1:jm,im);
TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
m_proj('Equidistant Cylindrical','lon',0,'lat',90);
%m_proj('Miller Cylindrical','lon',0,'lat',90);
%m_proj('mollweide','lon',[-180 180],'lat',[-90 90]);
% set the figure background to white
set(gcf,'color','white')
m_pcolor(Tlg,Tlt,log10(TsK1));
shading flat; % gets rid of the grid box edges
colormap(jet);%colormap(flipud(HSV));
m_coast('color',[0.3 0.3 0.3]);
m_grid('box','fancy','color',[0 0 0]);
caxis(log10(myscale))
%caxis([1000 1000000])%10409294
h=colorbar('EastOutside');
%ticks_wanted=unique([myscale(1),get(h,'YTick'),myscale(2)]);
ticks_wanted=10.^(0:3);
caxis(log10(myscale))
set(h,'YTick',log10(ticks_wanted));
set(h,'YTickLabel',ticks_wanted);
set(get(h,'title'),'string','mm/month');
title('Total GISS LSM runoff (original scheme) for 1983')
total_Rsu=total_Rs+total_Ru;
text(0,-0.85,['Land mean = ',num2str(total_Rsu,'%10.2f'),' mm/yr'])
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Time series plots
E_vsevap_time=E_tot_time-E_can_time-E_tran_time-E_bare_time;

figure(8)
plot(tlsm_mnth,E_tot_time,'r-*',tlsm_mnth,E_can_time,'k-<',...
     tlsm_mnth,E_tran_time,'g-^',tlsm_mnth,E_bare_time,'c->',...
     tlsm_mnth,E_vsevap_time,'r--');
axis([1986 1996 0 10]);
xr = xlabel('t [yr]');
yr = ylabel('Water Flux [kg/m^2/month]');
%legend('potential','total','canopy','transpiration','bare soil evaporation');
legend('total','canopy ev.','transpir.','bare soil ev','veg soil ev');
title('Evapotranspiration partitioning for 1986-1995')
%--------------------------------------------------------------------------

% Global plot of GPP
figure(9)
var_avg = GPP_mean_annual;
var_avg(find(land_mask<=0))=NaN;
TsK = rot90(var_avg);
TsK1 = TsK;
TsK1(1:jm,im+1) = TsK(1:jm,im);
TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
m_proj('Equidistant Cylindrical','lon',0,'lat',90);
%m_proj('Miller Cylindrical','lon',0,'lat',90);
%m_proj('mollweide','lon',[-180 180],'lat',[-90 90]);
% set the figure background to white
set(gcf,'color','white')
m_pcolor(Tlg,Tlt,TsK1);
%m_pcolor(Tlg,Tlt,log10(TsK1));
shading flat; % gets rid of the grid box edges
colormap(jet);%colormap(flipud(HSV));
m_coast('color',[0.3 0.3 0.3]);
m_grid('box','fancy','color',[0 0 0]);
%caxis(log10(myscale))
%caxis([0 3])
h=colorbar('EastOutside');
%ticks_wanted=unique([myscale(1),get(h,'YTick'),myscale(2)]);
%ticks_wanted=10.^(0:3);
%caxis(log10(myscale))
%set(h,'YTick',log10(ticks_wanted));
%set(h,'YTickLabel',ticks_wanted);
set(get(h,'title'),'string','kg C/m^2/yr');
title('GPP for 1983')
text(0,-0.95,['\Sigma GPP = ',num2str(total_gpp,'%10.2f'),' Pg C/yr'])
%--------------------------------------------------------------------------

% Global plot of NPP
figure(10)
var_avg = NPP_mean_annual;
var_avg(find(land_mask<=0))=NaN;
TsK = rot90(var_avg);
TsK1 = TsK;
TsK1(1:jm,im+1) = TsK(1:jm,im);
TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
m_proj('Equidistant Cylindrical','lon',0,'lat',90);
%m_proj('Miller Cylindrical','lon',0,'lat',90);
%m_proj('mollweide','lon',[-180 180],'lat',[-90 90]);
% set the figure background to white
set(gcf,'color','white')
m_pcolor(Tlg,Tlt,TsK1);
%m_pcolor(Tlg,Tlt,log10(TsK1));
shading flat; % gets rid of the grid box edges
colormap(jet);%colormap(flipud(HSV));
m_coast('color',[0.3 0.3 0.3]);
m_grid('box','fancy','color',[0 0 0]);
%caxis(log10(myscale))
%caxis([0 3])
h=colorbar('EastOutside');
%ticks_wanted=unique([myscale(1),get(h,'YTick'),myscale(2)]);
%ticks_wanted=10.^(0:3);
%caxis(log10(myscale))
%set(h,'YTick',log10(ticks_wanted));
%set(h,'YTickLabel',ticks_wanted);
set(get(h,'title'),'string','kg C/m^2/yr');
title('NPP for 1983')
text(0,-0.95,['\Sigma NPP = ',num2str(total_NPP,'%10.2f'),' Pg C/yr'])
%--------------------------------------------------------------------------

% Global plot of NEE
figure(11)
var_avg = NEE_mean_annual;
var_avg(find(land_mask<=0))=NaN;
TsK = rot90(var_avg);
TsK1 = TsK;
TsK1(1:jm,im+1) = TsK(1:jm,im);
TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
m_proj('Equidistant Cylindrical','lon',0,'lat',90);
%m_proj('Miller Cylindrical','lon',0,'lat',90);
%m_proj('mollweide','lon',[-180 180],'lat',[-90 90]);
% set the figure background to white
set(gcf,'color','white')
m_pcolor(Tlg,Tlt,TsK1);
%m_pcolor(Tlg,Tlt,log10(TsK1));
shading flat; % gets rid of the grid box edges
colormap(jet);%colormap(flipud(HSV));
m_coast('color',[0.3 0.3 0.3]);
m_grid('box','fancy','color',[0 0 0]);
%caxis(log10(myscale))
%caxis([0 3])
h=colorbar('EastOutside');
%ticks_wanted=unique([myscale(1),get(h,'YTick'),myscale(2)]);
%ticks_wanted=10.^(0:3);
%caxis(log10(myscale))
%set(h,'YTick',log10(ticks_wanted));
%set(h,'YTickLabel',ticks_wanted);
set(get(h,'title'),'string','kg C/m^2/yr');
title('NEE for 1983')
text(0,-0.95,['\Sigma NEE = ',num2str(total_NEE,'%10.2f'),' Pg C/yr'])

%--------------------------------------------------------------------------

% Global plots of bare soil moisture
figure(41)
var_avg = SW_bare_mean_annual;
var_avg(find(land_mask<=0))=NaN;
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
%caxis([0 3])
h=colorbar('EastOutside');
title('Soil moisture of the top 10 cm of bare soil for 1983')
text(0,-0.95,['<S_{b1}> = ',num2str(total_SW_bare,'%10.2f')])

figure(42)
var_avg = SW2_bare_mean_annual;
var_avg(find(land_mask<=0))=NaN;
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
%caxis([0 3])
h=colorbar('EastOutside');
title('Soil moisture of bare layer between 10 cm and 27 cm for 1983')
text(0,-0.95,['<S_{b2}> = ',num2str(total_SW2_bare,'%10.2f')])

figure(43)
var_avg = SW3_bare_mean_annual;
var_avg(find(land_mask<=0))=NaN;
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
%caxis([0 3])
h=colorbar('EastOutside');
title('Soil moisture of bare layer between 27 cm and 57 cm for 1983')
text(0,-0.95,['<S_{b3}> = ',num2str(total_SW3_bare,'%10.2f')])

figure(44)
var_avg = SW4_bare_mean_annual;
var_avg(find(land_mask<=0))=NaN;
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
%caxis([0 3])
h=colorbar('EastOutside');
title('Soil moisture of bare layer between 57 cm and 1.1 m for 1983')
text(0,-0.95,['<S_{b4}> = ',num2str(total_SW4_bare,'%10.2f')])

figure(45)
var_avg = SW5_bare_mean_annual;
var_avg(find(land_mask<=0))=NaN;
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
%caxis([0 3])
h=colorbar('EastOutside');
title('Soil moisture of bare layer between 1.1 m and 2 m for 1983')
text(0,-0.95,['<S_{b5}> = ',num2str(total_SW5_bare,'%10.2f')])

figure(46)
var_avg = SW6_bare_mean_annual;
var_avg(find(land_mask<=0))=NaN;
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
%caxis([0 3])
h=colorbar('EastOutside');
title('Soil moisture of bare layer between 2 m and 3.5 m for 1983')
text(0,-0.95,['<S_{b6}> = ',num2str(total_SW6_bare,'%10.2f')])
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
    temp(find(fearth==0))=NaN;
    VEG(:,:,i) = temp;
    skip = fread(fid,[4]);
end
fclose(fid);

%--------------------------------------------------------------------------
% Dryness Index Computations

% % Plot annual dryness index versus annual GPP
% figure(9)
% subplot(3,2,1)
% decid=VEG(:,:,6);
% temp_GPP=GPP_mean_annual.*(decid>1e-2);
% temp_DI=annual_dryness.*(decid>1e-2);
% temp_DI(find(land_mask<=0))=0;
% temp_GPP(find(land_mask<=0))=0;
% 
% GPP_vec1 = reshape(temp_GPP,[],1);
% DI_vec1 = reshape(temp_DI,[],1);
% scatter(DI_vec1,GPP_vec1)
% xr = xlabel('Dryness Index, D_I');
% yr = ylabel('GPP_{cell} [kg C/m^2/month]');
% title('Annual GPP versus D_I for deciduous forests')
% 
% subplot(3,2,2)
% bin_width = 0.2;
% bins =  (0:bin_width:20)';
% [n,bin] = histc(DI_vec1,bins);
% num_count = sum(n);
% stairs(bins,n/num_count/bin_width)
% xr = xlabel('Dryness Index, D_I');
% yr = ylabel('# in bin/total #/ bin width');
% title('Normalized Histogram of D_I for deciduous forests')
% axis([0 10 0 5]);
% 
% subplot(3,2,3)
% grass=VEG(:,:,3);
% temp_GPP=GPP_mean_annual.*(grass>1e-2);
% temp_DI=annual_dryness.*(grass>1e-2);
% temp_DI(find(land_mask<=0))=0;
% temp_GPP(find(land_mask<=0))=0;
% GPP_vec3 = reshape(temp_GPP,[],1);
% DI_vec3 = reshape(temp_DI,[],1);
% 
% 
% scatter(DI_vec3,GPP_vec3)
% xr = xlabel('Dryness Index, D_I');
% yr = ylabel('GPP_{cell} [kg C/m^2/month]');
% title('Annual GPP as a function of D_I for grasslands')
% 
% subplot(3,2,4)
% bin_width = 0.2;
% bins =  (0:bin_width:20)';
% [n,bin] = histc(DI_vec3,bins);
% num_count = sum(n);
% stairs(bins,n/num_count/bin_width)
% xr = xlabel('Dryness Index, D_I');
% yr = ylabel('# in bin/total #/ bin width');
% title('Normalized Histogram of D_I for grasslands')
% axis([0 10 0 5]);
% 
% subplot(3,2,5)
% evgreen=VEG(:,:,7);
% temp_GPP=GPP_mean_annual.*(evgreen>1e-2);
% temp_DI=annual_dryness.*(evgreen>1e-2);
% temp_DI(find(land_mask<=0))=0;
% temp_GPP(find(land_mask<=0))=0;
% 
% GPP_vec2 = reshape(temp_GPP,[],1);
% DI_vec2 = reshape(temp_DI,[],1);
% scatter(DI_vec2,GPP_vec2)
% xr = xlabel('Dryness Index, D_I');
% yr = ylabel('GPP_{cell} [kg C/m^2/month]');
% title('Annual GPP versus of D_I for rainforests')
% 
% subplot(3,2,6)
% bin_width = 0.2;
% bins =  (0:bin_width:20)';
% [n,bin] = histc(DI_vec2,bins);
% num_count = sum(n);
% stairs(bins,n/num_count/bin_width)
% xr = xlabel('Dryness Index, D_I');
% yr = ylabel('# in bin/total #/ bin width');
% title('Normalized Histogram of D_I for rainforests')
% axis([0 10 0 5]);
% 
% % Plot seasonal dryness index versus seasonal GPP
% figure(10)
% subplot(3,2,1)
% decid=VEG(:,:,6);
% temp_GPP=seasonal_GPP_avg(:,:,3).*(decid>1e-2).*N_hemisph;
% temp_DI=seasonal_dryness(:,:,3).*(decid>1e-2).*N_hemisph;
% GPP_vec = reshape(temp_GPP,[],1);
% DI_vec = reshape(temp_DI,[],1);
% scatter(DI_vec,GPP_vec)
% xr = xlabel('Dryness Index, D_I');
% yr = ylabel('GPP [kg C/m^2/month]');
% title('JJA GPP as a function of dryness index for deciduous forests')
% 
% subplot(3,2,2)
% bin_width = 0.1;
% bins =  (0:bin_width:20)';
% [n,bin] = histc(DI_vec,bins);
% num_count = sum(n);
% stairs(bins,n/num_count/bin_width)
% xr = xlabel('Dryness Index, D_I');
% yr = ylabel('# in bin/total #/ bin width');
% title('Normalized Histogram of D_I for deciduous forests')
% axis([0 2 0 8]);
% 
% subplot(3,2,3)
% grass=VEG(:,:,3);
% temp_GPP=seasonal_GPP_avg(:,:,3).*(grass>1e-2).*N_hemisph;
% temp_DI=seasonal_dryness(:,:,3).*(grass>1e-2).*N_hemisph;
% GPP_vec = reshape(temp_GPP,[],1);
% DI_vec = reshape(temp_DI,[],1);
% scatter(DI_vec,GPP_vec)
% xr = xlabel('Dryness Index, D_I');
% yr = ylabel('GPP [kg C/m^2/month]');
% title('JJA GPP as a function of dryness index for grasslands')
% 
% subplot(3,2,4)
% bin_width = 0.1;
% bins =  (0:bin_width:20)';
% [n,bin] = histc(DI_vec,bins);
% num_count = sum(n);
% stairs(bins,n/num_count/bin_width)
% xr = xlabel('Dryness Index, D_I');
% yr = ylabel('# in bin/total #/ bin width');
% title('Normalized Histogram of D_I for grasslands')
% axis([0 2 0 8]);
% 
% subplot(3,2,5)
% rainfor=VEG(:,:,8);
% temp_GPP=seasonal_GPP_avg(:,:,3).*(rainfor>1e-2).*N_hemisph;
% temp_DI=seasonal_dryness(:,:,3).*(rainfor>1e-2).*N_hemisph;
% GPP_vec = reshape(temp_GPP,[],1);
% DI_vec = reshape(temp_DI,[],1);
% scatter(DI_vec,GPP_vec)
% xr = xlabel('Dryness Index, D_I');
% yr = ylabel('GPP [kg C/m^2/month]');
% title('JJA GPP as a function of dryness index for rainforest')
% 
% subplot(3,2,6)
% bin_width = 0.1;
% bins =  (0:bin_width:20)';
% [n,bin] = histc(DI_vec,bins);
% num_count = sum(n);
% stairs(bins,n/num_count/bin_width)
% xr = xlabel('Dryness Index, D_I');
% yr = ylabel('# in bin/total #/ bin width');
% title('Normalized Histogram of D_I for rainforests')
% axis([0 2 0 8]);
% 
% % Plot seasonal dryness index versus seasonal GPP
% figure(11)
% subplot(3,1,1)
% decid=VEG(:,:,6);
% temp_GPP=seasonal_GPP_avg(:,:,1).*(decid>1e-2).*S_hemisph;
% temp_DI=seasonal_dryness(:,:,1).*(decid>1e-2).*S_hemisph;
% GPP_vec = reshape(temp_GPP,[],1);
% DI_vec = reshape(temp_DI,[],1);
% scatter(DI_vec,GPP_vec)
% xr = xlabel('Dryness Index');
% yr = ylabel('GPP [kg C/m^2/month]');
% title('DJF GPP as a function of dryness index for deciduous forests')
% 
% subplot(3,1,2)
% grass=VEG(:,:,3);
% temp_GPP=seasonal_GPP_avg(:,:,1).*(grass>1e-2).*S_hemisph;
% temp_DI=seasonal_dryness(:,:,1).*(grass>1e-2).*S_hemisph;
% GPP_vec = reshape(temp_GPP,[],1);
% DI_vec = reshape(temp_DI,[],1);
% scatter(DI_vec,GPP_vec)
% xr = xlabel('Dryness Index');
% yr = ylabel('GPP [kg C/m^2/month]');
% title('DJF GPP as a function of dryness index for grasslands')
% 
% subplot(3,1,3)
% rainfor=VEG(:,:,8);
% temp_GPP=seasonal_GPP_avg(:,:,1).*(rainfor>1e-2).*S_hemisph;
% temp_DI=seasonal_dryness(:,:,1).*(rainfor>1e-2).*S_hemisph;
% GPP_vec = reshape(temp_GPP,[],1);
% DI_vec = reshape(temp_DI,[],1);
% scatter(DI_vec,GPP_vec)
% xr = xlabel('Dryness Index');
% yr = ylabel('GPP [kg C/m^2/month]');
% title('DJF GPP as a function of dryness index for rainforest')
% 
% % Plot annual dryness index versus annual Rsoil
% figure(12)
% subplot(3,1,1)
% decid=VEG(:,:,6);
% temp_R_soil=R_soil_mean_annual.*(decid>0.2);
% temp_DI=annual_dryness.*(decid>0.2);
% R_soil_vec1 = reshape(temp_R_soil,[],1);
% DI_vec1 = reshape(temp_DI,[],1);
% scatter(DI_vec1,R_soil_vec1)
% xr = xlabel('Dryness Index, D_I');
% yr = ylabel('Rsoil_{cell} [kg C/m^2/month]');
% title('Annual R_soil versus D_I for deciduous forests')
% 
% subplot(3,1,2)
% grass=VEG(:,:,3);
% temp_R_soil=R_soil_mean_annual.*(grass>0.2);
% temp_DI=annual_dryness.*(grass>0.2);
% R_soil_vec2 = reshape(temp_R_soil,[],1);
% DI_vec2 = reshape(temp_DI,[],1);
% scatter(DI_vec2,R_soil_vec2)
% xr = xlabel('Dryness Index, D_I');
% yr = ylabel('Rsoil_{cell} [kg C/m^2/month]');
% title('Annual Rsoil as a function of D_I for grasslands')
% 
% subplot(3,1,3)
% evgreen=VEG(:,:,7);
% temp_R_soil=R_soil_mean_annual.*(evgreen>0.2);
% temp_DI=annual_dryness.*(evgreen>0.2);
% R_soil_vec3 = reshape(temp_R_soil,[],1);
% DI_vec3 = reshape(temp_DI,[],1);
% scatter(DI_vec3,R_soil_vec3)
% xr = xlabel('Dryness Index, D_I');
% yr = ylabel('Rsoil_{cell} [kg C/m^2/month]');
% title('Annual Rsoil versus of D_I for evergreen')
% 
% 
% 
% figure(13)
% subplot(3,1,1)
% 
% decid=VEG(:,:,6);
% temp_E_frac_can=(E_can_mean_annual./E_tot_mean_annual).*(decid>0.2);
% temp_E_frac_can(find(land_mask<=0))=0;
% E_frac_can_vec1 = reshape(temp_E_frac_can,[],1);
% 
% temp_E_frac_tran=E_frac_tran.*(decid>0.02);
% temp_E_frac_tran=(E_tran_mean_annual./E_tot_mean_annual).*(decid>0.2);
% temp_E_frac_tran(find(land_mask<=0))=0;
% E_frac_tran_vec1 = reshape(temp_E_frac_tran,[],1);
% 
% temp_E_frac_evapveg=E_frac_evapveg.*(decid>0.02);
% temp_E_frac_evapveg=(E_evapveg_mean_annual./E_tot_mean_annual).*(decid>0.2);
% temp_E_frac_evapveg(find(land_mask<=0))=0;
% E_frac_evapveg_vec1 = reshape(temp_E_frac_evapveg,[],1);
% 
% evgreen=VEG(:,:,7);
% temp_E_frac_can=(E_can_mean_annual./E_tot_mean_annual).*(evgreen>0.2);
% temp_E_frac_can(find(land_mask<=0))=0;
% E_frac_can_vec2 = reshape(temp_E_frac_can,[],1);
% 
% temp_E_frac_tran=(E_tran_mean_annual./E_tot_mean_annual).*(evgreen>0.2);
% temp_E_frac_tran(find(land_mask<=0))=0;
% E_frac_tran_vec2 = reshape(temp_E_frac_tran,[],1);
% 
% temp_E_frac_evapveg=(E_evapveg_mean_annual./E_tot_mean_annual).*(evgreen>0.2);
% temp_E_frac_evapveg(find(land_mask<=0))=0;
% E_frac_evapveg_vec2 = reshape(temp_E_frac_evapveg,[],1);
% 
% 
% grass=VEG(:,:,3);
% temp_E_frac_can=(E_can_mean_annual./E_tot_mean_annual).*(grass>0.2);
% temp_E_frac_can(find(land_mask<=0))=0;
% E_frac_can_vec3 = reshape(temp_E_frac_can,[],1);
% 
% temp_E_frac_tran=(E_tran_mean_annual./E_tot_mean_annual).*(grass>0.2);
% temp_E_frac_tran(find(land_mask<=0))=0;
% E_frac_tran_vec3 = reshape(temp_E_frac_tran,[],1);
% 
% temp_E_frac_evapveg=(E_evapveg_mean_annual./E_tot_mean_annual).*(grass>0.2);
% temp_E_frac_evapveg(find(land_mask<=0))=0;
% E_frac_evapveg_vec3 = reshape(temp_E_frac_evapveg,[],1);
% 
% 
% bin_width = 0.05;
% bins =  (0:bin_width:0.5)';
% hold on
% [n,bin] = histc(E_frac_can_vec,bins);
% num_count = sum(n);
% [nT,binT] = histc(E_frac_tran_vec,bins);
% num_countT = sum(nT);
% [nG,binG] = histc(E_frac_evapveg_vec,bins);
% num_countG = sum(nG);
% stairs(bins,n/num_count/bin_width,'g')
% stairs(bins,nT/num_countT/bin_width,'k')
% stairs(bins,nG/num_countG/bin_width,'b')
% hold off
% 
% xr = xlabel('ET component fraction');
% yr = ylabel('# in bin/total #/ bin width');
% title('Normalized Histogram of ET components')
% %legend('CanEvap','Transpir', 'SoilEvap')
% 
% 
% % Print global maps of vegetation fractions
% print_vegfrac=1;
% if (print_vegfrac == 1)
%     var_avg=VEG(:,:,1);
%     var_avg(find(var_avg<=0.2))=NaN;
%     figure(21)
%     TsK = rot90(var_avg);
%     TsK1 = TsK;
%     TsK1(1:jm,im+1) = TsK(1:jm,im);
%     TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
%     m_proj('Equidistant Cylindrical','lon',0,'lat',90);
%     % set the figure background to white
%     set(gcf,'color','white')
%     m_pcolor(Tlg,Tlt,TsK1);
%     shading flat; % gets rid of the grid box edges
%     colormap(jet);%colormap(flipud(HSV));
%     m_coast('color',[0.3 0.3 0.3]);
%     m_grid('box','fancy','color',[0 0 0]);
%     caxis([0 1])
%     h=colorbar('EastOutside');
%     title('Bright Bare Soil Fraction')
% 
%     var_avg=VEG(:,:,10);
%     var_avg(find(var_avg<=0.2))=NaN;
%     figure(22)
%     TsK = rot90(var_avg);
%     TsK1 = TsK;
%     TsK1(1:jm,im+1) = TsK(1:jm,im);
%     TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
%     m_proj('Equidistant Cylindrical','lon',0,'lat',90);
%     % set the figure background to white
%     set(gcf,'color','white')
%     m_pcolor(Tlg,Tlt,TsK1);
%     shading flat; % gets rid of the grid box edges
%     colormap(jet);%colormap(flipud(HSV));
%     m_coast('color',[0.3 0.3 0.3]);
%     m_grid('box','fancy','color',[0 0 0]);
%     caxis([0 1])
%     h=colorbar('EastOutside');
%     title('Dark Bare Soil Fraction')
% 
%     var_avg=VEG(:,:,2);
%     var_avg(find(var_avg<=0.2))=NaN;
%     figure(23)
%     TsK = rot90(var_avg);
%     TsK1 = TsK;
%     TsK1(1:jm,im+1) = TsK(1:jm,im);
%     TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
%     m_proj('Equidistant Cylindrical','lon',0,'lat',90);
%     % set the figure background to white
%     set(gcf,'color','white')
%     m_pcolor(Tlg,Tlt,TsK1);
%     shading flat; % gets rid of the grid box edges
%     colormap(jet);%colormap(flipud(HSV));
%     m_coast('color',[0.3 0.3 0.3]);
%     m_grid('box','fancy','color',[0 0 0]);
%     caxis([0 1])
%     h=colorbar('EastOutside');
%     title('Tundra Fraction')
% 
%     var_avg=VEG(:,:,3);
%     var_avg(find(var_avg<=0.2))=NaN;
%     figure(24)
%     TsK = rot90(var_avg);
%     TsK1 = TsK;
%     TsK1(1:jm,im+1) = TsK(1:jm,im);
%     TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
%     m_proj('Equidistant Cylindrical','lon',0,'lat',90);
%     % set the figure background to white
%     set(gcf,'color','white')
%     m_pcolor(Tlg,Tlt,TsK1);
%     shading flat; % gets rid of the grid box edges
%     colormap(jet);%colormap(flipud(HSV));
%     m_coast('color',[0.3 0.3 0.3]);
%     m_grid('box','fancy','color',[0 0 0]);
%     caxis([0 1])
%     h=colorbar('EastOutside');
%     title('Grassland Fraction')
% 
%     var_avg=VEG(:,:,4);
%     var_avg(find(var_avg<=0.2))=NaN;
%     figure(25)
%     TsK = rot90(var_avg);
%     TsK1 = TsK;
%     TsK1(1:jm,im+1) = TsK(1:jm,im);
%     TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
%     m_proj('Equidistant Cylindrical','lon',0,'lat',90);
%     % set the figure background to white
%     set(gcf,'color','white')
%     m_pcolor(Tlg,Tlt,TsK1);
%     shading flat; % gets rid of the grid box edges
%     colormap(jet);%colormap(flipud(HSV));
%     m_coast('color',[0.3 0.3 0.3]);
%     m_grid('box','fancy','color',[0 0 0]);
%     caxis([0 1])
%     h=colorbar('EastOutside');
%     title('Shrub/Grassland Fraction')
% 
%     var_avg=VEG(:,:,5);
%     var_avg(find(var_avg<=0.2))=NaN;
%     figure(26)
%     TsK = rot90(var_avg);
%     TsK1 = TsK;
%     TsK1(1:jm,im+1) = TsK(1:jm,im);
%     TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
%     m_proj('Equidistant Cylindrical','lon',0,'lat',90);
%     % set the figure background to white
%     set(gcf,'color','white')
%     m_pcolor(Tlg,Tlt,TsK1);
%     shading flat; % gets rid of the grid box edges
%     colormap(jet);%colormap(flipud(HSV));
%     m_coast('color',[0.3 0.3 0.3]);
%     m_grid('box','fancy','color',[0 0 0]);
%     caxis([0 1])
%     h=colorbar('EastOutside');
%     title('Tree/Grassland Fraction')
% 
%     var_avg=VEG(:,:,6);
%     var_avg(find(var_avg<=0.2))=NaN;
%     figure(27)
%     TsK = rot90(var_avg);
%     TsK1 = TsK;
%     TsK1(1:jm,im+1) = TsK(1:jm,im);
%     TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
%     m_proj('Equidistant Cylindrical','lon',0,'lat',90);
%     % set the figure background to white
%     set(gcf,'color','white')
%     m_pcolor(Tlg,Tlt,TsK1);
%     shading flat; % gets rid of the grid box edges
%     colormap(jet);%colormap(flipud(HSV));
%     m_coast('color',[0.3 0.3 0.3]);
%     m_grid('box','fancy','color',[0 0 0]);
%     caxis([0 1])
%     h=colorbar('EastOutside');
%     title('Deciduous Forest Fraction')
% 
%     var_avg=VEG(:,:,7);
%     var_avg(find(var_avg<=0.2))=NaN;
%     figure(28)
%     TsK = rot90(var_avg);
%     TsK1 = TsK;
%     TsK1(1:jm,im+1) = TsK(1:jm,im);
%     TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
%     m_proj('Equidistant Cylindrical','lon',0,'lat',90);
%     % set the figure background to white
%     set(gcf,'color','white')
%     m_pcolor(Tlg,Tlt,TsK1);
%     shading flat; % gets rid of the grid box edges
%     colormap(jet);%colormap(flipud(HSV));
%     m_coast('color',[0.3 0.3 0.3]);
%     m_grid('box','fancy','color',[0 0 0]);
%     caxis([0 1])
%     h=colorbar('EastOutside');
%     title('Evergreen Forest Fraction')
% 
%     var_avg=VEG(:,:,8);
%     var_avg(find(var_avg<=0.2))=NaN;
%     figure(29)
%     TsK = rot90(var_avg);
%     TsK1 = TsK;
%     TsK1(1:jm,im+1) = TsK(1:jm,im);
%     TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
%     m_proj('Equidistant Cylindrical','lon',0,'lat',90);
%     % set the figure background to white
%     set(gcf,'color','white')
%     m_pcolor(Tlg,Tlt,TsK1);
%     shading flat; % gets rid of the grid box edges
%     colormap(jet);%colormap(flipud(HSV));
%     m_coast('color',[0.3 0.3 0.3]);
%     m_grid('box','fancy','color',[0 0 0]);
%     caxis([0 1])
%     h=colorbar('EastOutside');
%     title('Rainforest Fraction')
% 
%     var_avg=VEG(:,:,9);
%     var_avg(find(var_avg<=0.2))=NaN;
%     figure(30)
%     TsK = rot90(var_avg);
%     TsK1 = TsK;
%     TsK1(1:jm,im+1) = TsK(1:jm,im);
%     TsK1(jm+1,1:im+1) = TsK1(jm,1:im+1);
%     m_proj('Equidistant Cylindrical','lon',0,'lat',90);
%     % set the figure background to white
%     set(gcf,'color','white')
%     m_pcolor(Tlg,Tlt,TsK1);
%     shading flat; % gets rid of the grid box edges
%     colormap(jet);%colormap(flipud(HSV));
%     m_coast('color',[0.3 0.3 0.3]);
%     m_grid('box','fancy','color',[0 0 0]);
%     caxis([0 1])
%     h=colorbar('EastOutside');
%     title('Cultivation Fraction')
% end