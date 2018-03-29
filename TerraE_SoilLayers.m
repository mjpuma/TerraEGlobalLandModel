%% Routine to extend soil layering using geometric progression
%   for the TerraE Global Land Model
clear all; close all; clc

%rootdir = '/Users/michaeljpuma/';
rootdir = '/Users/mpuma/';

% Open Current Soils file in ModelE
basefold = [rootdir 'ModelE_Feb42014/inputSGP/extractions/'];
filename = 'S144X900098M.ext.nc';
ncid = netcdf.open([basefold filename],'NC_NOWRITE');

% Get name and length of first dimension
[dimname, dimlen] = netcdf.inqDim(ncid,0);
dimIDs = netcdf.inqDimIDs(ncid);

% Open soil variables
dz_original_full = ncread([basefold filename],'dz');
q_original_full  = ncread([basefold filename],'q');
qk_original_full = ncread([basefold filename],'qk');


%% Single column values for original
dz_original = squeeze(dz_original_full);
zbottom_original = cumsum(dz_original);

z_original(1,1) = dz_original(1,1)/2;

depthsum = dz_original(1,1);
for i_layer = 2:6
    z_original(i_layer,1) = depthsum + dz_original(i_layer,1)/2;
    depthsum = depthsum+dz_original(i_layer,1);
end


% Original values in ModelE:
%     n, dz =  1,  9.99999642372131348E-2
%     n, dz =  2,  0.17254400253295898
%     n, dz =  3,  0.2977144718170166
%     n, dz =  4,  0.51368874311447144
%     n, dz =  5,  0.88633960485458374
%     n, dz =  6,  1.5293264389038086
common_ratio_original = dz_original(2)/dz_original(1);



%% Geometric progression for layering
% The n-th term of a geometric sequence with initial value a and common ratio r is given by
% a_n = a*r^{n-1}.
% Such a geometric sequence also follows the recursive relation
% a_n = r*a_{n-1} for every integer n\geq 1.

%%NOTES: 1) shallow configuration with 7 layers to a depth of ~3.676 m with
%%a first layer having a thickness of 5 cm.
%%2) deep configuration with 10 layers to a depth of ~21.43 m

common_ratio = 0;
a_layer1 = 0.05; % Thickness of first layer in meters
num_layers = 10;
zbottom_GISS=zeros(num_layers,1);

while zbottom_GISS(3,1) < .3
    common_ratio = common_ratio+.0001;
    dz_GISS  =zeros(num_layers,1);
    for i_layer = 1:num_layers
        dz_GISS(i_layer,1) = a_layer1*common_ratio^(i_layer-1);
    end
    total_depthGISS = sum(dz_GISS);
    
    z_GISS(1,1) = dz_GISS(1,1)/2;
    
    depthsum = dz_GISS(1,1);
    for i_layer = 2:num_layers
        z_GISS(i_layer,1) = depthsum + dz_GISS(i_layer,1)/2;
        depthsum = depthsum+dz_GISS(i_layer,1);
    end
    
    depthsum = 0;
    for i_layer = 1:num_layers
        depthsum = depthsum + dz_GISS(i_layer,1);
        zbottom_GISS(i_layer,1) = depthsum;
    end
    
end

%% FOR REFERENCE CLM4
%http://www.cesm.ucar.edu/models/cesm1.0/clm/CLM4_Tech_Note.pdf, p100
%exponential form is used to have more layers at the surface ... reasoning
%is assumed steep soil mositure gradient
% AU  - Lawrence, David M AU  - Oleson, Keith W TI  - Parameterization
% improvements and functional and structural advances in Version 4 of the
% Community Land Model JO  - Journal of Advances in Modeling Earth Systems
% JA  - J. Adv. Model. Earth Syst. VL  - 3 IS  - 1 SN  - 1942-2466 UR  -
% http://dx.doi.org/10.1029/2011MS00045
fs = 0.025;
num_layers_clm = 15;
for i_layer=1:num_layers_clm
    z_NCAR(i_layer,1)=fs*(exp(0.5*(i_layer-0.5))-1);
end

dz_NCAR(1,1) = (z_NCAR(1,1)+z_NCAR(2,1))/2;
for i_layer=2:num_layers_clm-1
    dz_NCAR(i_layer,1)=0.5*(z_NCAR(i_layer+1,1)+z_NCAR(i_layer-1,1));
end
dz_NCAR(num_layers_clm,1) = z_NCAR(num_layers_clm,1)-z_NCAR(num_layers_clm-1,1);

total_depthNCAR = sum(dz_NCAR);

%attvalue = ncreadatt([basefold filename],'/','ngm');

%% From CLM4.5
%http://nldr.library.ucar.edu/repository/collections/TECH-NOTE-000-000-000-870
% Compare to Table 6.1, p114

%%NOTES ON NCAR
%% 2.1.4. Soil/ground depth
% Nicolsky et al. (2007) and Alexeev et al. (2007) demonstrated
% that soil temperature dynamics cannot be accurately
% modeled with a shallow soil column and that a ground
% depth of at least 30 m is required for century-scale integrations.
% Therefore, in order to account for the thermal inertia
% of deep ground, the number of ground layers is extended in
% CLM4 from 10 to 15 layers, as in Lawrence et al. (2008).
% Layer thicknesses exponentially increase with depth, as
% before, ranging from a thickness of 0.018 m at the surface
% to 13.9 m for layer 15. The upper 10 layers are hydrologically
% active (i.e. the ?soil? layers) while the bottom five layers
% (3.8 m to 42 m depth) are thermal slabs that are not
% hydrologically active. The thermal conductivity for the deep
% ground layers is set at 3.0 W m21 K21, which is comparable
% to that reported for saturated granitic rock (Clauser and
% Huenges 1995), while the heat capacity is set to that of a
% generic rock (26106 J m23 K21) . The continued assumption
% of a globally uniform 3.8 m of hydrologically active soil
% remains unrealistic and is a deficiency of the model that
% requires attention in future development of the model.

%% Simplified bottom boundary condition for soil water
% equations
% In CLM3.5, the redistribution of water within the soil
% column/aquifer system takes place in two steps. In the first
% step, the soil hydrology equations are solved for the 10-layer
% soil column. Then, if the water table is deeper than the
% lowest soil layer, the aquifer recharge rate from the lowest
% soil layer to the unconfined aquifer is calculated. This twostep
% procedure decouples the water fluxes within the soil
% column from the flux of water between the lowest layer and
% the aquifer layer, leading on occasion to unrealistically large
% aquifer recharge rates.
% For CLM4, the aquifer is coupled directly to the soil
% column via the soil water equations, resulting in consistent
% moisture fluxes in the soil column / aquifer system. When
% the water table is within the soil column, a zero-flux
% boundary condition is applied at the bottom of the tenth
% layer, as in CLM3.5. When the water table drops from the
% lowest soil layer into the aquifer, an additional layer representing
% the portion of the aquifer between the bottom of the
% lowest layer and the water table is added to the system of soil
% water equations. The zero-flux boundary condition is then
% applied at the water table depth, rather than the bottom of
% the tenth layer.

%% Surface and subsurface runoff
% Surface runoff in CLM3.5 and CLM4 consists of overland
% flow due to saturation excess (Dunne runoff) and infiltration
% excess (Hortonian runoff) mechanisms. The saturation
% excess term is a function of the saturated fraction fsat of the
% soil column, which includes a dependence on the surface
% layer frozen soil impermeable area fraction ffrz,l (Niu and
% Yang 2006)
%      fsat~ 1{ffrz,1fmax expð{0:5fover z+Þzffrz,1 ð9Þ
% where fmax is the maximum saturated fraction, z+ is the
% water table depth, and fover is a decay factor. Subsurface
% runoff qdrai is calculated according to the following expression
% (Niu et al. 2005):
%      qdrai~ 1{fimp qdrai,max expð{fdraiz+Þ ð10Þ
% where fimp is the fraction of impermeable area determined
% from the ice content of the soil at depth, z+ is the water
% table depth, and fdrai is a decay factor. For CLM4, the decay
% factor fover and the maximum drainage qdrai,max when the
% water table is at the surface are adjusted through sensitivity
% analysis and comparison with observed runoff (fover5
% 2.5 in CLM3.5, fover50.5 in CLM4; qdrai,max 5 4.56
% 1024 kg m22 s21 in CLM3.5, qdrai,max 5 5.56
% 1023 kg m22 s21 in CLM4). The changes in these parameters
% help alleviate the wet soil bias detected in CLM3.5
% (Oleson et al. 2008c) and shifts the percentages of surface
% runoff and subsurface runoff from 30%:70% to 55%:45%.

