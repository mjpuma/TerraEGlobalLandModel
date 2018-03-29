% SINGLE CELL TAKEN FROM GLOBAL RUNS COUPLED LSM - ENT DGTEM (fort.9995)
clear all;  clc; close all;
%--------------------------------------------------------------------------
% Fluxnet sites for comparison/evaluation of model data
site_num = 2;

if (site_num ==1) % Vaira Ranch
    % Read in ecohydrologic data from Fluxnet sites for comparison/evaluation
    % of model data
    yr_start = 2001;
    tot_numyrs = 6; 
    num_leapyr = 0;
    num_spinup_leap = 1;
    num_spinup_nonleap = 4;
    plot_yr = 2006;
    plot_yr_num = 6;  
    num_day_plot_yr=365; 
    num_data_days=365;
    
    mnth_days = [31 28 31 30 31 30 31 31 30 31 30 31];
    tstep_per_day =48;    nee_steps = 48;
    
    i_yr = yr_start;
    for i=1:tot_numyrs
        numday=365;
        if( (i_yr/4-floor(i_yr/4)) == 0) 
            numday=366; 
        end
        num_yr(i) = tstep_per_day*numday;
        i_yr = i_yr+1;
    end
    
    tstep_to_day = tstep_per_day;
    totnum_days = (tot_numyrs-num_spinup_leap-num_spinup_nonleap)*365 + num_leapyr*366;
    tstep_to_yr = tstep_per_day*totnum_days;
    num_tsteps =totnum_days*tstep_per_day;  
    
    
    % Latent heat, sensible heat, and NEE data
    load /Users/mpuma/Documents/ET_partition_paper/Vaira2006_compdata.txt
%    LE_Wm2=Vaira2006_compdata(:,2);
%    H_Wm2=Vaira2006_compdata(:,3);
    NEE_umolm2s=Vaira2006_compdata(:,7); %Level 3, standarized, FC_umolm2s
    ncload  /NoBackup/FLUXNET_reto/selected/netcdf/ameriflux.Vaira.2006.nc
    LE_Wm2=LE;
    H_Wm2=H;
    %NEE_umolm2s=NEE;
    time=(0.5:num_yr(plot_yr_num)/tstep_per_day-.5)';
    tdata_day(1:num_day_plot_yr,1) = plot_yr+time/num_day_plot_yr;
    t_hr = ((1/nee_steps):(1/nee_steps):num_data_days)/num_day_plot_yr; 
    t_hr = plot_yr+t_hr';
    % Soil moisture data    
%     theta_10cm =Vaira2006_compdata(:,4);
%     theta_20cm =Vaira2006_compdata(:,5);
%     SWC30cm_Vaira = (theta_10cm+theta_20cm)/2;
    soilmoist1 = NaN;%SWC30cm_Vaira;
    soilmoist2 = NaN;%SWC30cm_Vaira;
    soilmoist3 = NaN;%SWC30cm_Vaira;
    t_sm = 1/nee_steps:(1/nee_steps):num_data_days*(1/nee_steps); 
    t_sm = t_sm';

    % Create time vector
    i_yr = yr_start;
    i_start = 0; i_stop = 0;
    for i=1:tot_numyrs
        if(i_start == 0)
            i_start = 1;
        else
            i_start = i_start+num_yr(i-1);
        end
        i_stop = i_stop+num_yr(i);
        
        numday=365;
        if( (i_yr/4-floor(i_yr/4)) == 0) 
            numday=366; 
        end
        
        time = 1/tstep_per_day:(1/tstep_per_day):numday; 
        time = time';      
        t_all(i_start:i_stop,1) = i_yr+time/(60*60*24)/numday;
        
        % advance year
        i_yr = i_yr+1;
    end   
    % Read in model output and assign only last year of data to vectors 
    load /NoBackup/LSM_single_out/Vaira_nosoilevap.9995;
    fort=Vaira_nosoilevap; clear Vaira_nosoilevap
    fort = fort(i_start:i_stop,:); 
    load /NoBackup/LSM_single_out/Vaira_soilevap_or.9995;
    fort2=Vaira_soilevap_or; clear Vaira_soilevap_or;
    fort2=fort2(i_start:i_stop,:);
    load /NoBackup/LSM_single_out/Vaira_soilevap_new.9995;
    fort3=Vaira_soilevap_new; clear Vaira_soilevap_new;
    fort3=fort3(i_start:i_stop,:);
    
    % Create time vector for the plotted year
    t = t_all(i_start:i_stop,:);

    % Daily time series
    i_start = 1;
    i_stop = num_yr(plot_yr_num)/tstep_per_day;
    time=(1.5:num_yr(plot_yr_num)/tstep_per_day+.5)';
    t_day(i_start:i_stop,1) = plot_yr+time/num_day_plot_yr;    
    
%--------------------------------------------------------------------------    
elseif(site_num ==2) % Morgan Monroe State Forest
    % Read in ecohydrologic data from Fluxnet sites for comparison/evaluation
    % of model data
    yr_start = 1983;
    tot_numyrs = 1; 
    num_leapyr = 0;
    num_spinup_leap = 0;
    num_spinup_nonleap = 0;
    plot_yr = 1983;
    plot_yr_num = 1;  
    num_day_plot_yr=365; 
    num_data_days=365;
    
    mnth_days = [31 28 31 30 31 30 31 31 30 31 30 31];
    tstep_per_day =48;    nee_steps = 24;
    
    i_yr = yr_start;
    for i=1:tot_numyrs
        numday=365;
        if( (i_yr/4-floor(i_yr/4)) == 0) 
            numday=366; 
        end
        num_yr(i) = tstep_per_day*numday;
        i_yr = i_yr+1;
    end
    
    tstep_to_day = tstep_per_day;
    totnum_days = (tot_numyrs-num_spinup_leap-num_spinup_nonleap)*365 + num_leapyr*366;
    tstep_to_yr = tstep_per_day*totnum_days;
    num_tsteps =totnum_days*tstep_per_day; 
    t_subday = ((1/tstep_per_day):(1/tstep_per_day):num_data_days)/num_day_plot_yr; 
    t_subday = plot_yr+t_subday';

%%%% Latent heat, sensible heat, and NEE data
%     load /Users/mpuma/Documents/ET_partition_paper/MMSF2005_compdata.txt
%     LE_Wm2=MMSF2005_compdata(:,2);
%     H_Wm2=MMSF2005_compdata(:,3);
%     NEE_umolm2s=MMSF2005_compdata(:,5);
    ncload('/NoBackup/FLUXNET_reto/selected/netcdf/ameriflux.Morgan_Monroe_State_Forest.2003.nc',...
        'LE','H','NEE','SoilTemp1','SoilMoist1','DoY','Hour','SWdown','time','Tair')
    t_fluxnet=time; clear time;
    LE_Wm2=LE;
    H_Wm2=H;
    NEE_umolm2s=NEE;
    time=(0.5:num_yr(plot_yr_num)/tstep_per_day-.5)';
    tdata_day(1:num_day_plot_yr,1) = plot_yr+time/num_day_plot_yr;
    t_hr = ((1/nee_steps):(1/nee_steps):num_data_days)/num_day_plot_yr; 
    t_hr = plot_yr+t_hr';
    
    % Soil moisture data
    load /Users/mpuma/Documents/Ent_standalone/MMSF/MMSF2005_soilmoist30cm.csv
    SWC30cm_MMSF = MMSF2005_soilmoist30cm(:,1);
    soilmoist1 = SWC30cm_MMSF;
    soilmoist2 = SWC30cm_MMSF;
    t_sm = 1/nee_steps:(1/nee_steps):num_day_plot_yr*(1/nee_steps); 
    t_sm = t_sm';
    
    % Create time vector
    i_yr = yr_start;
    i_start = 0; i_stop = 0;
    for i=1:tot_numyrs
        if(i_start == 0)
            i_start = 1;
        else
            i_start = i_start+num_yr(i-1);
        end
        i_stop = i_stop+num_yr(i);
        
        numday=365;
        if( (i_yr/4-floor(i_yr/4)) == 0) 
            numday=366; 
        end
        
        time = 1/tstep_per_day:(1/tstep_per_day):numday; 
        time = time';      
        t_all(i_start:i_stop,1) = i_yr+time/(60*60*24)/numday;
        
        % advance year
        i_yr = i_yr+1;
    end

    % Read in model output and assign only last year of data to vectors
    load /NoBackup/LSM_single_out/Global_MMSFnointerp.9995;
    fort=Global_MMSFnointerp; clear Global_MMSFnointerp
    fort = fort(i_start:i_stop,:); 
    load /NoBackup/LSM_single_out/Global_MMSFlininterp.9995;
    fort2=Global_MMSFlininterp; clear Global_MMSFlininterp;
    fort2=fort2(i_start:i_stop,:);
    load /NoBackup/LSM_single_out/Global_MMSF.9995;
    fort3=Global_MMSF; clear Global_MMSF;
    fort3=fort3(i_start:i_stop,:);
    
    % Create time vector for the plotted year
    t = t_all(i_start:i_stop,:);

    % Daily time series
    i_start = 1;
    i_stop = num_yr(plot_yr_num)/tstep_per_day;
    time=(1.5:num_yr(plot_yr_num)/tstep_per_day+.5)';
    t_day(i_start:i_stop,1) = plot_yr+time/num_day_plot_yr;
    
%--------------------------------------------------------------------------   
elseif(site_num ==3) % Hyytiala, Finland
    % Read in ecohydrologic data from Fluxnet sites for comparison/evaluation
    % of model data
    yr_start = 1997;
    tot_numyrs = 9; 
    num_leapyr = 0;
    num_spinup_leap = 2;
    num_spinup_nonleap = 6;
    plot_yr = 2005;
    plot_yr_num = 9;  
    num_day_plot_yr=365; 
    num_data_days=365;
    
    mnth_days = [31 28 31 30 31 30 31 31 30 31 30 31];
    tstep_per_day =48;    nee_steps = 48;
    
    i_yr = yr_start;
    for i=1:tot_numyrs
        numday=365;
        if( (i_yr/4-floor(i_yr/4)) == 0) 
            numday=366; 
        end
        num_yr(i) = tstep_per_day*numday;
        i_yr = i_yr+1;
    end
    
    tstep_to_day = tstep_per_day;
    totnum_days = (tot_numyrs-num_spinup_leap-num_spinup_nonleap)*365 + num_leapyr*366;
    tstep_to_yr = tstep_per_day*totnum_days;
    num_tsteps =totnum_days*tstep_per_day; 

    % Latent heat, sensible heat, and NEE data
%     load /Users/mpuma/Documents/Hyytiala_data/1788959_part02/CEIP_EC_L3_FIHyy_2005_notitle.txt
%     LE_Wm2=CEIP_EC_L3_FIHyy_2005_notitle(:,15);
%     H_Wm2=CEIP_EC_L3_FIHyy_2005_notitle(:,14);
%     NEE_umolm2s=CEIP_EC_L3_FIHyy_2005_notitle(:,10); %Level 3, standarized
    ncload  /NoBackup/FLUXNET_reto/selected/netcdf/carboeurope.Hyytiala.2005.nc
    LE_Wm2=LE;
    H_Wm2=H;
    NEE_umolm2s=FCO2;
    time=(0.5:num_yr(plot_yr_num)/tstep_per_day-.5)';
    tdata_day(1:num_day_plot_yr,1) = plot_yr+time/num_day_plot_yr;
    t_hr = ((1/nee_steps):(1/nee_steps):num_data_days)/num_day_plot_yr; 
    t_hr = plot_yr+t_hr';
    
    % Soil moisture data
    load /Users/mpuma/Documents/Hyytiala_data/1788959_part02/CEIP_EC_L3_FIHyy_2005_notitle.txt
    % FLUXNET site measurements (Hyyti?l?: ?4 to ?30 cm and ?30 to ?68 cm)
    % Grainer(spell???) 2007
    SWC1 = CEIP_EC_L3_FIHyy_2005_notitle(:,31)/100;
    SWC2 = CEIP_EC_L3_FIHyy_2005_notitle(:,32)/100;
    t_sm = 1/nee_steps:(1/nee_steps):num_data_days*(1/nee_steps); 
    t_sm = t_sm';
    soilmoist1 = SWC1;
    soilmoist2 = SWC1;
    soilmoist3 = SWC2;
    
    % Create time vector
    i_yr = yr_start;
    i_start = 0; i_stop = 0;
    for i=1:tot_numyrs
        if(i_start == 0)
            i_start = 1;
        else
            i_start = i_start+num_yr(i-1);
        end
        i_stop = i_stop+num_yr(i);
        
        numday=365;
        if( (i_yr/4-floor(i_yr/4)) == 0) 
            numday=366; 
        end
        
        time = 1/tstep_per_day:(1/tstep_per_day):numday; 
        time = time';      
        t_all(i_start:i_stop,1) = i_yr+time/(60*60*24)/numday;
        
        % advance year
        i_yr = i_yr+1;
    end   
    % Read in model output and assign only last year of data to vectors 
    load /NoBackup/LSM_single_out/Hyytiala_nosoilevap.9995;
    fort=Hyytiala_nosoilevap; clear Hyytiala_nosoilevap
    fort = fort(i_start:i_stop,:); 
    load /NoBackup/LSM_single_out/Hyytiala_soilevap_or.9995;
    fort2=Hyytiala_soilevap_or; clear Hyytiala_soilevap_or;
    fort2=fort2(i_start:i_stop,:);
    load /NoBackup/LSM_single_out/Hyytiala_soilevap_new.9995;
    fort3=Hyytiala_soilevap_new; clear Hyytiala_soilevap_new;
    fort3=fort3(i_start:i_stop,:);
    
    % Create time vector for the plotted year
    t = t_all(i_start:i_stop,:);

    % Daily time series
    i_start = 1;
    i_stop = num_yr(plot_yr_num)/tstep_per_day;
    time=(1.5:num_yr(plot_yr_num)/tstep_per_day+.5)';
    t_day(i_start:i_stop,1) = plot_yr+time/num_day_plot_yr;

%--------------------------------------------------------------------------
elseif(site_num == 4) % Santarem km 83, Brazil

    % Read in ecohydrologic data from Fluxnet sites for comparison/evaluation
    % of model data
    yr_start = 2001;
    tot_numyrs = 3; 
    num_leapyr = 0;
    num_spinup_leap = 0;
    num_spinup_nonleap = 2;
    plot_yr = 2003;
    plot_yr_num = 3;  
    num_day_plot_yr=365; 
    num_data_days=365;
    
    mnth_days = [31 28 31 30 31 30 31 31 30 31 30 31];
    tstep_per_day =24;    nee_steps = 24;
    
    i_yr = yr_start;
    for i=1:tot_numyrs
        numday=365;
        if( (i_yr/4-floor(i_yr/4)) == 0) 
            numday=366; 
        end
        num_yr(i) = tstep_per_day*numday;
        i_yr = i_yr+1;
    end
    
    tstep_to_day = tstep_per_day;
    totnum_days = (tot_numyrs-num_spinup_leap-num_spinup_nonleap)*365 + num_leapyr*366;
    tstep_to_yr = tstep_per_day*totnum_days;
    num_tsteps =totnum_days*tstep_per_day; 

    % Latent heat, sensible heat, and NEE data
    ncload /NoBackup/FLUXNET_reto/selected/netcdf/lba.Santarem_KM83.2003.nc
    LE_Wm2=LE;
    H_Wm2=H;
    NEE_umolm2s=NEE;  
    time=(0.5:num_yr(plot_yr_num)/tstep_per_day-.5)';
    tdata_day(1:num_day_plot_yr,1) = plot_yr+time/num_day_plot_yr;
    t_hr = ((1/nee_steps):(1/nee_steps):num_data_days)/num_day_plot_yr; 
    t_hr = plot_yr+t_hr';

    
    % Soil moisture data
    %load ????
    soilmoist1 = NaN;
    soilmoist2 = NaN;
    soilmoist3 = NaN;

    t_sm = 1/nee_steps:(1/nee_steps):num_data_days*(1/nee_steps); 
    t_sm = t_sm';
    
    % Create time vector
    i_yr = yr_start;
    i_start = 0; i_stop = 0;
    for i=1:tot_numyrs
        if(i_start == 0)
            i_start = 1;
        else
            i_start = i_start+num_yr(i-1);
        end
        i_stop = i_stop+num_yr(i);
        
        numday=365;
        if( (i_yr/4-floor(i_yr/4)) == 0) 
            numday=366; 
        end
        
        time = 1/tstep_per_day:(1/tstep_per_day):numday; 
        time = time';      
        t_all(i_start:i_stop,1) = i_yr+time/(60*60*24)/numday;
        
        % advance year
        i_yr = i_yr+1;
    end   
    % Read in model output and assign only last year of data to vectors 
    load /NoBackup/LSM_single_out/Santarem_83_nosoilevap.9995;
    fort=Santarem_83_nosoilevap; clear Santarem_83_nosoilevap
    fort = fort(i_start:i_stop,:); 
    load /NoBackup/LSM_single_out/Santarem_83_soilevap_or.9995;
    fort2=Santarem_83_soilevap_or; clear Santarem_83_soilevap_or;
    fort2=fort2(i_start:i_stop,:);
    load /NoBackup/LSM_single_out/Santarem_83_soilevap_new.9995;
    fort3=Santarem_83_soilevap_new; clear Santarem_83_soilevap_new;
    fort3=fort3(i_start:i_stop,:);
    
    % Create time vector for the plotted year
    t = t_all(i_start:i_stop,:);

    % Daily time series
    i_start = 1;
    i_stop = num_yr(plot_yr_num)/tstep_per_day;
    time=(1.5:num_yr(plot_yr_num)/tstep_per_day+.5)';
    t_day(i_start:i_stop,1) = plot_yr+time/num_day_plot_yr;    

%--------------------------------------------------------------------------
elseif(site_num == 5) % Harvard Forest, MA, USA
    % Read in ecohydrologic data from Fluxnet sites for comparison/evaluation
    % of model data
    yr_start = 1983;
    tot_numyrs = 1; 
    num_leapyr = 0;
    num_spinup_leap = 0;
    num_spinup_nonleap = 0;
    plot_yr = 1983;
    plot_yr_num = 1;  
    num_day_plot_yr=365; 
    num_data_days=365;
    
    mnth_days = [31 28 31 30 31 30 31 31 30 31 30 31];
    tstep_per_day =48;    nee_steps = 24;
    
    i_yr = yr_start;
    for i=1:tot_numyrs
        numday=365;
        if( (i_yr/4-floor(i_yr/4)) == 0) 
            numday=366; 
        end
        num_yr(i) = tstep_per_day*numday;
        i_yr = i_yr+1;
    end
    
    tstep_to_day = tstep_per_day;
    totnum_days = (tot_numyrs-num_spinup_leap-num_spinup_nonleap)*365 + num_leapyr*366;
    tstep_to_yr = tstep_per_day*totnum_days;
    num_tsteps =totnum_days*tstep_per_day; 
    t_subday = ((1/tstep_per_day):(1/tstep_per_day):num_data_days)/num_day_plot_yr; 
    t_subday = plot_yr+t_subday';
    
    % Latent heat, sensible heat, and NEE data
    filename = '/NoBackup/FLUXNET_reto/selected/netcdf/ameriflux.Harvard_Forest.2002.nc';
    ncload (filename,'LE'); ncload (filename,'H'); ncload (filename,'NEE');
    ncload (filename,'Rnet');
    filename = '/NoBackup/FLUXNET_reto/selected/alma/Harvard_Forest/Harvard_Forest_1999.nc';
    ncload (filename,'SWdown');
    LE_Wm2=LE;
    H_Wm2=H;
    NEE_umolm2s=NEE;  
    time=(0.5:num_yr(plot_yr_num)/tstep_per_day-.5)';
    tdata_day(1:num_day_plot_yr,1) = plot_yr+time/num_day_plot_yr;
    t_hr = ((1/nee_steps):(1/nee_steps):num_data_days)/num_day_plot_yr; 
    t_hr = plot_yr+t_hr';
    
    % Soil moisture data
    %load ????
    soilmoist1 = zeros(1,size(t_subday,1));
    soilmoist2 = zeros(1,size(t_subday,1));
    soilmoist3 = zeros(1,size(t_subday,1));

    t_sm = 1/nee_steps:(1/nee_steps):num_data_days*(1/nee_steps); 
    t_sm = t_sm';
    
    % Create time vector
    i_yr = yr_start;
    i_start = 0; i_stop = 0;
    for i=1:tot_numyrs
        if(i_start == 0)
            i_start = 1;
        else
            i_start = i_start+num_yr(i-1);
        end
        i_stop = i_stop+num_yr(i);
        
        numday=365;
        if( (i_yr/4-floor(i_yr/4)) == 0) 
            numday=366; 
        end
        
        time = 1/tstep_per_day:(1/tstep_per_day):numday; 
        time = time';      
        t_all(i_start:i_stop,1) = i_yr+time/(60*60*24)/numday;
        
        % advance year
        i_yr = i_yr+1;
    end   
    % Read in model output and assign only last year of data to vectors 
    load /NoBackup/LSM_single_out/Global_harvard.9995;
    fort=Global_harvard; clear Global_harvard
    fort = fort(i_start:i_stop,:); 
    load /NoBackup/LSM_single_out/Global_harvard2.9995;
    fort2=Global_harvard2; clear Global_harvard2;
    fort2=fort2(i_start:i_stop,:);
    load /NoBackup/LSM_single_out/Global_Harvard3.9995;
    fort3=Global_Harvard3; clear Global_Harvard3;
    fort3=fort3(i_start:i_stop,:);
    
    % Create time vector for the plotted year
    t = t_all(i_start:i_stop,:);

    % Daily time series
    i_start = 1;
    i_stop = num_yr(plot_yr_num)/tstep_per_day;
    time=(1.5:num_yr(plot_yr_num)/tstep_per_day+.5)';
    t_day(i_start:i_stop,1) = plot_yr+time/num_day_plot_yr;  
  
end

t_mon = 1:12;

%--------------------------------------------------------------------------
% Parameters, constants, and conversion factors
dz=[0.09999;0.1725;0.2977;0.5137;0.8864;1.529];% meters
lambdaLE = 2454000;	% latent heat of vaporization [J kg-1 @ 20 C]
rhow = 1000; % water density kg/m3
kg_m2_to_mm= 1/rhow*1000; %kg/m2=> m => mm
kg_to_umol = 1000/12.011*1000000;
umolm2s_to_molm2d=1/1000000*60*60*24;

%--------------------------------------------------------------------------
% Net Ecosystem Exchange data 
NEE_data=NEE_umolm2s /1000000*86400; %Convert from umol/m2/s to mol/m2/day
NEE_data(NEE_umolm2s==-9999)=-9999;
num_thrsteps = size(NEE_umolm2s,1);

i=1; j=1; jcount=0; NEE_data_d=zeros(num_data_days,1); 
n_count = zeros(num_data_days,1);
while (i<=size(t_hr,1))
    if (jcount<nee_steps)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    if(NEE_data(i)>-9999)
        NEE_data_d(j)=NEE_data_d(j)+NEE_data(i);
        n_count(j)=n_count(j)+1;
    end
    i = i+1;
end
NEE_data_d = NEE_data_d./n_count;
NEE_data_mnth=avg_day_to_mnth(mnth_days, NEE_data_d, num_data_days);

NEE_miss =NEE_umolm2s; NEE_miss(NEE_umolm2s~=-9999)=NaN; NEE_miss(NEE_umolm2s==-9999)=0;
NEE_data(NEE_umolm2s==-9999)=NaN;

%--------------------------------------------------------------------------
% Sensible heat data
i=1; j=1; jcount=0; H_Wm2_d=zeros(num_data_days,1);
n_count = zeros(num_data_days,1);
while (i<=size(t_hr,1))
    if (jcount<nee_steps)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    if(H_Wm2(i)>-9999)
        H_Wm2_d(j)=H_Wm2_d(j)+H_Wm2(i);
        n_count(j)=n_count(j)+1;
    end
    i = i+1;
end
H_Wm2_d = H_Wm2_d./n_count;
H_Wm2_mnth=avg_day_to_mnth(mnth_days, H_Wm2_d, num_data_days);

H_Wm2_miss =H_Wm2; H_Wm2_miss(H_Wm2~=-9999)=NaN; H_Wm2_miss(H_Wm2==-9999)=0;
H_Wm2(H_Wm2==-9999)=NaN;
%H_Wm2(find(H_Wm2==-9999))=NaN;
%--------------------------------------------------------------------------
% Latent heat data
i=1; j=1; jcount=0; LE_Wm2_d=zeros(num_data_days,1);
n_count = zeros(num_data_days,1);
while (i<=size(t_hr,1))
    if (jcount<nee_steps)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    if(LE_Wm2(i)>-9999)
        LE_Wm2_d(j)=LE_Wm2_d(j)+LE_Wm2(i);
        n_count(j)=n_count(j)+1;
    end
    i = i+1;
end
LE_Wm2_d = LE_Wm2_d./n_count;
LE_Wm2_mnth=avg_day_to_mnth(mnth_days, LE_Wm2_d, num_data_days);

LE_Wm2_miss =H_Wm2; LE_Wm2_miss(LE_Wm2~=-9999)=NaN; LE_Wm2_miss(LE_Wm2==-9999)=0;
LE_Wm2(LE_Wm2==-9999)=NaN;

%LE_Wm2(find(LE_Wm2==-9999))=NaN;
%LE_Wm2_d(find(n_count<nee_steps/2))=NaN;

% Convert LE to evapotranspiration
ET_data = LE_Wm2_d/lambdaLE*kg_m2_to_mm*24*60*60;
ET_data_mnth=avg_day_to_mnth(mnth_days, ET_data, num_data_days);
%--------------------------------------------------------------------------

 %kg C/m2/30min->kg/m2/s  -> g/m2/s -> mol/m2/s ->mol/m2/day
agpp  = fort(:,21)/(30*60) * 1000   / 12.011 * 86400;
i=1; j=1; jcount=0; agpp_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    agpp_d(j)=agpp_d(j)+agpp(i);
    i = i+1;
end
agpp_d = agpp_d/tstep_per_day;
agpp_mnth=avg_day_to_mnth(mnth_days, agpp_d, num_data_days);

 % kg C/m2/30min  ->kg/m2/s -> g/m2/s-> mol/m2/s ->mol/m2/day
arauto  = fort(:,22)/(30*60) * 1000   / 12.011 * 86400;
i=1; j=1; jcount=0; arauto_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    arauto_d(j)=arauto_d(j)+arauto(i);
    i = i+1;
end
arauto_d = arauto_d/tstep_per_day;
arauto_mnth=avg_day_to_mnth(mnth_days, arauto_d, num_data_days);

 % kg C/m2/30min  ->kg/m2/s -> g/m2/s-> mol/m2/s ->mol/m2/day
asoilresp = fort(:,23)/(30*60) * 1000   / 12.011 * 86400;
i=1; j=1; jcount=0; asoilresp_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    asoilresp_d(j)=asoilresp_d(j)+asoilresp(i);
    i = i+1;
end
asoilresp_d = asoilresp_d/tstep_per_day;
asoilresp_mnth=avg_day_to_mnth(mnth_days, asoilresp_d, num_data_days);


anpp = agpp-arauto;
anpp_d = agpp_d-arauto_d;
anpp_mnth=avg_day_to_mnth(mnth_days, anpp_d, num_data_days);

anee=-(agpp-arauto-asoilresp);
i=1; j=1; jcount=0; anee_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    anee_d(j)=anee_d(j)+anee(i);
    i = i+1;
end
anee_d = anee_d/tstep_per_day;
anee_mnth=avg_day_to_mnth(mnth_days, anee_d, num_data_days);


aclab  = fort(:,31);
i=1; j=1; jcount=0; aclab_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aclab_d(j)=aclab_d(j)+aclab(i);
    i = i+1;
end
aclab_d = aclab_d/tstep_per_day;
aclab_mnth=avg_day_to_mnth(mnth_days, aclab_d, num_data_days);

asoilCpoolsum  = fort(:,32);
i=1; j=1; jcount=0; asoilCpoolsum_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    asoilCpoolsum_d(j)=asoilCpoolsum_d(j)+asoilCpoolsum(i);
    i = i+1;
end
asoilCpoolsum_d = asoilCpoolsum_d/tstep_per_day;

aClivepool_leaf  = fort(:,85);
i=1; j=1; jcount=0; aClivepool_leaf_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aClivepool_leaf_d(j)=aClivepool_leaf_d(j)+aClivepool_leaf(i);
    i = i+1;
end
aClivepool_leaf_d = aClivepool_leaf_d/tstep_per_day;

aClivepool_froot = fort(:,86);
i=1; j=1; jcount=0; aClivepool_froot_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aClivepool_froot_d(j)=aClivepool_froot_d(j)+aClivepool_froot(i);
    i = i+1;
end
aClivepool_froot_d = aClivepool_froot_d/tstep_per_day;

aClivepool_wood = fort(:,87);
i=1; j=1; jcount=0; aClivepool_wood_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aClivepool_wood_d(j)=aClivepool_wood_d(j)+aClivepool_wood(i);
    i = i+1;
end
aClivepool_wood_d = aClivepool_wood_d/tstep_per_day;


theta_v1  = fort(:,44)/dz(1);
i=1; j=1; jcount=0; theta_v1_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    theta_v1_d(j)=theta_v1_d(j)+theta_v1(i);
    i = i+1;
end
theta_v1_d = theta_v1_d/tstep_per_day;
theta_v1_mnth=avg_day_to_mnth(mnth_days, theta_v1_d, num_data_days);

theta_v2  = fort(:,45)/dz(2);
i=1; j=1; jcount=0; theta_v2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    theta_v2_d(j)=theta_v2_d(j)+theta_v2(i);
    i = i+1;
end
theta_v2_d = theta_v2_d/tstep_per_day;
theta_v2_mnth=avg_day_to_mnth(mnth_days, theta_v2_d, num_data_days);

theta_v3 = fort(:,46)/dz(3);
i=1; j=1; jcount=0; theta_v3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    theta_v3_d(j)=theta_v3_d(j)+theta_v3(i);
    i = i+1;
end
theta_v3_d = theta_v3_d/tstep_per_day;
theta_v3_mnth=avg_day_to_mnth(mnth_days, theta_v3_d, num_data_days);

theta_v4 = fort(:,47)/dz(4);
i=1; j=1; jcount=0; theta_v4_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    theta_v4_d(j)=theta_v4_d(j)+theta_v4(i);
    i = i+1;
end
theta_v4_d = theta_v4_d/tstep_per_day;
theta_v4_mnth=avg_day_to_mnth(mnth_days, theta_v4_d, num_data_days);


theta_v5 = fort(:,48)/dz(5);
i=1; j=1; jcount=0; theta_v5_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    theta_v5_d(j)=theta_v5_d(j)+theta_v5(i);
    i = i+1;
end
theta_v5_d = theta_v5_d/tstep_per_day;
theta_v5_mnth=avg_day_to_mnth(mnth_days, theta_v5_d, num_data_days);

theta_v6 = fort(:,49)/dz(6);
i=1; j=1; jcount=0; theta_v6_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    theta_v6_d(j)=theta_v6_d(j)+theta_v6(i);
    i = i+1;
end
theta_v6_d = theta_v6_d/tstep_per_day;
theta_v6_mnth=avg_day_to_mnth(mnth_days, theta_v6_d, num_data_days);


aruns  = fort(:,19)*kg_m2_to_mm*tstep_to_day; %rate mm/d for each 
i=1; j=1; jcount=0; aruns_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aruns_d(j)=aruns_d(j)+aruns(i);
    i = i+1;
end
aruns_d = aruns_d/tstep_per_day; %average rate per day (mm/d) 
aruns_mnth=avg_day_to_mnth(mnth_days, aruns_d, num_data_days);

arunu = fort(:,20)*kg_m2_to_mm*tstep_to_day;
i=1; j=1; jcount=0; arunu_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    arunu_d(j)=arunu_d(j)+arunu(i);
    i = i+1;
end
arunu_d = arunu_d/tstep_per_day;
arunu_mnth=avg_day_to_mnth(mnth_days, arunu_d, num_data_days);

abetad = fort(:,24);
i=1; j=1; jcount=0; abetad_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    abetad_d(j)=abetad_d(j)+abetad(i);
    i = i+1;
end
abetad_d = abetad_d/tstep_per_day;
abetad_mnth=avg_day_to_mnth(mnth_days, abetad_d, num_data_days);

% Input meteorological drivers variables
precip  = fort(:,1);
i=1; j=1; jcount=0; precip_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j = j+1; jcount = 1;
    end
    precip_d(j)=precip_d(j)+precip(i);
    i = i+1;
end
precip_d = precip_d/tstep_per_day;
precip_d = precip_d*1000*(24*60*60);

precip = precip*1000*(24*60*60);
precip_mnth=avg_day_to_mnth(mnth_days, precip_d, num_data_days);


srheat  = fort(:,5);
i=1; j=1; jcount=0; srheat_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    srheat_d(j)=srheat_d(j)+srheat(i);
    i = i+1;
end
srheat_d = srheat_d/tstep_per_day;
srheat_mnth=avg_day_to_mnth(mnth_days, srheat_d, num_data_days);

trheat  = fort(:,6);
i=1; j=1; jcount=0; trheat_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    trheat_d(j)=trheat_d(j)+trheat(i);
    i = i+1;
end
trheat_d = trheat_d/tstep_per_day;
trheat_mnth=avg_day_to_mnth(mnth_days, trheat_d, num_data_days);

ts = fort(:,7); ts = ts-273.15;
i=1; j=1; jcount=0; ts_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    ts_d(j)=ts_d(j)+ts(i);
    i = i+1;
end
ts_d = ts_d/tstep_per_day;
ts_mnth=avg_day_to_mnth(mnth_days, ts_d, num_data_days);

qs  = fort(:,8);
i=1; j=1; jcount=0; qs_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    qs_d(j)=qs_d(j)+qs(i);
    i = i+1;
end
qs_d = qs_d/tstep_per_day;
qs_mnth=avg_day_to_mnth(mnth_days, qs_d, num_data_days);

ps = fort(:,9);
i=1; j=1; jcount=0; ps_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    ps_d(j)=ps_d(j)+ps(i);
    i = i+1;
end
ps_d = ps_d/tstep_per_day;
ps_mnth=avg_day_to_mnth(mnth_days, ps_d, num_data_days);

rhosurf = fort(:,10);
i=1; j=1; jcount=0; rhosurf_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    rhosurf_d(j)=rhosurf_d(j)+rhosurf(i);
    i = i+1;
end
rhosurf_d = rhosurf_d/tstep_per_day;
rhosurf_mnth=avg_day_to_mnth(mnth_days, rhosurf_d, num_data_days);

cdh = fort(:,11);
i=1; j=1; jcount=0; cdh_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    cdh_d(j)=cdh_d(j)+cdh(i);
    i = i+1;
end
cdh_d = cdh_d/tstep_per_day;
cdh_mnth=avg_day_to_mnth(mnth_days, cdh_d, num_data_days);


ws = fort(:,13);
i=1; j=1; jcount=0; ws_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    ws_d(j)=ws_d(j)+ws(i);
    i = i+1;
end
ws_d = ws_d/tstep_per_day;
ws_mnth=avg_day_to_mnth(mnth_days, ws_d, num_data_days);

% Evapotranspiration partitioning: values are total per dt
aevap  = fort(:,14)*tstep_to_day*kg_m2_to_mm; % mm/d
i=1; j=1; jcount=0; aevap_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aevap_d(j)=aevap_d(j)+aevap(i);
    i = i+1;
end
aevap_d = aevap_d/tstep_per_day;
aevap_mnth=avg_day_to_mnth(mnth_days, aevap_d, num_data_days);

aepp = fort(:,33)*tstep_to_day*kg_m2_to_mm; % mm/d
i=1; j=1; jcount=0; aepp_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aepp_d(j)=aepp_d(j)+aepp(i);
    i = i+1;
end
aepp_d = aepp_d/tstep_per_day;
aepp_mnth=avg_day_to_mnth(mnth_days, aepp_d, num_data_days);


aepc  = fort(:,37)*kg_m2_to_mm*tstep_to_day;
i=1; j=1; jcount=0; aepc_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aepc_d(j)=aepc_d(j)+aepc(i);
    i = i+1;
end
aepc_d = aepc_d/tstep_per_day;
aepc_mnth=avg_day_to_mnth(mnth_days, aepc_d, num_data_days);

aevapw  = fort(:,15)*kg_m2_to_mm*tstep_to_day;
i=1; j=1; jcount=0; aevapw_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aevapw_d(j)=aevapw_d(j)+aevapw(i);
    i = i+1;
end
aevapw_d = aevapw_d/tstep_per_day;
aevapw_mnth=avg_day_to_mnth(mnth_days, aevapw_d, num_data_days);

aevapd  = fort(:,16)*kg_m2_to_mm*tstep_to_day;
i=1; j=1; jcount=0; aevapd_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aevapd_d(j)=aevapd_d(j)+aevapd(i);
    i = i+1;
end
aevapd_d = aevapd_d/tstep_per_day;
aevapd_mnth=avg_day_to_mnth(mnth_days, aevapd_d, num_data_days);

aevapb  = fort(:,17)*kg_m2_to_mm*tstep_to_day;
i=1; j=1; jcount=0; aevapb_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aevapb_d(j)=aevapb_d(j)+aevapb(i);
    i = i+1;
end
aevapb_d = aevapb_d/tstep_per_day;
aevapb_mnth=avg_day_to_mnth(mnth_days, aevapb_d, num_data_days);

aevapvg  = fort(:,25)*kg_m2_to_mm*tstep_to_day;
i=1; j=1; jcount=0; aevapvg_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aevapvg_d(j)=aevapvg_d(j)+aevapvg(i);
    i = i+1;
end
aevapvg_d = aevapvg_d/tstep_per_day;
aevapvg_mnth=avg_day_to_mnth(mnth_days, aevapvg_d, num_data_days);

aevapvs  = fort(:,26)*kg_m2_to_mm*tstep_to_day;
i=1; j=1; jcount=0; aevapvs_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aevapvs_d(j)=aevapvs_d(j)+aevapvs(i);
    i = i+1;
end
aevapvs_d = aevapvs_d/tstep_per_day;
aevapvs_mnth=avg_day_to_mnth(mnth_days, aevapvs_d, num_data_days);

aevapvb = fort(:,27)*kg_m2_to_mm*tstep_to_day;
i=1; j=1; jcount=0; aevapvb_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aevapvb_d(j)=aevapvb_d(j)+aevapvb(i);
    i = i+1;
end
aevapvb_d = aevapvb_d/tstep_per_day;
aevapvb_mnth=avg_day_to_mnth(mnth_days, aevapvb_d, num_data_days);

aintercep = fort(:,18)*kg_m2_to_mm*tstep_to_day;
i=1; j=1; jcount=0; aintercep_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aintercep_d(j)=aintercep_d(j)+aintercep(i);
    i = i+1;
end
aintercep_d = aintercep_d/tstep_per_day;
aintercep_mnth=avg_day_to_mnth(mnth_days, aintercep_d, num_data_days);

lai = fort(:,97);
i=1; j=1; jcount=0; lai_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    lai_d(j)=lai_d(j)+lai(i);
    i = i+1;
end
lai_d = lai_d/tstep_per_day;

canopy_wat = fort(:,98);
i=1; j=1; jcount=0; canopy_wat_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    canopy_wat_d(j)=canopy_wat_d(j)+canopy_wat(i);
    i = i+1;
end
canopy_wat_d = canopy_wat_d/tstep_per_day;

canopy_heat = fort(:,99);
i=1; j=1; jcount=0; canopy_heat_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    canopy_heat_d(j)=canopy_heat_d(j)+canopy_heat(i);
    i = i+1;
end
canopy_heat_d = canopy_heat_d/tstep_per_day;

GCANOPY= fort(:,100);
i=1; j=1; jcount=0; GCANOPY_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    GCANOPY_d(j)=GCANOPY_d(j)+GCANOPY(i);
    i = i+1;
end
GCANOPY_d = GCANOPY_d/tstep_per_day;
GCANOPY_mnth=avg_day_to_mnth(mnth_days, GCANOPY_d, num_data_days);


ashg= fort(:,35)/((24/tstep_per_day)*60*60); %(J/m2 summed for dtsec)/dtsec=J/m2/s=W/m2
i=1; j=1; jcount=0; ashg_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    ashg_d(j)=ashg_d(j)+ashg(i);
    i = i+1;
end
ashg_d = ashg_d/tstep_per_day;
ashg_mnth=avg_day_to_mnth(mnth_days, ashg_d, num_data_days);

% len=size(canopy_heat,1);
% ht = canopy_heat;
% w=canopy_wat;
% rhow = 1000; %kg/m3 density of water
% lhm = 334590; %J/kg latent heat of melting at 0 deg C
% shw_kg = 4185; %J/kg C heat capacity of water at 20 C
% shi_kg = 2060; %J/kg heat capacity of pure ice at 0C
% ht_cap = (0.01 + 0.002*LAI_site+ 0.001*LAI_site^2)*shw_kg*rhow;
% lhmv=lhm*rhow; shwv=shw_kg*rhow; shiv=shi_kg*rhow; % volumetric quantities
% tp = ht./( ht_cap + w.*shwv );


Qf = fort(:,84);
i=1; j=1; jcount=0; Qf_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    Qf_d(j)=Qf_d(j)+Qf(i);
    i = i+1;
end
Qf_d = Qf_d/tstep_per_day;
Qf_mnth=avg_day_to_mnth(mnth_days, Qf_d, num_data_days);

tcanopy = fort(:,103);
i=1; j=1; jcount=0; tcanopy_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    tcanopy_d(j)=tcanopy_d(j)+tcanopy(i);
    i = i+1;
end
tcanopy_d = tcanopy_d/tstep_per_day;
tcanopy_mnth=avg_day_to_mnth(mnth_days, tcanopy_d, num_data_days);


IPARtot = fort(:,101);
i=1; j=1; jcount=0; IPARtot_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    IPARtot_d(j)=IPARtot_d(j)+IPARtot(i);
    i = i+1;
end
IPARtot_d = IPARtot_d/tstep_per_day;
IPARtot_mnth=avg_day_to_mnth(mnth_days, IPARtot_d, num_data_days);


IPARdir = fort(:,102);
i=1; j=1; jcount=0; IPARdir_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    IPARdir_d(j)=IPARdir_d(j)+IPARdir(i);
    i = i+1;
end
IPARdir_d = IPARdir_d/tstep_per_day;
IPARdir_mnth=avg_day_to_mnth(mnth_days, IPARdir_d, num_data_days);

% Snow - water equivalent depth [m]  MULTIPLY BY SNOW FRACTION!!!!!
wsn_vtot  = fort(:,71)+fort(:,72)+fort(:,73);
i=1; j=1; jcount=0; wsn_vtot_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    wsn_vtot_d(j)=wsn_vtot_d(j)+wsn_vtot(i);
    i = i+1;
end
wsn_vtot_d = wsn_vtot_d/tstep_per_day;
wsn_vtot_mnth=avg_day_to_mnth(mnth_days, wsn_vtot_d, num_data_days);

%--------------------------------------------------------------------------
clear fort;
%--------------------------------------------------------------------------
Qf2 = fort2(:,84);
i=1; j=1; jcount=0; Qf2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    Qf2_d(j)=Qf2_d(j)+Qf2(i);
    i = i+1;
end
Qf2_d = Qf2_d/tstep_per_day;
Qf2_mnth=avg_day_to_mnth(mnth_days, Qf2_d, num_data_days);

% Input meteorological drivers variables
precip2  = fort2(:,1);
i=1; j=1; jcount=0; precip2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j = j+1; jcount = 1;
    end
    precip2_d(j)=precip2_d(j)+precip2(i);
    i = i+1;
end
precip2_d = precip2_d/tstep_per_day;
precip2_d = precip2_d*1000*(24*60*60);

precip2 = precip2*1000*(24*60*60);
precip2_mnth=avg_day_to_mnth(mnth_days, precip2_d, num_data_days);


srheat2  = fort2(:,5);
i=1; j=1; jcount=0; srheat2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    srheat2_d(j)=srheat2_d(j)+srheat2(i);
    i = i+1;
end
srheat2_d = srheat2_d/tstep_per_day;
srheat2_mnth=avg_day_to_mnth(mnth_days, srheat2_d, num_data_days);

trheat2  = fort2(:,6);
i=1; j=1; jcount=0; trheat2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    trheat2_d(j)=trheat2_d(j)+trheat2(i);
    i = i+1;
end
trheat2_d = trheat2_d/tstep_per_day;
trheat2_mnth=avg_day_to_mnth(mnth_days, trheat2_d, num_data_days);

ts2 = fort2(:,7); ts2 = ts2-273.15;
i=1; j=1; jcount=0; ts2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    ts2_d(j)=ts2_d(j)+ts2(i);
    i = i+1;
end
ts2_d = ts2_d/tstep_per_day;
ts2_mnth=avg_day_to_mnth(mnth_days, ts2_d, num_data_days);

qs2  = fort2(:,8);
i=1; j=1; jcount=0; qs2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    qs2_d(j)=qs2_d(j)+qs2(i);
    i = i+1;
end
qs2_d = qs2_d/tstep_per_day;
qs2_mnth=avg_day_to_mnth(mnth_days, qs2_d, num_data_days);

ps2 = fort2(:,9);
i=1; j=1; jcount=0; ps2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    ps2_d(j)=ps2_d(j)+ps2(i);
    i = i+1;
end
ps2_d = ps2_d/tstep_per_day;
ps2_mnth=avg_day_to_mnth(mnth_days, ps2_d, num_data_days);

rhosurf2 = fort2(:,10);
i=1; j=1; jcount=0; rhosurf2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    rhosurf2_d(j)=rhosurf2_d(j)+rhosurf2(i);
    i = i+1;
end
rhosurf2_d = rhosurf2_d/tstep_per_day;
rhosurf2_mnth=avg_day_to_mnth(mnth_days, rhosurf2_d, num_data_days);

cdh2 = fort2(:,11);
i=1; j=1; jcount=0; cdh2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    cdh2_d(j)=cdh2_d(j)+cdh2(i);
    i = i+1;
end
cdh2_d = cdh2_d/tstep_per_day;
cdh2_mnth=avg_day_to_mnth(mnth_days, cdh2_d, num_data_days);


ws2 = fort2(:,13);
i=1; j=1; jcount=0; ws2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    ws2_d(j)=ws2_d(j)+ws2(i);
    i = i+1;
end
ws2_d = ws2_d/tstep_per_day;
ws2_mnth=avg_day_to_mnth(mnth_days, ws2_d, num_data_days);


%kg C/m2/30min->kg/m2/s  -> g/m2/s -> mol/m2/s ->mol/m2/day
agpp2  = fort2(:,21)/(30*60) * 1000   / 12.011 * 86400;
i=1; j=1; jcount=0; agpp2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    agpp2_d(j)=agpp2_d(j)+agpp2(i);
    i = i+1;
end
agpp2_d = agpp2_d/tstep_per_day;
agpp2_mnth=avg_day_to_mnth(mnth_days, agpp2_d, num_data_days);

 % kg C/m2/30min  ->kg/m2/s -> g/m2/s-> mol/m2/s ->mol/m2/day
arauto2  = fort2(:,22)/(30*60) * 1000   / 12.011 * 86400;
i=1; j=1; jcount=0; arauto2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    arauto2_d(j)=arauto2_d(j)+arauto2(i);
    i = i+1;
end
arauto2_d = arauto2_d/tstep_per_day;

 % kg C/m2/30min  ->kg/m2/s -> g/m2/s-> mol/m2/s ->mol/m2/day
asoilresp2 = fort2(:,23)/(30*60) * 1000   / 12.011 * 86400;
i=1; j=1; jcount=0; asoilresp2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    asoilresp2_d(j)=asoilresp2_d(j)+asoilresp2(i);
    i = i+1;
end
asoilresp2_d = asoilresp2_d/tstep_per_day;

anpp2 = agpp2-arauto2;
anpp2_d = agpp2_d-arauto2_d; 

anee2=-(agpp2-arauto2-asoilresp2);
i=1; j=1; jcount=0; anee2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    anee2_d(j)=anee2_d(j)+anee2(i);
    i = i+1;
end
anee2_d = anee2_d/tstep_per_day;

aclab2  = fort2(:,31);
i=1; j=1; jcount=0; aclab2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aclab2_d(j)=aclab2_d(j)+aclab2(i);
    i = i+1;
end
aclab2_d = aclab2_d/tstep_per_day;

asoilCpoolsum2  = fort2(:,32);
i=1; j=1; jcount=0; asoilCpoolsum2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    asoilCpoolsum2_d(j)=asoilCpoolsum2_d(j)+asoilCpoolsum2(i);
    i = i+1;
end
asoilCpoolsum2_d = asoilCpoolsum2_d/tstep_per_day;


aClivepool_leaf2  = fort2(:,85);
i=1; j=1; jcount=0; aClivepool_leaf2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aClivepool_leaf2_d(j)=aClivepool_leaf2_d(j)+aClivepool_leaf2(i);
    i = i+1;
end
aClivepool_leaf2_d = aClivepool_leaf2_d/tstep_per_day;

aClivepool_froot2 = fort2(:,86);
i=1; j=1; jcount=0; aClivepool_froot2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aClivepool_froot2_d(j)=aClivepool_froot2_d(j)+aClivepool_froot2(i);
    i = i+1;
end
aClivepool_froot2_d = aClivepool_froot2_d/tstep_per_day;

aClivepool_wood2 = fort2(:,87);
i=1; j=1; jcount=0; aClivepool_wood2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aClivepool_wood2_d(j)=aClivepool_wood2_d(j)+aClivepool_wood2(i);
    i = i+1;
end
aClivepool_wood2_d = aClivepool_wood2_d/tstep_per_day;

theta2_v1  = fort2(:,44)/dz(1);
i=1; j=1; jcount=0; theta2_v1_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    theta2_v1_d(j)=theta2_v1_d(j)+theta2_v1(i);
    i = i+1;
end
theta2_v1_d = theta2_v1_d/tstep_per_day;

theta2_v2  = fort2(:,45)/dz(2);
i=1; j=1; jcount=0; theta2_v2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    theta2_v2_d(j)=theta2_v2_d(j)+theta2_v2(i);
    i = i+1;
end
theta2_v2_d = theta2_v2_d/tstep_per_day;

theta2_v3 = fort2(:,46)/dz(3);
i=1; j=1; jcount=0; theta2_v3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    theta2_v3_d(j)=theta2_v3_d(j)+theta2_v3(i);
    i = i+1;
end
theta2_v3_d = theta2_v3_d/tstep_per_day;

theta2_v4 = fort2(:,47)/dz(4);
i=1; j=1; jcount=0; theta2_v4_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    theta2_v4_d(j)=theta2_v4_d(j)+theta2_v4(i);
    i = i+1;
end
theta2_v4_d = theta2_v4_d/tstep_per_day;

theta2_v5 = fort2(:,48)/dz(5);
i=1; j=1; jcount=0; theta2_v5_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    theta2_v5_d(j)=theta2_v5_d(j)+theta2_v5(i);
    i = i+1;
end
theta2_v5_d = theta2_v5_d/tstep_per_day;

theta2_v6 = fort2(:,49)/dz(6);
i=1; j=1; jcount=0; theta2_v6_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    theta2_v6_d(j)=theta2_v6_d(j)+theta2_v6(i);
    i = i+1;
end
theta2_v6_d = theta2_v6_d/tstep_per_day;

aruns2  = fort2(:,19)*kg_m2_to_mm*tstep_to_day;
i=1; j=1; jcount=0; aruns2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aruns2_d(j)=aruns2_d(j)+aruns2(i);
    i = i+1;
end
aruns2_d = aruns2_d/tstep_per_day;

arunu2 = fort2(:,20)*kg_m2_to_mm*tstep_to_day;
i=1; j=1; jcount=0; arunu2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    arunu2_d(j)=arunu2_d(j)+arunu2(i);
    i = i+1;
end
arunu2_d = arunu2_d/tstep_per_day;

abetad2 = fort2(:,24);
i=1; j=1; jcount=0; abetad2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    abetad2_d(j)=abetad2_d(j)+abetad2(i);
    i = i+1;
end
abetad2_d = abetad2_d/tstep_per_day;

% Evapotranspiration partitioning: values are total per 30 min
aevap2  = fort2(:,14)*tstep_to_day*kg_m2_to_mm; % mm/d
i=1; j=1; jcount=0; aevap2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aevap2_d(j)=aevap2_d(j)+aevap2(i);
    i = i+1;
end
aevap2_d = aevap2_d/tstep_per_day;
aevap2_mnth=avg_day_to_mnth(mnth_days, aevap2_d, num_data_days);

aepp2 = fort2(:,33)*tstep_to_day*kg_m2_to_mm; % mm/d
i=1; j=1; jcount=0; aepp2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aepp2_d(j)=aepp2_d(j)+aepp2(i);
    i = i+1;
end
aepp2_d = aepp2_d/tstep_per_day;
aepp2_mnth=avg_day_to_mnth(mnth_days, aepp2_d, num_data_days);

aepc2  = fort2(:,37)*kg_m2_to_mm*tstep_to_day;
i=1; j=1; jcount=0; aepc2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aepc2_d(j)=aepc2_d(j)+aepc2(i);
    i = i+1;
end
aepc2_d = aepc2_d/tstep_per_day;
aepc2_mnth=avg_day_to_mnth(mnth_days, aepc2_d, num_data_days);

aevapw2  = fort2(:,15)*kg_m2_to_mm*tstep_to_day;
i=1; j=1; jcount=0; aevapw2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aevapw2_d(j)=aevapw2_d(j)+aevapw2(i);
    i = i+1;
end
aevapw2_d = aevapw2_d/tstep_per_day;
aevapw2_mnth=avg_day_to_mnth(mnth_days, aevapw2_d, num_data_days);

aevapd2  = fort2(:,16)*kg_m2_to_mm*tstep_to_day;
i=1; j=1; jcount=0; aevapd2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aevapd2_d(j)=aevapd2_d(j)+aevapd2(i);
    i = i+1;
end
aevapd2_d = aevapd2_d/tstep_per_day;
aevapd2_mnth=avg_day_to_mnth(mnth_days, aevapd2_d, num_data_days);

aevapb2  = fort2(:,17)*kg_m2_to_mm*tstep_to_day;
i=1; j=1; jcount=0; aevapb2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aevapb2_d(j)=aevapb2_d(j)+aevapb2(i);
    i = i+1;
end
aevapb2_d = aevapb2_d/tstep_per_day;
aevapb2_mnth=avg_day_to_mnth(mnth_days, aevapb2_d, num_data_days);

aevapvg2  = fort2(:,25)*kg_m2_to_mm*tstep_to_day;
i=1; j=1; jcount=0; aevapvg2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aevapvg2_d(j)=aevapvg2_d(j)+aevapvg2(i);
    i = i+1;
end
aevapvg2_d = aevapvg2_d/tstep_per_day;
aevapvg2_mnth=avg_day_to_mnth(mnth_days, aevapvg2_d, num_data_days);

aevapvs2  = fort2(:,26)*kg_m2_to_mm*tstep_to_day;
i=1; j=1; jcount=0; aevapvs2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aevapvs2_d(j)=aevapvs2_d(j)+aevapvs2(i);
    i = i+1;
end
aevapvs2_d = aevapvs2_d/tstep_per_day;
aevapvs2_mnth=avg_day_to_mnth(mnth_days, aevapvs2_d, num_data_days);

aevapvb2 = fort2(:,27)*kg_m2_to_mm*tstep_to_day;
i=1; j=1; jcount=0; aevapvb2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aevapvb2_d(j)=aevapvb2_d(j)+aevapvb2(i);
    i = i+1;
end
aevapvb2_d = aevapvb2_d/tstep_per_day;
aevapvb2_mnth=avg_day_to_mnth(mnth_days, aevapvb2_d, num_data_days);

aintercep2 = fort2(:,18)*kg_m2_to_mm*tstep_to_day;
i=1; j=1; jcount=0; aintercep2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aintercep2_d(j)=aintercep2_d(j)+aintercep2(i);
    i = i+1;
end
aintercep2_d = aintercep2_d/tstep_per_day;
aintercep2_mnth=avg_day_to_mnth(mnth_days, aintercep2_d, num_data_days);

lai2 = fort2(:,97);
i=1; j=1; jcount=0; lai2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    lai2_d(j)=lai2_d(j)+lai2(i);
    i = i+1;
end
lai2_d = lai2_d/tstep_per_day;

canopy_wat2 = fort2(:,98);
i=1; j=1; jcount=0; canopy_wat2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    canopy_wat2_d(j)=canopy_wat2_d(j)+canopy_wat2(i);
    i = i+1;
end
canopy_wat2_d = canopy_wat2_d/tstep_per_day;

canopy_heat2 = fort2(:,99);
i=1; j=1; jcount=0; canopy_heat2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    canopy_heat2_d(j)=canopy_heat2_d(j)+canopy_heat2(i);
    i = i+1;
end
canopy_heat2_d = canopy_heat2_d/tstep_per_day;

GCANOPY2= fort2(:,100);
i=1; j=1; jcount=0; GCANOPY2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    GCANOPY2_d(j)=GCANOPY2_d(j)+GCANOPY2(i);
    i = i+1;
end
GCANOPY2_d = GCANOPY2_d/tstep_per_day;

ashg2= fort2(:,35)/((24/tstep_per_day)*60*60); %(J/m2 summed for dtsec)/dtsec=J/m2/s=W/m2
i=1; j=1; jcount=0; ashg2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    ashg2_d(j)=ashg2_d(j)+ashg2(i);
    i = i+1;
end
ashg2_d = ashg2_d/tstep_per_day;
ashg2_mnth=avg_day_to_mnth(mnth_days, ashg2_d, num_data_days);

% % computes the average LAI for the vegetation 
% % (average of min and max)
% % averaged over all pfts in cell
% len=size(canopy_heat2,1);
% ht = canopy_heat2;
% w=canopy_wat2;
% tp2 = ht./( ht_cap + w.*shwv );

tcanopy2 = fort2(:,103);
i=1; j=1; jcount=0; tcanopy2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    tcanopy2_d(j)=tcanopy2_d(j)+tcanopy2(i);
    i = i+1;
end
tcanopy2_d = tcanopy2_d/tstep_per_day;
tcanopy2_mnth=avg_day_to_mnth(mnth_days, tcanopy2_d, num_data_days);


IPARtot2 = fort2(:,101);
i=1; j=1; jcount=0; IPARtot2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    IPARtot2_d(j)=IPARtot2_d(j)+IPARtot2(i);
    i = i+1;
end
IPARtot2_d = IPARtot2_d/tstep_per_day;
IPARtot2_mnth=avg_day_to_mnth(mnth_days, IPARtot2_d, num_data_days);


IPARdir2 = fort2(:,102);
i=1; j=1; jcount=0; IPARdir2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    IPARdir2_d(j)=IPARdir2_d(j)+IPARdir2(i);
    i = i+1;
end
IPARdir2_d = IPARdir2_d/tstep_per_day;
IPARdir2_mnth=avg_day_to_mnth(mnth_days, IPARdir2_d, num_data_days);

% Snow - water equivalent depth [m]
wsn_vtot2  = fort2(:,71)+fort2(:,72)+fort2(:,73);
i=1; j=1; jcount=0; wsn_vtot2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    wsn_vtot2_d(j)=wsn_vtot2_d(j)+wsn_vtot2(i);
    i = i+1;
end
wsn_vtot2_d = wsn_vtot2_d/tstep_per_day;
wsn_vtot2_mnth=avg_day_to_mnth(mnth_days, wsn_vtot2_d, num_data_days);
%--------------------------------------------------------------------------
clear fort2;
%--------------------------------------------------------------------------
% Input meteorological drivers variables
precip3  = fort3(:,1);
i=1; j=1; jcount=0; precip3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j = j+1; jcount = 1;
    end
    precip3_d(j)=precip3_d(j)+precip3(i);
    i = i+1;
end
precip3_d = precip3_d/tstep_per_day;
precip3_d = precip3_d*1000*(24*60*60);

precip3 = precip3*1000*(24*60*60);
precip3_mnth=avg_day_to_mnth(mnth_days, precip3_d, num_data_days);


srheat3  = fort3(:,5);
i=1; j=1; jcount=0; srheat3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    srheat3_d(j)=srheat3_d(j)+srheat3(i);
    i = i+1;
end
srheat3_d = srheat3_d/tstep_per_day;
srheat3_mnth=avg_day_to_mnth(mnth_days, srheat3_d, num_data_days);

trheat3  = fort3(:,6);
i=1; j=1; jcount=0; trheat3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    trheat3_d(j)=trheat3_d(j)+trheat3(i);
    i = i+1;
end
trheat3_d = trheat3_d/tstep_per_day;
trheat3_mnth=avg_day_to_mnth(mnth_days, trheat3_d, num_data_days);

ts3 = fort3(:,7); ts3 = ts3-273.15;
i=1; j=1; jcount=0; ts3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    ts3_d(j)=ts3_d(j)+ts3(i);
    i = i+1;
end
ts3_d = ts3_d/tstep_per_day;
ts3_mnth=avg_day_to_mnth(mnth_days, ts3_d, num_data_days);

qs3  = fort3(:,8);
i=1; j=1; jcount=0; qs3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    qs3_d(j)=qs3_d(j)+qs3(i);
    i = i+1;
end
qs3_d = qs3_d/tstep_per_day;
qs3_mnth=avg_day_to_mnth(mnth_days, qs3_d, num_data_days);

ps3 = fort3(:,9);
i=1; j=1; jcount=0; ps3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    ps3_d(j)=ps3_d(j)+ps3(i);
    i = i+1;
end
ps3_d = ps3_d/tstep_per_day;
ps3_mnth=avg_day_to_mnth(mnth_days, ps3_d, num_data_days);

rhosurf3 = fort3(:,10);
i=1; j=1; jcount=0; rhosurf3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    rhosurf3_d(j)=rhosurf3_d(j)+rhosurf3(i);
    i = i+1;
end
rhosurf3_d = rhosurf3_d/tstep_per_day;
rhosurf3_mnth=avg_day_to_mnth(mnth_days, rhosurf3_d, num_data_days);

cdh3 = fort3(:,11);
i=1; j=1; jcount=0; cdh3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    cdh3_d(j)=cdh3_d(j)+cdh3(i);
    i = i+1;
end
cdh3_d = cdh3_d/tstep_per_day;
cdh3_mnth=avg_day_to_mnth(mnth_days, cdh3_d, num_data_days);

ws3 = fort3(:,13);
i=1; j=1; jcount=0; ws3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    ws3_d(j)=ws3_d(j)+ws3(i);
    i = i+1;
end
ws3_d = ws3_d/tstep_per_day;
ws3_mnth=avg_day_to_mnth(mnth_days, ws3_d, num_data_days);

%kg C/m2/30min->kg/m2/s  -> g/m2/s -> mol/m2/s ->mol/m2/day
agpp3  = fort3(:,21)/(30*60) * 1000   / 12.011 * 86400;
i=1; j=1; jcount=0; agpp3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    agpp3_d(j)=agpp3_d(j)+agpp3(i);
    i = i+1;
end
agpp3_d = agpp3_d/tstep_per_day;
agpp3_mnth=avg_day_to_mnth(mnth_days, agpp3_d, num_data_days);

 % kg C/m2/30min  ->kg/m2/s -> g/m2/s-> mol/m2/s ->mol/m2/day
arauto3  = fort3(:,22)/(30*60) * 1000   / 12.011 * 86400;
i=1; j=1; jcount=0; arauto3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    arauto3_d(j)=arauto3_d(j)+arauto3(i);
    i = i+1;
end
arauto3_d = arauto3_d/tstep_per_day;

 % kg C/m2/30min  ->kg/m2/s -> g/m2/s-> mol/m2/s ->mol/m2/day
asoilresp3 = fort3(:,23)/(30*60) * 1000   / 12.011 * 86400;
i=1; j=1; jcount=0; asoilresp3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    asoilresp3_d(j)=asoilresp3_d(j)+asoilresp3(i);
    i = i+1;
end
asoilresp3_d = asoilresp3_d/tstep_per_day;

anpp3 = agpp3-arauto3;
anpp3_d = agpp3_d-arauto3_d; 

anee3=-(agpp3-arauto3-asoilresp3);
i=1; j=1; jcount=0; anee3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    anee3_d(j)=anee3_d(j)+anee3(i);
    i = i+1;
end
anee3_d = anee3_d/tstep_per_day;

aclab3  = fort3(:,31);
i=1; j=1; jcount=0; aclab3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aclab3_d(j)=aclab3_d(j)+aclab3(i);
    i = i+1;
end
aclab3_d = aclab3_d/tstep_per_day;

asoilCpoolsum3  = fort3(:,32);
i=1; j=1; jcount=0; asoilCpoolsum3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    asoilCpoolsum3_d(j)=asoilCpoolsum3_d(j)+asoilCpoolsum3(i);
    i = i+1;
end
asoilCpoolsum3_d = asoilCpoolsum3_d/tstep_per_day;


aClivepool_leaf3  = fort3(:,85);
i=1; j=1; jcount=0; aClivepool_leaf3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aClivepool_leaf3_d(j)=aClivepool_leaf3_d(j)+aClivepool_leaf3(i);
    i = i+1;
end
aClivepool_leaf3_d = aClivepool_leaf3_d/tstep_per_day;

aClivepool_froot3 = fort3(:,86);
i=1; j=1; jcount=0; aClivepool_froot3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aClivepool_froot3_d(j)=aClivepool_froot3_d(j)+aClivepool_froot3(i);
    i = i+1;
end
aClivepool_froot3_d = aClivepool_froot3_d/tstep_per_day;

aClivepool_wood3 = fort3(:,87);
i=1; j=1; jcount=0; aClivepool_wood3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aClivepool_wood3_d(j)=aClivepool_wood3_d(j)+aClivepool_wood3(i);
    i = i+1;
end
aClivepool_wood3_d = aClivepool_wood3_d/tstep_per_day;

theta3_v1  = fort3(:,44)/dz(1);
i=1; j=1; jcount=0; theta3_v1_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    theta3_v1_d(j)=theta3_v1_d(j)+theta3_v1(i);
    i = i+1;
end
theta3_v1_d = theta3_v1_d/tstep_per_day;

theta3_v2  = fort3(:,45)/dz(2);
i=1; j=1; jcount=0; theta3_v2_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    theta3_v2_d(j)=theta3_v2_d(j)+theta3_v2(i);
    i = i+1;
end
theta3_v2_d = theta3_v2_d/tstep_per_day;

theta3_v3 = fort3(:,46)/dz(3);
i=1; j=1; jcount=0; theta3_v3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    theta3_v3_d(j)=theta3_v3_d(j)+theta3_v3(i);
    i = i+1;
end
theta3_v3_d = theta3_v3_d/tstep_per_day;

theta3_v4 = fort3(:,47)/dz(4);
i=1; j=1; jcount=0; theta3_v4_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    theta3_v4_d(j)=theta3_v4_d(j)+theta3_v4(i);
    i = i+1;
end
theta3_v4_d = theta3_v4_d/tstep_per_day;

theta3_v5 = fort3(:,48)/dz(5);
i=1; j=1; jcount=0; theta3_v5_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    theta3_v5_d(j)=theta3_v5_d(j)+theta3_v5(i);
    i = i+1;
end
theta3_v5_d = theta3_v5_d/tstep_per_day;

theta3_v6 = fort3(:,49)/dz(6);
i=1; j=1; jcount=0; theta3_v6_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    theta3_v6_d(j)=theta3_v6_d(j)+theta3_v6(i);
    i = i+1;
end
theta3_v6_d = theta3_v6_d/tstep_per_day;

aruns3  = fort3(:,19)*kg_m2_to_mm*tstep_to_day;
i=1; j=1; jcount=0; aruns3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aruns3_d(j)=aruns3_d(j)+aruns3(i);
    i = i+1;
end
aruns3_d = aruns3_d/tstep_per_day;
aruns3_mnth=avg_day_to_mnth(mnth_days, aruns3_d, num_data_days);

arunu3 = fort3(:,20)*kg_m2_to_mm*tstep_to_day;
i=1; j=1; jcount=0; arunu3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    arunu3_d(j)=arunu3_d(j)+arunu3(i);
    i = i+1;
end
arunu3_d = arunu3_d/tstep_per_day;
arunu3_mnth=avg_day_to_mnth(mnth_days, arunu3_d, num_data_days);

abetad3 = fort3(:,24);
i=1; j=1; jcount=0; abetad3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    abetad3_d(j)=abetad3_d(j)+abetad3(i);
    i = i+1;
end
abetad3_d = abetad3_d/tstep_per_day;
abetad3_mnth=avg_day_to_mnth(mnth_days, abetad3_d, num_data_days);

% Evapotranspiration partitioning: values are total per 30 min
aevap3  = fort3(:,14)*tstep_to_day*kg_m2_to_mm; % mm/d
i=1; j=1; jcount=0; aevap3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aevap3_d(j)=aevap3_d(j)+aevap3(i);
    i = i+1;
end
aevap3_d = aevap3_d/tstep_per_day;
aevap3_mnth=avg_day_to_mnth(mnth_days, aevap3_d, num_data_days);

aepp3 = fort3(:,33)*tstep_to_day*kg_m2_to_mm; % mm/d
i=1; j=1; jcount=0; aepp3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aepp3_d(j)=aepp3_d(j)+aepp3(i);
    i = i+1;
end
aepp3_d = aepp3_d/tstep_per_day;
aepp3_mnth=avg_day_to_mnth(mnth_days, aepp3_d, num_data_days);

aepc3  = fort3(:,37)*kg_m2_to_mm*tstep_to_day;
i=1; j=1; jcount=0; aepc3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aepc3_d(j)=aepc3_d(j)+aepc3(i);
    i = i+1;
end
aepc3_d = aepc3_d/tstep_per_day;
aepc3_mnth=avg_day_to_mnth(mnth_days, aepc3_d, num_data_days);

aevapw3  = fort3(:,15)*kg_m2_to_mm*tstep_to_day;
i=1; j=1; jcount=0; aevapw3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aevapw3_d(j)=aevapw3_d(j)+aevapw3(i);
    i = i+1;
end
aevapw3_d = aevapw3_d/tstep_per_day;
aevapw3_mnth=avg_day_to_mnth(mnth_days, aevapw3_d, num_data_days);

aevapd3  = fort3(:,16)*kg_m2_to_mm*tstep_to_day;
i=1; j=1; jcount=0; aevapd3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aevapd3_d(j)=aevapd3_d(j)+aevapd3(i);
    i = i+1;
end
aevapd3_d = aevapd3_d/tstep_per_day;
aevapd3_mnth=avg_day_to_mnth(mnth_days, aevapd3_d, num_data_days);

aevapb3  = fort3(:,17)*kg_m2_to_mm*tstep_to_day;
i=1; j=1; jcount=0; aevapb3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aevapb3_d(j)=aevapb3_d(j)+aevapb3(i);
    i = i+1;
end
aevapb3_d = aevapb3_d/tstep_per_day;
aevapb3_mnth=avg_day_to_mnth(mnth_days, aevapb3_d, num_data_days);

aevapvg3  = fort3(:,25)*kg_m2_to_mm*tstep_to_day;
i=1; j=1; jcount=0; aevapvg3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aevapvg3_d(j)=aevapvg3_d(j)+aevapvg3(i);
    i = i+1;
end
aevapvg3_d = aevapvg3_d/tstep_per_day;
aevapvg3_mnth=avg_day_to_mnth(mnth_days, aevapvg3_d, num_data_days);

aevapvs3  = fort3(:,26)*kg_m2_to_mm*tstep_to_day;
i=1; j=1; jcount=0; aevapvs3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aevapvs3_d(j)=aevapvs3_d(j)+aevapvs3(i);
    i = i+1;
end
aevapvs3_d = aevapvs3_d/tstep_per_day;
aevapvs3_mnth=avg_day_to_mnth(mnth_days, aevapvs3_d, num_data_days);

aevapvb3 = fort3(:,27)*kg_m2_to_mm*tstep_to_day;
i=1; j=1; jcount=0; aevapvb3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aevapvb3_d(j)=aevapvb3_d(j)+aevapvb3(i);
    i = i+1;
end
aevapvb3_d = aevapvb3_d/tstep_per_day;
aevapvb3_mnth=avg_day_to_mnth(mnth_days, aevapvb3_d, num_data_days);

aintercep3 = fort3(:,18)*kg_m2_to_mm*tstep_to_day;
i=1; j=1; jcount=0; aintercep3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    aintercep3_d(j)=aintercep3_d(j)+aintercep3(i);
    i = i+1;
end
aintercep3_d = aintercep3_d/tstep_per_day;
aintercep3_mnth=avg_day_to_mnth(mnth_days, aintercep3_d, num_data_days);

lai3 = fort3(:,97);
i=1; j=1; jcount=0; lai3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    lai3_d(j)=lai3_d(j)+lai3(i);
    i = i+1;
end
lai3_d = lai3_d/tstep_per_day;

canopy_wat3 = fort3(:,98);
i=1; j=1; jcount=0; canopy_wat3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    canopy_wat3_d(j)=canopy_wat3_d(j)+canopy_wat3(i);
    i = i+1;
end
canopy_wat3_d = canopy_wat3_d/tstep_per_day;

canopy_heat3 = fort3(:,99);
i=1; j=1; jcount=0; canopy_heat3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    canopy_heat3_d(j)=canopy_heat3_d(j)+canopy_heat3(i);
    i = i+1;
end
canopy_heat3_d = canopy_heat3_d/tstep_per_day;

GCANOPY3= fort3(:,100);
i=1; j=1; jcount=0; GCANOPY3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    GCANOPY3_d(j)=GCANOPY3_d(j)+GCANOPY3(i);
    i = i+1;
end
GCANOPY3_d = GCANOPY3_d/tstep_per_day;


ashg3= fort3(:,35)/((24/tstep_per_day)*60*60); %(J/m2 summed for dtsec)/dtsec=J/m2/s=W/m2
i=1; j=1; jcount=0; ashg3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    ashg3_d(j)=ashg3_d(j)+ashg3(i);
    i = i+1;
end
ashg3_d = ashg3_d/tstep_per_day;
ashg3_mnth=avg_day_to_mnth(mnth_days, ashg3_d, num_data_days);

% % computes the average LAI for the vegetation 
% % (average of min and max)
% % averaged over all pfts in cell
% len=size(canopy_heat3,1);
% ht = canopy_heat3;
% w=canopy_wat3;
% tp3 = ht./( ht_cap + w.*shwv );

tcanopy3 = fort3(:,103);
i=1; j=1; jcount=0; tcanopy3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    tcanopy3_d(j)=tcanopy3_d(j)+tcanopy3(i);
    i = i+1;
end
tcanopy3_d = tcanopy3_d/tstep_per_day;
tcanopy3_mnth=avg_day_to_mnth(mnth_days, tcanopy3_d, num_data_days);

IPARtot3 = fort3(:,101);
i=1; j=1; jcount=0; IPARtot3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    IPARtot3_d(j)=IPARtot3_d(j)+IPARtot3(i);
    i = i+1;
end
IPARtot3_d = IPARtot3_d/tstep_per_day;
IPARtot3_mnth=avg_day_to_mnth(mnth_days, IPARtot3_d, num_data_days);


IPARdir3 = fort3(:,102);
i=1; j=1; jcount=0; IPARdir3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    IPARdir3_d(j)=IPARdir3_d(j)+IPARdir3(i);
    i = i+1;
end
IPARdir3_d = IPARdir3_d/tstep_per_day;
IPARdir3_mnth=avg_day_to_mnth(mnth_days, IPARdir3_d, num_data_days);

% Snow - water equivalent depth [m]
wsn_vtot3  = fort3(:,71)+fort3(:,72)+fort3(:,73);
i=1; j=1; jcount=0; wsn_vtot3_d=zeros(totnum_days,1);
while (i<=size(t,1))
    if (jcount<tstep_per_day)
        jcount = jcount+1;
    else
        j=j+1; jcount=1;
    end
    wsn_vtot3_d(j)=wsn_vtot3_d(j)+wsn_vtot3(i);
    i = i+1;
end
wsn_vtot3_d = wsn_vtot3_d/tstep_per_day;
wsn_vtot3_mnth=avg_day_to_mnth(mnth_days, wsn_vtot3_d, num_data_days);
%--------------------------------------------------------------------------
clear fort3;
%--------------------------------------------------------------------------



plot_all = 1;

if (plot_all == 0)
    t_subday = (t_subday-1983)*365;
    t_subday_data= (DoY+(Hour+5)/24-1);
    figure(1)
    plot(t_subday,IPARtot,'k',...
        t_subday,IPARtot2,'b--',...
        t_subday,IPARtot3,'g-.',...
        t_subday_data,SWdown,'--r+');
    xlabel('t ','FontSize',9)
    ylabel(' Flux [W/m^2]','FontSize',9)
    legend('model1','model2','model3','Measured')
    title('')
    
    figure(2)
    plot(t_subday,srheat,t_subday,srheat2,t_subday_data,SWdown,'r-*')
    xlabel('t ','FontSize',9)
    ylabel(' Incoming SW raidation [W/m^2]','FontSize',9)
    legend('3point','linear','Measured')
    axis([0 5 0 500])
    title('')
    
    figure(3)
    plot(t_subday,ts,t_subday,ts2,t_subday_data,Tair-273.15,'r-*')
    xlabel('t ','FontSize',9)
    ylabel(' Temperature [deg C]','FontSize',9)
    legend('3point','linear','Measured')
    axis([0 5 -20 10])
    title('')
    
elseif(plot_all==1)
    %--------------------------------------------------------------------------
    % PLOTS OF GISSCLIM_LSM DIAGNOSTICS IN THE FORMAT OF KIANG FOR
    % ENT_STANDALONE: FORCINGS
    figure(1)
    subplot(3,2,1)
    title('J0')

    subplot(3,2,2)
    title('J1')

    subplot(3,2,3)
    plot(t_subday,precip,'b',...,
        t_subday,precip2,'k--')
    %axis([0 366 0.01 0.02])
    xlabel('t [yr]','FontSize',9);
    ylabel('Total Precip [mm/d]','FontSize',9);
    legend('constant','linear')

    subplot(3,2,4)
        title('TairC')
    plot(t_subday,ts,'b',...
        t_subday,ts2,'k--')
    %axis([0 400 -0.1 0.3])
    xlabel('t [yr]','FontSize',9)
    ylabel('TairC [deg C]','FontSize',9)
    legend('constant','linear')


    subplot(3,2,5)
    title('TcanopyC')
    plot(t_subday,tcanopy,'b',...
        t_subday,tcanopy2,'k--')
    %axis([0 400 -0.1 0.3])
    xlabel('t [yr]','FontSize',9)
    ylabel('Canopy temperature [deg C]','FontSize',9)
    legend('constant','linear')
    

    subplot(3,2,6)
    title('Qf')
    plot(t_subday,Qf,'b',...
        t_subday,Qf2,'k--')
    %axis([0 366 0.01 0.02])
    xlabel('t [yr]','FontSize',9)
    ylabel('Canopy air humidity, Q_f [kg/kg]','FontSize',9)
    legend('constant','linear')
    

    %--------------------------------------------------------------------------

    figure(2)
    subplot(3,2,1)
     title('P_mbar')
    plot(t_subday,ps,'b',...
        t_subday,ps2,'k--')
    %axis([0 366 0.01 0.02])
    xlabel('t [yr]','FontSize',9)
    ylabel('Pressure [mb]','FontSize',9)
    legend('constant','linear')
   

    subplot(3,2,2)
    title('Ca')


    subplot(3,2,3)
    title('Ch')
    plot(t_subday,cdh,'b',...
        t_subday,cdh2,'k--')
    %axis([0 366 0.01 0.02])
    xlabel('t [day]','FontSize',9)
    ylabel('Atmos. Transfer Coef. [-]','FontSize',9)
    legend('constant','linear')
    

    subplot(3,2,4)
    title('U')
    plot(t_subday,ws,'b',...
         t_subday,ws2,'k--')
    %axis([0 366 0.01 0.02])
    xlabel('t [day]','FontSize',9)
    ylabel('Wind Speed [m/s]','FontSize',9)
    legend('constant','linear')
    

    subplot(3,2,5)
    title('IPARdif')
    plot(t_subday,IPARtot-IPARdir,'b',...
        t_subday,IPARtot2-IPARdir2,'k--')
    %axis([0 366 0.01 0.02])
    xlabel('t [yr]','FontSize',9)
    ylabel('Diffuse IPAR [W/m^2]','FontSize',9)
    legend('constant','linear')
    

    subplot(3,2,6)
    title('IPARdir')
    plot(t_subday,IPARdir,'b',...
        t_subday,IPARdir2,'k--')
    %axis([0 366 0.01 0.02])
    xlabel('t [yr]','FontSize',9)
    ylabel('Direct IPAR [W/m^2]','FontSize',9)
    legend('constant','linear')

    %--------------------------------------------------------------------------

    figure(3)
    subplot(3,2,1)
    title('Sm1')
    plot(t_day,theta_v1_d,'k',...
        t_day,theta2_v1_d,'b--')
    xlabel('t [day]','FontSize',9)
    ylabel('Vol. Soil Moisture, 0 to 0.1 m [-]','FontSize',9)
    legend('constant','linear')
    

    subplot(3,2,2)
    plot(t_day,theta_v2_d,'k',...
        t_day,theta2_v2_d,'b--')
    xlabel('t [day]','FontSize',9)
    ylabel('Vol. Soil Moisture, 0.1 to 0.27 [-]','FontSize',9)
    legend('constant','linear')
    title('Sm2')

    subplot(3,2,3)
    plot(t_day,theta_v3_d,'k',...
        t_day,theta2_v3_d,'b--')
    xlabel('t [day]','FontSize',9)
    ylabel('Vol. Soil Moisture, 0.27 to 0.57 m [-]','FontSize',9)
    legend('constant','linear')
    title('Sm3')

    subplot(3,2,4)
    title('Sm4')
    plot(t_day,theta_v4_d,'k',...
        t_day,theta2_v4_d,'b--')
    xlabel('t [day]','FontSize',9)
    ylabel('Vol. Soil Moisture, 0.57 to 1.1 m [-]','FontSize',9)
    legend('constant','linear')
    

    subplot(3,2,5)
    title('Sm5')
    plot(t_day,theta_v5_d,'k',...
        t_day,theta2_v5_d,'b--')
    xlabel('t [day]','FontSize',9)
    ylabel('Vol. Soil Moisture, 1.1 to 2.0 m [-]','FontSize',9)
    legend('constant','linear')
    

    subplot(3,2,6)
    title('Sm6')
    plot(t_day,theta_v6_d,'k',...
        t_day,theta2_v6_d,'b--')
    xlabel('t [day]','FontSize',9)
    ylabel('Vol. Soil Moisture, 2.0 to 3.5 m [-]','FontSize',9)
    legend('constant','linear')
    
    %--------------------------------------------------------------------------

%     figure(4)
%     subplot(3,2,1)
%     plot(t_day,tsoil_lay1_d,'k',...
%         t_day,tsoil2_lay1_d,'b--',...
%         t_day,tsoil3_lay1_d,'g-.')
%     xlabel('t [yr]','FontSize',9)
%     ylabel('Soil temperature [deg C], 0 to 0.1 m [-]','FontSize',9)
%     legend('no soilevap','soil evap','grid scale')
%     title('St1')
% 
%     subplot(3,2,2)
%     plot(t_day,tsoil_lay2_d,'k',...
%         t_day,tsoil2_lay2_d,'b--',...
%         t_day,tsoil3_lay2_d,'g-.',t_obs)
%     xlabel('t [yr]','FontSize',9)
%     ylabel('Soil temperature [deg C], 0.1 to 0.27 [-]','FontSize',9)
%     legend('no soilevap','soil evap','grid scale')
%     title('St2')
% 
%     subplot(3,2,3)
%     plot(t_day,tsoil_lay3_d,'k',...
%         t_day,tsoil2_lay3_d,'b--',...
%         t_day,tsoil3_lay3_d,'g-.')
%     xlabel('t [yr]','FontSize',9)
%     ylabel('Soil temperature [deg C], 0.27 to 0.57 m [-]','FontSize',9)
%     legend('no soilevap','soil evap','grid scale')
%     title('St3')
% 
%     subplot(3,2,4)
%     plot(t_day,tsoil_lay4_d,'k',...
%         t_day,tsoil2_lay4_d,'b--',...
%         t_day,tsoil3_lay4_d,'g-.')
%     xlabel('t [yr]','FontSize',9)
%     ylabel('Vol. Soil Moisture, 0.57 to 1.1 m [-]','FontSize',9)
%     legend('no soilevap','soil evap','grid scale')
%     title('St4')
% 
%     subplot(3,2,5)
%     plot(t_day,tsoil_lay5_d,'k',...
%         t_day,tsoil2_lay5_d,'b--',...
%         t_day,tsoil3_lay5_d,'g-.')
%     xlabel('t [yr]','FontSize',9)
%     ylabel('Vol. Soil Moisture, 1.1 to 2.0 m [-]','FontSize',9)
%     legend('no soilevap','soil evap','grid scale')
%     title('St5')
% 
%     subplot(3,2,6)
%     plot(t_day,tsoil_lay6_d,'k',...
%         t_day,tsoil2_lay6_d,'b--',...
%         t_day,tsoil3_lay6_d,'g-.')
%     xlabel('t [yr]','FontSize',9)
%     ylabel('Vol. Soil Moisture, 2.0 to 3.5 m [-]','FontSize',9)
%     legend('no soilevap','soil evap','grid scale')
%     title('St6')
    %--------------------------------------------------------------------------

    figure(5)
    subplot(3,2,1)
    title('CosZen')

    subplot(3,2,2)


    subplot(3,2,3)
    plot(t_subday,lai,'k',...
        t_subday,lai2,'b--')
    xlabel('t [yr]','FontSize',9)
    ylabel('LAI [-]','FontSize',9)
    legend('constant','linear')
    title('LAI')

    subplot(3,2,3)

    subplot(3,2,4)

    subplot(3,2,5)

    subplot(3,2,6)

    
    figure(24)
    subplot(3,2,1)
    plot(t_day,agpp2_d,'k',...
        t_day,arauto2_d,'b--',...
        t_day,asoilresp2_d,'g-.',...
        t_day,anpp2_d,'r--')
    xlabel('t [yr]','FontSize',9)
    ylabel('mol C/m^2/day','FontSize',9)
    legend('GPP','R_{auto}','R_{soil}','NPP')
    title('C FLUXES - linear interp')

    subplot(3,2,2)
    title('GPP')
    plot(t_day,agpp_d,'k',...
        t_day,agpp2_d,'b--')
    xlabel('t [yr]','FontSize',9)
    ylabel('GPP [mol C/m^2/day]','FontSize',9)
    legend('constant','linear')
    
    subplot(3,2,3)
    title('Rauto')
    plot(t_day,arauto_d,'k',...
        t_day,arauto2_d,'b--')
    xlabel('t [yr]','FontSize',9)
    ylabel('Autotropic respiration [mol C/m^2/day]','FontSize',9)
    legend('constant','linear')

    subplot(3,2,4)
    title('Soilresp')
    plot(t_day,asoilresp_d,'k',...
        t_day,asoilresp2_d,'b--')
    xlabel('t [yr]','FontSize',9)
    ylabel('Soil respiration [mol C/m^2/day]','FontSize',9)
    legend('constant','linear')


    subplot(3,2,5)
    plot(t_day,anpp_d,'k',...
        t_day,anpp2_d,'b--')
    xlabel('t [yr]','FontSize',9)
    ylabel('NPP [mol C/m^2/day]','FontSize',9)
    legend('constant','linear')
    title('NPP')

    subplot(3,2,6)
    title('CO2flux')
    plot(t_day,anee_d,'k',...
        t_day,anee2_d,'b--')
    xlabel('t [yr]','FontSize',9)
    ylabel('NEE [mol C/m^2/day]','FontSize',9)
    legend('constant','linear')

else

    %--------------------------------------------------------------------------
    % Input meteorological drivers variables
    figure(1)
    subplot(3,3,1);
    plot(t_subday,precip,t_day,precip_d);
    %axis([0 366 0.01 0.02])
    xlabel('t [day]','FontSize',9);
    ylabel('Total Precip [mm/d]','FontSize',9);
    legend('half-hourly','daily')

    subplot(3,3,2);
    plot(t_subday,srheat,t_day,srheat_d)
    %axis([0 366 0.01 0.02])
    xlabel('t [day]','FontSize',9)
    ylabel('Incoming SW [W/m^2]','FontSize',9)
    legend('half-hourly','daily')

    subplot(3,3,3);
    plot(t_subday,trheat,t_day,trheat_d)
    %axis([0 366 0.01 0.02])
    xlabel('t [d]','FontSize',9)
    ylabel('Incoming LW [W/m^2]','FontSize',9)
    legend('half-hourly','daily')

    subplot(3,3,4);
    plot(t_subday,ts,t_day,ts_d)
    %axis([0 400 -0.1 0.3])
    xlabel('t [day]','FontSize',9)
    ylabel('Air Temp [deg C]','FontSize',9)
    legend('half-hourly','daily')

    subplot(3,3,5);
    plot(t_subday,qs,t_day,qs_d)
    %axis([0 366 0.01 0.02])
    xlabel('t [day]','FontSize',9)
    ylabel('Humidity ratio [kg/kg]','FontSize',9)
    legend('half-hourly','daily')

    subplot(3,3,6);
    plot(t_subday,ps,t_day,ps_d)
    %axis([0 366 0.01 0.02])
    xlabel('t [day]','FontSize',9)
    ylabel('Pressure [mb]','FontSize',9)
    legend('half-hourly','daily')

    subplot(3,3,7);
    plot(t_subday,rhosurf,t_day,rhosurf_d)
    %axis([0 366 0.01 0.02])
    xlabel('t [day]','FontSize',9)
    ylabel('Air density [kg/m^3]','FontSize',9)
    legend('half-hourly','daily')

    subplot(3,3,8);
    plot(t_subday,cdh,t_day,cdh_d)
    %axis([0 366 0.01 0.02])
    xlabel('t [day]','FontSize',9)
    ylabel('Atmos. Transfer Coef. [-]','FontSize',9)
    legend('half-hourly','daily')

    subplot(3,3,9);
    plot(t_subday,ws,t_day,ws_d)
    %axis([0 366 0.01 0.02])
    xlabel('t [day]','FontSize',9)
    ylabel('Wind Speed [m/s]','FontSize',9)
    legend('half-hourly','daily')


    figure(2)
    subplot(1,3,1)
    plot(t_subday,GCANOPY,'b+',...
        t_day,GCANOPY_d,'k--',...
        plot_yr+t_mon/12,GCANOPY_mnth,'g')
    xlabel('t [d]','FontSize',9)
    ylabel('Canopy conductance [m/s]','FontSize',9)
    legend('hourly','daily','monthly')

    subplot(1,3,2)
    plot(t_subday,Qf,'b+',...
        t_day,Qf_d,'k--',...
        plot_yr+t_mon/12,Qf_mnth,'g')
    xlabel('t [d]','FontSize',9)
    ylabel('Qf [kg/kg]','FontSize',9)
    legend('hourly','daily','monthly')

    subplot(1,3,3)
    plot(t_subday,agpp,'b+',...
        t_day,agpp_d,'k--',...
        plot_yr+t_mon/12,agpp_mnth,'g')
    xlabel('t [d]','FontSize',9)
    ylabel('GPP [mol C/m^2/day]','FontSize',9)
    legend('hourly','daily','monthly')

    figure(3)
    subplot(1,3,1)
    plot(t_subday,cdh.*ws,'b+',...
        t_day,cdh_d.*ws_d,'k--',...
        plot_yr+t_mon/12,cdh_mnth.*ws_mnth,'g')
    xlabel('time [yr]','FontSize',9)
    ylabel('C_h * U [m/s]','FontSize',9)
    legend('30 min','daily','monthly')

    subplot(1,3,2)
    plot(t_subday,ws,'b+',...
        t_day,ws_d,'k--',...
        plot_yr+t_mon/12,ws_mnth,'g')
    xlabel('time [yr]','FontSize',9)
    ylabel('U [m/s]','FontSize',9)
    legend('30 min','daily','monthly')

    subplot(1,3,3)
    plot(t_subday,cdh,'b+',...
        t_day,cdh_d,'k--',...
        plot_yr+t_mon/12,cdh_mnth,'g')
    xlabel('time [yr]','FontSize',9)
    ylabel('C_h [-]','FontSize',9)
    legend('30 min','daily','monthly')



    figure(4)
    subplot(3,3,1);
    plot(t_subday,IPARtot,'b+',...
        t_day,IPARtot_d,'k--',...
        t_mon/12+plot_yr,IPARtot_mnth,'g-.')
    %axis([0 366 0.01 0.02])
    xlabel('t [yr]','FontSize',9)
    ylabel('Total IPAR [W/m^2]','FontSize',9)
    legend('half-hourly','daily','monthly')

    subplot(3,3,2);
    plot(t_subday,IPARdir,'b+',...
        t_day,IPARdir_d,'k--',...
        t_mon/12+plot_yr,IPARdir_mnth,'g-.')
    %axis([0 366 0.01 0.02])
    xlabel('t [yr]','FontSize',9)
    ylabel('Direct IPAR [W/m^2]','FontSize',9)
    legend('half-hourly','daily','monthly')

    subplot(3,3,3);
    plot(t_subday,IPARtot-IPARdir,'b+',...
        t_day,IPARtot_d-IPARdir_d,'k--',...
        t_mon/12+plot_yr,IPARtot_mnth-IPARdir_mnth,'g-.')
    %axis([0 366 0.01 0.02])
    xlabel('t [yr]','FontSize',9)
    ylabel('Diffuse IPAR [W/m^2]','FontSize',9)
    legend('half-hourly','daily','monthly')

    subplot(3,3,4);
    plot(t_subday,ts,'b+',...
        t_day,ts_d,'k--',...
        t_mon/12+plot_yr,ts_mnth,'g-.')
    %axis([0 400 -0.1 0.3])
    xlabel('t [yr]','FontSize',9)
    ylabel('Air Temp [deg C]','FontSize',9)
    legend('half-hourly','daily','monthly')

    subplot(3,3,5);
    plot(t_subday,tcanopy,'b+',...
        t_day,tcanopy_d,'k--',...
        t_mon/12+plot_yr,tcanopy_mnth,'g-.')
    %axis([0 400 -0.1 0.3])
    xlabel('t [yr]','FontSize',9)
    ylabel('Canopy temperature [deg C]','FontSize',9)
    legend('half-hourly','daily','monthly')

    subplot(3,3,6);
    plot(t_subday,Qf,'b+',...
        t_day,Qf_d,'k--',...
        t_mon/12+plot_yr,Qf_mnth,'g-.')
    %axis([0 366 0.01 0.02])
    xlabel('t [yr]','FontSize',9)
    ylabel('Canopy air humidity, Q_f [kg/kg]','FontSize',9)
    legend('half-hourly','daily','monthly')

    subplot(3,3,7);
    plot(t_subday,ps,'b+',...
        t_day,ps_d,'k--',...
        t_mon/12+plot_yr,ps_mnth,'g-.')
    %axis([0 366 0.01 0.02])
    xlabel('t [yr]','FontSize',9)
    ylabel('Pressure [mb]','FontSize',9)
    legend('half-hourly','daily','monthly')

    subplot(3,3,8);
    plot(t_subday,cdh,'b+',...
        t_day,cdh_d,'k--',...
        t_mon/12+plot_yr,cdh_mnth,'g-.')
    %axis([0 366 0.01 0.02])
    xlabel('t [yr]','FontSize',9)
    ylabel('Atmos. Transfer Coef. [-]','FontSize',9)
    legend('half-hourly','daily','monthly')

    subplot(3,3,9);
    plot(t_subday,ws,'b+',...
        t_day,ws_d,'k--',...
        t_mon/12+plot_yr,ws_mnth,'g-.')
    %axis([0 366 0.01 0.02])
    xlabel('t [yr]','FontSize',9)
    ylabel('Wind Speed [m/s]','FontSize',9)
    legend('half-hourly','daily','monthly')


    figure
    subplot(3,1,1)
    plot(t_hr, NEE_data, 'r.', ...
        t_hr, NEE_miss,'g+')
    xlabel('t [yr]','FontSize',9)
    ylabel('NEE [mol C/m^2/day]','FontSize',9)
    legend('data','missing data')
    title('FLUXNET DATA')

    subplot(3,1,2)
    plot(t_hr, LE_Wm2 , 'r.', ...
        t_hr, LE_Wm2_miss,'g+')
    xlabel('t [yr]','FontSize',9)
    ylabel('LE [W/m^2]','FontSize',9)
    legend('data','missing data')

    subplot(3,1,3)
    plot(t_hr, H_Wm2, 'r.', ...
        t_hr, H_Wm2_miss,'g+')
    xlabel('t [yr]','FontSize',9)
    ylabel('H [W/m^2]','FontSize',9)
    legend('data','missing data')

end