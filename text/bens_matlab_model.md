# The original model code

Ben kindly supplied me with the Matlab code used by them for their 2012 sea level rise paper. Hereafter I'll try to understand the code with close attention to the scaling parameters and time scales.

```matlab
function model_cru401(region_number)

L_V_scaling_error=1;      % in percent/100 (i.e., 1 = 100 %)
A_V_scaling_error=0.4;      % in percent/100 (i.e., 0.4 = 40 %)

area_error_rgi=5;  % in percent

tau_L_error=5;            % in percent/100 (i.e., 5 = 500 %)
tau_A_error=5;            % in percent/100 (i.e., 5 = 500 %)
```

At the start of the `model_cru401` function some relative error thresholds are set.

```matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% read in rgi data

eval(['load rgi_v6/region' region_number '.mat;']);

nan_idx=zeros(size(rgi_area));

nan_idx(isnan(rgi_lat) | isnan(rgi_lon) | isnan(rgi_year) | isnan (rgi_zmax) | isnan(rgi_zmean) | isnan(rgi_zmin))=1;

rgi_year(rgi_year>2012)=2012;
rgi_year(rgi_year<1903)=1903;

%%% end read in rgi data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
```

Each RGI region has it's own input routine in a separate `*.mat` file, which is called here given the `region_number`. @TODO: eventually look into these file(s).

```matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% start model

load cross_validation_cru401 t_star model_bias number_closest_glaciers length_of_climatology lon_mb lat_mb cru_years lapse_precip temp_prec_solid temp_melt mse_mu_star_cross

cru_years=1901:2016;
```

Reading a lot of other data files?! @TODO: I'll probably have to run the Matlab code to understand...

```matlab
% Bahr et al. 1997, The physical basis...
gamma_glacier=1.375;
% c_a_glacier=0.191 m^(3-2gamma), Bahr 1997, Global distributions of glacier properties: A stochastic scaling paradigm
c_a_glacier=0.0340;
% Bahr et al. 1997, The physical basis...
q_glacier=2.2;
% c_l_glacier=4.55 m^(3-q), Radic, Hock, Oerlemans 2008: Analysis of scaling methods 
c_l_glacier=0.018;

% Radic & Hock (2010), Regional and global volumes of glaciers derived from...
gamma_icecap=1.25;
% =1.7026 m^(3-2gamma), Radic & Hock (2010), Regional and global volumes of glaciers derived from...
c_a_icecap=0.0538;
% from H = 3.4/(1000^0.5) L^0.5, 3.4 (Patterson 1994) for L and H in m changed to km...
q_icecap=2.5;
% and V = 2/3 pi H L^2 
c_l_icecap=0.2252;

```

**Important:** This is a key section, since the scaling parameters (in their specific units) are defined here. As far as I can tell, my considerations were right and the change of units results in a change of the scaling constants (`c_a` and `c_l`). The parameters for ice caps are not of interest, at least for now...

```matlab
mean_rmse=mean(mse_mu_star_cross.^(1/2));

time_cru=(1901+(1/24)):(1/12):(2017-(1/24));

height_cru_glacier=nan(1,length(rgi_area));
height_cru_clim_glacier=nan(1,length(rgi_area));
lapse_temp=nan(1,length(rgi_area));
rho_lapse_temp=nan(1,length(rgi_area));
p_lapse_temp=nan(1,length(rgi_area));
dist=nan(length(rgi_area),length(lon_mb));
mean_dist=nan(1,length(rgi_area));
t_star_rgi=nan(1,length(rgi_area));
model_bias_rgi=nan(1,length(rgi_area));
mu_rgi=nan(1,length(rgi_area));
temp_range=nan(1,length(rgi_area));
turnover=nan(1,length(rgi_area));
temp_glacier_tongue_clim=nan(1,12);
precip_glacier_clim=nan(1,12);
A=nan(length(rgi_area),length(time_cru)/12);
L=nan(length(rgi_area),length(time_cru)/12);
V=nan(length(rgi_area),length(time_cru)/12);
dL=nan(length(rgi_area),length(time_cru)/12);
dV=nan(length(rgi_area),length(time_cru)/12-1);
mb_modeled=nan(length(rgi_area),length(time_cru)/12-1);
mb_monthly_sust=nan(length(rgi_area),length(time_cru)/12-1,12);
mb_monthly_total=nan(length(rgi_area),length(time_cru)/12-1,12);
n_months_above_freezing=nan(length(rgi_area),length(time_cru)/12);
A_error=nan(length(rgi_area),length(time_cru)/12);
L_error=nan(length(rgi_area),length(time_cru)/12);
V_error=nan(length(rgi_area),length(time_cru)/12);
dL_cum_error=nan(length(rgi_area),length(time_cru)/12);
mb_error=nan(length(rgi_area),length(time_cru)/12);
dV_cum_error=nan(length(rgi_area),length(time_cru)/12-1);
```

Well...

