% fill up temp and precip data with nans to allow calculation of
% climatologies
precip_filled_nan=[nan(12*(length_of_climatology-1)/2,1);precip_cru_glacier;nan(12*(length_of_climatology-1)/2,1)];

% select data from time window
precip=precip_filled_nan((t_star_rgi_idx-1)*12+1:(t_star_rgi_idx+(length_of_climatology-1))*12);

% calculate climatologies
for k=1:12
    precip_glacier_clim(k)=nanmean(precip(k:12:length(precip)));
end

