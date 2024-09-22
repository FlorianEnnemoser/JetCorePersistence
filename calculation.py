import numpy as np
import xarray as xr
import scipy as sp

SEASON_MAPPING = {1: "DJF", 2: "DJF", 3: "MAM", 4: "MAM", 5: "MAM", 6: "JJA", 7: "JJA", 8: "JJA", 9: "SON", 10: "SON", 11: "SON", 12: "DJF"}

def find_latitudinal_max(data : xr.DataArray) -> xr.DataArray:
    lat_maxima = data.max(dim='latitude')    
    return data.where(data == lat_maxima)

def jetcore_latitudes(data : xr.DataArray) -> xr.DataArray:
    return find_latitudinal_max(data).idxmax("latitude").rename("jetcore")

def JLI(data : xr.DataArray) -> xr.DataArray:
    return find_latitudinal_max(data).idxmax("latitude").mean("longitude").rename("JLI")

def JLI_Woollings(data : xr.DataArray) -> xr.DataArray:
    return data.mean("longitude").idxmax("latitude").rename("jetcore_woolings")

def JCO(data : xr.DataArray) -> xr.DataArray:
    return find_latitudinal_max(data).count("time").rename("JCO")

def calc_EKE(data_u : xr.DataArray , data_v : xr.DataArray) -> xr.DataArray:
    group_data_u = data_u.groupby("time.dayofyear")
    group_data_v = data_v.groupby("time.dayofyear")
    u_mean = group_data_u.mean('time')
    v_mean = group_data_v.mean('time')
    EKE = ( 0.5 * ((group_data_u-u_mean)**2 + (group_data_v-v_mean)**2) ).rename('EKE')
    EKE = EKE.drop('dayofyear').squeeze(drop=True)
    return EKE


def _new_linregress(x, y) -> np.ndarray:
    return np.array([*sp.stats.linregress(x, y)])

def lin_regress_xr(x : xr.DataArray, y : xr.DataArray, dim : str) -> xr.Dataset:
    """dim -> über die es applied wird. also x, y wird korreliert über dim (zb 'time')"""
    return xr.apply_ufunc(_new_linregress, x, y,
                        input_core_dims=[[dim], [dim]],
                        output_core_dims=[["parameter"]],
                        vectorize=True,
                        dask="parallelized",
                        output_dtypes=['float64'],
                        output_sizes={"parameter": 5},
                        ).assign_coords({"parameter":["k","d","r","p","err"]}).to_dataset(dim="parameter")

def seasonal_mapping(first_day,last_day) -> str:
    event_range = xr.date_range(start=str(first_day.values),end=str(last_day.values),freq="D")
    unique_months, unique_counts = np.unique(event_range.month, return_counts=True)
    return SEASON_MAPPING.get(unique_months[np.argmax(unique_counts)])

def persistence(data: xr.DataArray, lat_limit: float = 2.5, consecutive_days: int = 3) -> xr.DataArray:
    data = data.squeeze()
    counter, first_day = 0, 0
    counter_arr, average_lat, save_average_lat, season_info  = [], [], [], []

    for daily_difference in np.abs(data.diff("time")):
        if daily_difference <= lat_limit:
            if counter == 0:
                first_day = daily_difference.time
            counter += 1
            average_lat.append(data.sel(time=daily_difference.time))
        else:
            if counter <= consecutive_days:
                counter, first_day = 0, 0
                average_lat = []
            else:
                counter_arr.append(counter)
                save_average_lat.append(np.mean(average_lat))
                season_info.append(seasonal_mapping(first_day,daily_difference.time))
                counter, average_lat = 0, []

    return xr.DataArray(counter_arr,coords={"latitude":save_average_lat},name="consecutivedays").assign_coords(seasons=("latitude",season_info))



