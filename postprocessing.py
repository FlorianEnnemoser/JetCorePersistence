import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cartopy.crs as ccrs
import cartopy.feature as cf
import xarray as xr
import pathlib
from preprocessing import Region

def plotting_simple(dataarray : xr.DataArray,
                    levels : int = 11,
                    cmap : str ="jet",
                    title : str ="",
                    cbartitle : str ="",
                    mark_area : Region = None,
                    proj = ccrs.PlateCarree()
                    ):
    
    fig = plt.figure(dpi=120)
    ax = plt.axes(projection=proj)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=0.2, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'fontfamily': 'serif'}
    gl.ylabel_style = {'fontfamily': 'serif'}
    ax.set_xlabel("Longitude (°E)",fontfamily='serif')
    ax.set_ylabel("Latitude (°N)",fontfamily='serif')
    #ax.set_extent([-180, 180, 0, -90], crs=ccrs.PlateCarree()) #lon1,lon2,lat1,lat2
    im = dataarray.plot.pcolormesh(ax=ax,
                            cmap=cmap,
                            levels=levels,
                            transform=ccrs.PlateCarree(),add_colorbar=False)        
                        #cbar_kwargs={'orientation': 'horizontal','label':cbartitle})
    
    ax.coastlines(resolution='50m',color="k",lw=0.5)

    if mark_area:
        ax.plot(mark_area.lon_rect, mark_area.lat_rect, color='red', linewidth=2, transform=ccrs.PlateCarree())        

    cb = fig.colorbar(im, ax=ax,
                              orientation="horizontal")
    cb.set_label(label=cbartitle,fontsize=12,fontfamily='serif')
    plt.title(title)
    for ax in fig.get_axes():
        for xlabel in ax.get_xticklabels():
            xlabel.set_fontfamily('serif')
        for ylabel in ax.get_yticklabels():
            ylabel.set_fontfamily('serif')
    return fig