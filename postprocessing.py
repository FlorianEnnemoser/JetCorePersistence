import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.transforms as mtransforms
import cartopy.crs as ccrs
import cartopy.feature as cf
import numpy as np
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

def create_data_and_plot_map_correlation(linregression,
                                         levels : int =32,
                                         cbar_title : str ="cbar_title (missing)",
                                         var_display : str ="slope",
                                         figurename : str=None,
                                         norm : np.ndarray =None,
                                         cmap : str="RdBu_r"):

    def add_number(ax,label,extratext):
        trans = mtransforms.ScaledTranslation(5/72, 12/72, fig.dpi_scale_trans)
        axs[ax].text(0.0, 1.0, f"{label}) {extratext}", transform=axs[ax].transAxes + trans,
                    fontsize='medium', verticalalignment='top', fontfamily='serif')
    
    proj = ccrs.PlateCarree()
    fig, axs = plt.subplot_mosaic("""
                                ab
                                cd
                              """,figsize=(5*phi,5) ,layout='constrained',sharex=True,sharey=True,
                                  subplot_kw={'projection': proj})

    for ax,s in zip(axs,seasons_idx):
        add_number(ax,ax,s)
        fig_corr(axs[ax],[-20,40,32.5,75],fig)
        if var_display == "r2":
            print("correlation^2")
            im = (linregression.sel(seasons=s)["rvalue"]**2 *100) \
                .plot.contourf(ax=axs[ax],levels=levels,cmap=cmap,norm=norm,add_colorbar=False) #rvalue slope pvalue
        else:
            im = (linregression.sel(seasons=s)[var_display]) \
                .plot.contourf(ax=axs[ax],levels=levels,cmap=cmap,norm=norm,add_colorbar=False) #rvalue slope pvalue
            

        (linregression.sel(seasons=s).pvalue <= 0.1) \
            .plot.contourf(ax=axs[ax],hatches=[None,"....","...."],colors="None",add_colorbar=False,levels=[0,0.5,1])
        axs[ax].set_title("")
        
    cbar_ax = fig.add_axes([0.15, -0.10, 0.7, 0.05]) #left,bot,breite,höhe
    cb = fig.colorbar(im, cax=cbar_ax,
                                  orientation="horizontal",extend="both",
                                  fraction=0.10) 
    cb.set_label(label=cbar_title,fontfamily='serif')
    set_correct_font(fig)
    
    if figurename:
        plt.savefig(savepath(figurename),bbox_inches='tight',facecolor="none", edgecolor='none')


def Figure_One_AvgPers_vs_Lat_and_Boxplot(figure_title : str = "Fig1_AvgPers_vs_Latitude_UnderestimationRegion", persistence_over_lat : xr.DataArray = None):
    plot_path = pathlib.Path(r"")
    h, w = figaspect(16/9)
    fig, axs = plt.subplot_mosaic("""a""", layout='constrained',sharex=True,sharey=True, figsize=(w,h))
    models = select_only_models(persistence_over_lat)["mean"]#.mean("longitude")
    obs = select_only_obs(persistence_over_lat)["mean"]#.mean("longitude")
    models_counts = select_only_models(persistence_over_lat)#.mean("longitude")
    obs_counts = select_only_obs(persistence_over_lat)#.mean("longitude")

    ax = "a"
    models.mean("model").plot(ax=axs[ax],c="b",marker="x")
    obs.plot(ax=axs[ax],c="k",marker="x")
    axs[ax].grid(ls="--",c="grey")
    axs[ax].set_title("")
    axs[ax].set_xlabel("")
    axs[ax].set_ylabel("")
    annotate_lines(axs[ax],obs_counts,weight='bold')
    annotate_lines(axs[ax],models_counts.mean("model"),shiftx=5,shifty=5)
    ax_box = axs[ax].inset_axes([1.15, 0, 0.2, 1], transform=axs[ax].transAxes)
    boxprops=dict(facecolor='b')
    whis=(5, 95)
    medianprops = dict(linestyle='-.', linewidth=2.5, color='orange')
    ax_box.boxplot(models.mean("latitude_bins"), vert=True, widths=[0.2], patch_artist=True, boxprops=boxprops,whis=whis,medianprops=medianprops)    
    ax_box.grid(ls="--",c="grey")
    axs[ax].set_ylim(-0.75, 15.75)
    axs[ax].set_yticks([0,5,10,15])
    ax_box.set_ylim(5,8)
    ax_box.set_yticks([5,6,7,8])
    ax_box.set_yticklabels(ax_box.get_yticklabels(),fontfamily='serif')
    ax_box.set_xticklabels(["MMS"], fontfamily='serif')
    ax_box.axhline(models.mean(["model","latitude_bins"]),c="b",ls="--")
    ax_box.axhline(obs.mean("latitude_bins"),c="k")    
    ax_box.set_title('')            
    fig.supylabel("Mean persistence of jet core (days)", fontfamily='serif')
    fig.supxlabel("Latitude (°N)", fontfamily='serif')
    _set_correct_font(fig)
    plt.savefig(plot_path / f"{figure_title}.pdf", bbox_inches = 'tight',facecolor="none", edgecolor='none',dpi=200)
    plt.close()    