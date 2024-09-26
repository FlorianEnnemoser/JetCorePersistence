import parse_config as parse

def main() -> None:
    parse.execute(filepath="configurations/fig1.toml")

def fig1():
    active_region = Region("underestimRegion_30_65Lat_-60_0Lon",lat=(30,65),lon=(-60,0))  
    u = Analysis(save_folder="u_new").load_as_xr("u_new").u
    
    jli = jet(u).JLI()    
    pers_over_lat = []
    pers_raw = []
    for i,g in tqdm(jli.groupby("model"),desc="Calculating Persistence "):
        pers_data = persistence_without_season(g.squeeze()).sortby("latitude")
        pers_raw.append(pers_data.mean())
        pers_over_lat.append(binned_over_latitudes(pers_data))

    pers_over_lat = xr.concat(pers_over_lat,dim=jli.model)

    
    save_path = pathlib.Path(r"F:\MT_nice_python\save_data")    
    pers_over_lat = xr.load_dataset(save_path / f"ALL_{active_region.name}_BinnedPersistence.nc")
    print(pers_over_lat)
    Figure_One_AvgPers_vs_Lat_and_Boxplot(persistence_over_lat=pers_over_lat,figure_title="Fig1_new")        

def figS3():
    jet_pr_pers_wet_data = [correlation_past_jc_pers_pr_wet,
                            correlation_past_jc_lat_pr_wet,
                            correlation_past_jc_velo_pr_wet]
    fignames = ["MAP_HIST_c_p01_jc_pers_pr_wet.svg",
                "MAP_HIST_c_p01_jc_lat_pr_wet.svg",
                "MAP_HIST_c_p01_jc_velo_pr_wet.svg"]

    for d, name in zip(jet_pr_pers_wet_data,fignames):
        plt.close()
        create_data_and_plot_map_correlation(d,
                                        levels=np.linspace(-1,1,11),
                                        var_display="rvalue",
                                        cbar_title=r"Correlation (-)",
                                        figurename=name)
        plt.show()

if __name__ == "__main__":
    main()

