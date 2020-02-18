import argparse
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import cmocean


def convert_pr_units(darray):
    """Convert kg m-2 s-1 to mm day-1.
    
    Args:
      darray (xarray.DataArray): Precipitation data
    
    """
    
    darray.data = darray.data*86400
    darray.attrs['units'] = 'mm/day'
    
    return darray


def create_plot(clim, model_name, season, gridlines=False, levels=None):
    """Plot the precipitation climatology.
    
    Args:
      clim (xarray.DataArray): Precipitation climatology data
      model_name (str): Name of the climate model
      season (str): Season
      
    Kwargs:
      gridlines (bool): Select whether to plot gridlines    
    
    """

    if not levels:
        levels = np.arange(0, 13.5, 1.5)    
    
    fig = plt.figure(figsize=[12,5])
    ax = fig.add_subplot(111, projection=ccrs.PlateCarree(central_longitude=180))
    clim.sel(season=season).plot.contourf(ax=ax,
                                          levels=levels,
                                          extend='max',
                                          transform=ccrs.PlateCarree(),
                                          cbar_kwargs={'label': clim.units},
                                          cmap=cmocean.cm.haline_r)
    ax.coastlines()
    if gridlines:
        plt.gca().gridlines()

    title = '%s precipitation climatology (%s)' %(model_name, season)
    plt.title(title)


def main(inargs):
    """Run the program."""

    dset = xr.open_dataset(inargs.pr_file)
    
    clim = dset['pr'].groupby('time.season').mean('time', keep_attrs=True)
    clim = convert_pr_units(clim)

    create_plot(clim, dset.attrs['model_id'], inargs.season, gridlines=inargs.gridlines, levels=inargs.cbar_levels)
    plt.savefig(inargs.output_file, dpi=200)


if __name__ == '__main__':
#说明：在命令行下依次输入pr_file，season，output_file，并将这三个字符串整合到一个变量args中，使整个脚本显得更加干净。
    description='Plot the precipitation climatology for a given season.'
    parser = argparse.ArgumentParser(description=description)
    
    parser.add_argument("pr_file", type=str, help="Precipitation data file")
    parser.add_argument("season", type=str, choices=['DJF', 'MAM', 'JJA', 'SON'],help="Season to plot")
    parser.add_argument("output_file", type=str, help="Output file name")
    parser.add_argument("--gridlines", action="store_true", default=False,help="Include gridlines on the plot")#store_true表示当在终端输入时不输入这项，则为默认值，当终端输入时输入这项，则更改为True
    parser.add_argument("--cbar_levels", type=float, nargs='*', default=None, help='list of levels / tick marks to appear on the colorbar')
    args = parser.parse_args()
    
    main(args)
    
