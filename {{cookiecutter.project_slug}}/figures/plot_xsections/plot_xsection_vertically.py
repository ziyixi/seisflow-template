"""
Plot xsection for the PPM model. (netcdf file) The model variables should be an array with the axis order of (lon,lat,dep).
"""
import numpy as np 
from scipy.io.netcdf import netcdf_file
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
import tqdm
from numpy.linalg import norm
from obspy.taup import TauPyModel

model = TauPyModel(model="iasp91")

def load_netcdf(netcdf_fname,variable):
    """
    Load netcdf files
    """
    result=None
    with netcdf_file(netcdf_fname,"r") as f:
        result=f.variables[variable][:].copy()
        longitude=f.variables["longitude"][:].copy()
        latitude=f.variables["latitude"][:].copy()
        depth=f.variables["depth"][:].copy()
    return result,longitude,latitude,depth

def interpolate_grids(grid_data,longitude,latitude,depth):
    """
    Interpolate the given grids, return the interpolated function.
    """
    interp_func=RegularGridInterpolator((longitude,latitude,depth),grid_data,method="linear")
    return interp_func


def generate_figure_mesh(start_point,end_point,depth_range,resolution,interp_func,used_coorname):
    """
    Generate the mesh points to interpolate.
    """
    start_lon,start_lat=start_point
    end_lon,end_lat=end_point
    start_dep,end_dep=depth_range
    h_res,v_res=resolution
    # we don't consider the geological feature, just generate some even grids.
    lons=np.linspace(start_lon,end_lon,h_res+1)
    lats=np.linspace(start_lat,end_lat,h_res+1)
    deps=np.linspace(start_dep,end_dep,v_res+1)
    array_to_interpolate=np.zeros((h_res+1,v_res+1,3))
    for ih in range(h_res+1):
        for iv in range(v_res+1):
            array_to_interpolate[ih,iv,:]=np.array([lons[ih],lats[ih],deps[iv]])
    Z=interp_func(array_to_interpolate)
    if(used_coorname=="latitude"):
        x=lats
    elif(used_coorname=="longitude"):
        x=lons
    else:
        raise Exception("no such coordinate name")
    y=deps
    X,Y=np.meshgrid(x,y,indexing="ij")
    return X,Y,Z

def plot_figure(X,Y,Z,auto_cmap=True,vmin=None,vmax=None,level=100,vround=2,plot_type="perturbation",used_coorname="latitude"):
    """
    Plot the contourf map.
    """
    # ! here we use %
    vround=0
    Z=Z*100
    if(auto_cmap):
        vmin_round = round(np.min(Z), 2)
        if(vmin_round < np.min(Z)):
            vmin = vmin_round
        else:
            vmin = vmin_round-10**(-vround)
        vmax_round = round(np.max(Z), 2)
        if(vmax_round > np.max(Z)):
            vmax = vmax_round
        else:
            vmax = vmax_round+10**(-vround)
    else:
        # make all the values above vmax and below vmin as vmax and vmin
        Z[Z>vmax]=vmax
        Z[Z<vmin]=vmin
    v = np.arange(vmin, vmax+10**(-vround), 10**(-vround))
    plt.figure(figsize=(10,4))
    levels = np.linspace(vmin, vmax, level+1)
    # ! here we use %
    plt.contourf(X, Y, Z,levels=levels, cmap=plt.cm.seismic_r)
    cbar=plt.colorbar(ticks=v, label=plot_type,extend="both",orientation="horizontal",fraction=0.16,pad=0.19)
    if(plot_type=="perturbation"):
        # put for temporarily useage
        cbar.set_label(r'$\delta lnV_s$(%)',size=20)
    plt.xlabel(used_coorname+" (°)",fontsize=20)
    plt.ylabel("depth (km)",fontsize=20)
    plt.gca().invert_yaxis()
    ax=plt.gca()
    # * if plot earthquakes
    ax.tick_params(axis='both', which='major', labelsize=15)


def plot_earthquakes(tolerance,event_list,start_point,end_point,used_coorname,depth_range):
    """
    Plot earthquakes near the cross section. 
    """
    events=np.loadtxt(event_list)
    p1=np.array(start_point)
    p2=np.array(end_point)
    used_events=[]
    for row in events:
        lat,lon,dep=row
        p3=np.array([lon,lat])
        d = norm(np.cross(p2-p1, p1-p3))/norm(p2-p1)
        if(d<=tolerance):
            minlon=min(start_point[0],end_point[0])
            maxlon=max(start_point[0],end_point[0])
            minlat=min(start_point[1],end_point[1])
            maxlat=max(start_point[1],end_point[1])
            if((minlat<lat<maxlat) & (minlon<lon<maxlon)):
                used_events.append(row)
    all_used_events=np.array(used_events)
    if(used_coorname=="latitude"):
        scatter_x=all_used_events[:,0]
        xlim=[start_point[1],end_point[1]]
    elif(used_coorname=="longitude"):
        scatter_x=all_used_events[:,1]
        xlim=[start_point[0],end_point[0]]
    scatter_y=all_used_events[:,2]
    plt.scatter(scatter_x,scatter_y,color="k",s=4)
    plt.ylim([depth_range[1],depth_range[0]])
    plt.xlim(xlim)

def plot_events_and_stations():
    used_events=[31.61,138.17]
    used_stations=[
        [45.6263, 131.8099],
        [45.7422, 131.0790],
        [46.3296, 132.1489],
        [46.6589, 131.0930],
        [46.7325, 130.3245],
        [47.3528, 130.2364],
        [47.5960, 130.8286],
        [49.0412, 129.8660],
        [49.0813, 128.9159]
    ]
    used_stations=np.array(used_stations)
    stations_scatter_y=np.zeros_like(used_stations[:,0])
    stations_scatter_y[:]=-10
    if(used_coorname=="latitude"):
        scatter_x=used_events[0]
        plt.scatter(used_stations[:,0],stations_scatter_y,marker="^",color="b",s=30,clip_on=False)
    elif(used_coorname=="longitude"):
        scatter_x=used_events[1]
        plt.scatter(used_stations[:,1],stations_scatter_y,marker="^",color="b",s=30,clip_on=False)
    scatter_y=351.95
    plt.scatter(scatter_x,scatter_y,marker="*",color="y",s=160)
    for each_station in used_stations:
        arrivals = model.get_ray_paths_geo(source_depth_in_km=scatter_y, source_latitude_in_deg=used_events[0],
            source_longitude_in_deg=used_events[1],receiver_latitude_in_deg=each_station[0],receiver_longitude_in_deg=each_station[1],phase_list=["S"])
        dist_km=arrivals[0].path["dist"]
        dist_ratio=dist_km/dist_km[-1]
        lat_plot_start=used_events[0]
        lat_plot_end=each_station[0]
        dist_lat=np.zeros_like(dist_ratio)
        for index in range(len(dist_ratio)):
            dist_lat[index]=lat_plot_start+(lat_plot_end-lat_plot_start)*dist_ratio[index]
        dist_depth=arrivals[0].path["depth"]
        plt.plot(dist_lat,dist_depth,color="k",linewidth=.5)
    

if __name__ == "__main__":
    netcdf_fname="./hybrid_noeli_notopo_no410.nc"
    variable="vsv"
    # netcdf_fname="./EARA2014_dlnVs_h50v10_1000km_EAsia_East.nc"
    # variable="vs"
    start_point=(138,31)
    end_point=(132,50)
    depth_range=(0,800)
    resolution=(200,30)
    used_coorname="latitude"
    auto_cmap=False
    # vmin,vmax=-0.07,0.07
    vmin,vmax=-6,6
    level=500
    vround=2
    plot_type="perturbation"
    tolerance=0.5
    event_list="./China_Tomo_seis_19600102_20081231_EHB"
    # code part
    grid_data,longitude,latitude,depth=load_netcdf(netcdf_fname,variable)
    interp_func=interpolate_grids(grid_data,longitude,latitude,depth)
    X,Y,Z=generate_figure_mesh(start_point,end_point,depth_range,resolution,interp_func,used_coorname)
    plot_figure(X,Y,Z,auto_cmap=auto_cmap,vmin=vmin,vmax=vmax,level=level,vround=vround,plot_type=plot_type,used_coorname=used_coorname)
    plot_earthquakes(tolerance,event_list,start_point,end_point,used_coorname,depth_range)
    plot_events_and_stations()
    plt.show()