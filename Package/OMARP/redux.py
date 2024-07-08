import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tabulate import tabulate
from astropy import units as u
from astropy.io import fits
from astropy.visualization import ZScaleInterval
from astropy.coordinates import SkyCoord

def main():
    
    Redux("Sky Chart","sloan_u","20230209/Images",use_bias=True,dark_exp_time=60)

def tables(folder):
    """Función para crear un archivo de texto con una tabla con los datos de todos los fits"""
    name_list=sorted(glob.glob(f"{folder}/*.fit"))
    a=np.empty([0,5])
    for i in range(len(name_list)):
        f=fits.open(name_list[i])
        im_type=f[0].header['IMAGETYP']
        im_exptime=f[0].header['EXPTIME']
        im_filter=f[0].header['FILTER']
        if im_type=='Light Frame':
            im_object=f[0].header['OBJECT']
        else:
            im_object='None'
        row=[im_type,im_object,im_exptime,im_filter,name_list[i]]
        a=np.vstack((a,row))
    tab=tabulate(a,headers=["Type","Object","Exp_Time_(s)","Filter","File_Name"],tablefmt='plain') 
    with open('imagelist.txt', 'w') as f:
        f.write(tab)
    return print(tab)    

def Fits_Array(list, extract_exptime=False, x_axis=4096, y_axis=4096):
    """Función para generar un array tridimensional de imágenes, siendo 'list'
    una lista con los pathnames de las imágenes que se quiere estackear"""
    cube=np.empty([x_axis,y_axis,0])
    expt_list=[]
    for i in range(len(list)):
        data=fits.open(list[i])
        cube=np.dstack((cube,np.array(data[0].data)))
        if extract_exptime==True:
            expt_list.append(data[0].header['EXPTIME'])
    if extract_exptime==True:
        return cube,expt_list
    else:
        return cube

def plots(cube,contr,dir,im_name):
    """Función para plotear las imágenes de un cubo"""
    fig=plt.figure(figsize=(20,5*round((cube.shape[2]+1)/4)),layout='constrained')
    for i in range(cube.shape[2]):
        zscale=ZScaleInterval(contrast=contr)
        minv,maxv=zscale.get_limits(cube[:,:,i])
        fig.add_subplot(round((cube.shape[2]+1)/4),4,i+1).imshow(cube[:,:,i],vmin=minv,vmax=maxv,cmap='gray',origin='lower')
        plt.axis('off')
    plt.suptitle(f'{im_name}')
    plt.savefig(f"{dir}/{im_name}.png",facecolor='white')
    plt.close()

def s_plot(im,contr,dir,im_name,size=20):
    fig=plt.figure(figsize=(size,size))
    zscale=ZScaleInterval(contrast=contr)
    minv,maxv=zscale.get_limits(im)
    plt.imshow(im,vmin=minv,vmax=maxv,cmap='gray',origin='lower')
    plt.axis('off')
    plt.title(f'{im_name}')
    plt.savefig(f'{dir}/{im_name}.png',facecolor='white')
    plt.close()

def Redux(obj,ff,output_folder,x_axis=4096,y_axis=4096,contrast=0.02,use_bias=False,dark_exp_time=100):
    """obj:object,ff=filter,dn=directory name"""
    im_table=pd.read_fwf('imagelist.txt')

    #Listas de archivos de lights y flats
    print('Extracting light and flat file names')
    light_list=im_table[(im_table["Type"]=='Light Frame')&(im_table["Object"]==f"{obj}")&\
                        (im_table["Filter"]==f"{ff}")]["File_Name"].tolist()
    if light_list==[]:
        return print('Missing light images!')
    flat_list=im_table[(im_table["Type"]=='Flat Field')&(im_table["Filter"]==f"{ff}")]["File_Name"].tolist()
    if flat_list==[]:
        return print('Missing flats!')

    #Bias
    if use_bias==True:
        print('Working on bias frames')
        bias_list=im_table[im_table["Type"]=='Bias Frame']["File_Name"].tolist() #Lista de archivos
        if bias_list==[]:
            return print('Missing bias frames!')
        bias=Fits_Array(bias_list) #Cubo
        Master_Bias=np.nanmedian(bias,axis=2) #Creación de master bias
        plots(bias,0.25,output_folder,'Bias')
        s_plot(Master_Bias,0.25,output_folder,'Master_bias')
        print('Bias frames and master bias plotted')
        del bias_list, bias
        print('Done with bias step')

    #Cubos de datos y listas de tiempos de exposición
    lights=Fits_Array(light_list,extract_exptime=True) #Tupla cubo,exptimes
    flats=Fits_Array(flat_list,extract_exptime=True) #Tupla cubo,exptimes
    del flat_list

    #Lista reducida de tiempos de exposición
    redux_explist=list(set(lights[1]))
    redux_explist.extend(list(set(flats[1])))

    #Darks
    print('Working on dark frames')
    if use_bias==True:
        dark_list=im_table[(im_table["Type"]=='Dark Frame')&(im_table["Exp_Time_(s)"]==dark_exp_time)]["File_Name"].tolist()
        if dark_list==[]:
            return print(f'Missing {dark_exp_time}sec dark frames!')
        dark_cube=Fits_Array(dark_list)
        plots(dark_cube,0.25,output_folder,f'Darks_{dark_exp_time}sec')
        print(f'{dark_exp_time}sec dark frames plotted')
        for i in range(len(dark_list)):
            dark_cube[:,:,i]=dark_cube[:,:,i]-Master_Bias
        Dark_Current=np.nanmedian(dark_cube,axis=2)/dark_exp_time
        s_plot(Dark_Current,0.25,output_folder,f'Dark_current')
        print('Dark current ready')
        del dark_list, dark_cube
    else:
        Master_Dark_cube=np.empty([x_axis,y_axis,0])
        for i in redux_explist:
            dark_list=im_table[(im_table["Type"]=='Dark Frame')&(im_table["Exp_Time_(s)"]==i)]["File_Name"].tolist()
            if dark_list==[]:
                return print(f'Missing {i}sec dark frames!')
            dark_cube=Fits_Array(dark_list)
            master_dark=np.nanmedian(dark_cube,axis=2)
            plots(dark_cube,0.25,output_folder,f'Darks_{i}sec')
            s_plot(master_dark,0.25,output_folder,f'Master_dark_{i}sec')
            print(f'{i}sec darks and master dark plotted')
            Master_Dark_cube=np.dstack([Master_Dark_cube,master_dark])
        del dark_list, dark_cube, master_dark
    print('Done with dark stage')

    #Flats
    print('Working on flat frames')
    redu_flat_cube=np.empty([x_axis,y_axis,0])
    if use_bias==True:
        for i in range(len(flats[1])):
            reduced_flat=flats[0][:,:,i]-flats[1][i]*Dark_Current-Master_Bias
            redu_flat_cube=np.dstack([redu_flat_cube,reduced_flat])
        Norm_master_flat=np.nanmedian(redu_flat_cube,axis=2)/np.nanmedian(np.nanmedian(redu_flat_cube,axis=2))    
        plots(flats[0],0.5,output_folder,f"Flats_{ff}")
        print(f'Flats {ff} plotted')
    else:
        for i in range(len(flats[1])):
            reduced_flat=flats[0][:,:,i]-Master_Dark_cube[:,:,np.array(redux_explist)==flats[1][i]].squeeze(axis=2)
            redu_flat_cube=np.dstack([redu_flat_cube,reduced_flat])
        Norm_master_flat=np.nanmedian(redu_flat_cube,axis=2)/np.nanmedian(np.nanmedian(redu_flat_cube,axis=2))
        plots(flats[0],0.5,output_folder,f"Flats_{ff}")
        print(f'Filter {ff} flats plotted')
    s_plot(Norm_master_flat,0.5,output_folder,f'Master_flat_{ff}')
    print(f'Normalized Master Flat filter {ff} plotted')
    del flats, redu_flat_cube, reduced_flat
    print('Done with flat stage')

    #Reduccion
    print('Reduction stage started')
    reducc=np.empty([x_axis,y_axis,0])
    if use_bias==True:
        for i in range(len(lights[1])):
            oper=lights[1][i]*(lights[0][:,:,i]-lights[1][i]*Dark_Current-Master_Bias)\
            /Norm_master_flat
            reducc=np.dstack((reducc,oper))
    else:    
        for i in range(len(lights[1])):
            oper=lights[1][i]*((lights[0][:,:,i]-Master_Dark_cube[:,:,np.array(redux_explist)==lights[1][i]].squeeze(axis=2)))\
            /Norm_master_flat
            reducc=np.dstack((reducc,oper))
    BPM=fits.open('BPM.fit')[0].data
    prod=(np.sum(reducc,axis=2)/np.sum(np.array(lights[1])))*BPM
    plots(lights[0],0.25,output_folder,f"Raw_images_{obj}_{ff}")
    print(f'Raw {obj} images plotted')
    s_plot(prod,contrast,output_folder,f'Reduced_{obj}_{ff}')
    print(f'Reduced {obj} image plotted')
    del lights

    #Fits
    print('Creating FITS file')
    phdu=fits.PrimaryHDU(header=fits.open(light_list[0])[0].header,data=prod)
    del light_list
    phdu.header['PIXSCALE'] = 0.36
    phdu.header['CTYPE1'] = 'RA---TAN'
    phdu.header['CTYPE2'] = 'DEC--TAN'
    phdu.header['CRPIX1'] = int(phdu.header['NAXIS1']/2)
    phdu.header['CRPIX2'] = int(phdu.header['NAXIS2']/2)
    coord=SkyCoord(f'{phdu.header['OBJCTRA']} {phdu.header['OBJCTDEC']}', unit=(u.hourangle, u.deg))
    phdu.header['CRVAL1'] = coord.ra.value
    phdu.header['CRVAL2'] = coord.dec.value
    phdu.header['CD1_1'] = -0.0001
    phdu.header['CD1_2'] = 0
    phdu.header['CD2_1'] = 0
    phdu.header['CD2_2'] = 0.0001

    phdu.writeto(f'{output_folder}/Reduced_{obj}_{ff}.fits',overwrite=True)
    
    return print('Done!')

main()