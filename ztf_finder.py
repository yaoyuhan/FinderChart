""" Finder chart script, using ztfquery """

import numpy as np
from astropy.io import fits as pyfits
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import sys
import astropy.wcs
from astropy.io import fits
from scipy.ndimage.filters import gaussian_filter
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord
from ztfquery import query,marshal
from ztfquery.io import download_single_url
from ztfquery.io import get_cookie


def deg2hour(ra, dec, sep=":"):
    '''
    Transforms the coordinates in degrees into HH:MM:SS DD:MM:SS with the requested separator.
    '''
    if ( type(ra) is str and type(dec) is str ):
        return ra, dec
    c = SkyCoord(ra, dec, frame='icrs', unit='deg')
    ra = c.ra.to_string(unit=u.hourangle, sep=sep, precision=2, pad=True)
    dec = c.dec.to_string(sep=sep, precision=2, alwayssign=True, pad=True)
    return str(ra), str(dec)


def hour2deg(ra, dec):
    '''
    Transforms string HH:MM:SS DD:MM:SS coordinates into degrees (floats).
    '''
    try:
        ra = float(ra)
        dec = float(dec)
        
    except:
        c = SkyCoord(ra, dec, frame='icrs', unit=(u.hourangle, u.deg))
        
        ra = c.ra.deg
        dec = c.dec.deg
    
    return ra, dec


def get_offset(ra1, dec1, ra2, dec2):
    '''
    Code from Nadia
    Computes the offset in arcsec between two coordinates.
    The offset is from (ra1, dec1) which is the offset star
    to ra2, dec2 which is the fainter target
    '''
    bright_star = SkyCoord(ra1, dec1, frame='icrs', unit=(u.deg, u.deg))
    target = SkyCoord(ra2, dec2, frame='icrs', unit=(u.deg, u.deg))
    dra, ddec = bright_star.spherical_offsets_to(target)
    return dra.to(u.arcsec).value, ddec.to(u.arcsec).value 


def get_pos(name):
    """
    Get position of target, if you need to
    This is very slow, so it's probably better
    to supply an RA and DEC
    """
    m = marshal.MarshalAccess()
    m.load_target_sources()
    coords = m.get_target_coordinates(name)
    ra = coords.ra.values[0]
    dec = coords.dec.values[0]
    return ra, dec


def get_lc(name):
    """ Get light curve of target """
    m = marshal.MarshalAccess()
    m.download_lightcurve(name) 
    lcpath = './Data/marshal/lightcurves/'+name+'/marshal_plot_lc_lightcurve_'+name+'.csv'
    lc = pd.read_csv(lcpath)
    # lc = marshal.get_local_lightcurves(name)
    # lc_dict = [lc[key] for key in lc.keys()][0]
    return lc#_dict


def get_refstars(xpos, ypos, cat):
    """
    Select reference stars.
    
    Parameters
    ----------
    xpos: x position of target
    ypos: y position of targethttps://github.com/yaoyuhan/FinderChart.git
    cat: ZTF source catalog
    """
    sep_pix = np.sqrt(
            (sourceinfo['xpos']-allsourcesinfo['xpos'])**2 + \
            (sourceinfo['ypos']-allsourcesinfo['ypos'])**2)
    
    # should be separated by at least 10 pixels
    crit_a = np.logical_and(sep_pix > 10, allsourcesinfo['flags']==0)
    crit_b = np.logical_and(
            allsourcesinfo['chi'] < 2, allsourcesinfo['snr'] > 10)
    crit_c = allsourcesinfo['sharp'] < 0.3
    crit_ab = np.logical_and(crit_a, crit_b)
    crit = np.logical_and(crit_ab, crit_c)
    
    # should be bright
    mag_crit = np.logical_and(
            allsourcesinfo['mag'] >= 15, allsourcesinfo['mag'] <= 19)
    choose_ind = np.where(np.logical_and(crit, mag_crit))
    
    # choose the closest three stars
    nref = 3
    order = np.argsort(sep_pix[choose_ind])
    
    # hack for now
    filt = 'R'
    colors = ['red', 'orange', 'purple']
    shapes = ['box', 'hexagon', 'circle']
    
    refstars = []
    for i in range(0,nref):
        refstars.append({
                          'name':  'S%s' %i,
                          'color':  colors[i],
                          'shape':  shapes[i],
                          'dist':  sep_pix[choose_ind][order][i] * 1.0, # pixel scale
                          'x_sub':  allsourcesinfo['xpos'][choose_ind][order][i],
                          'y_sub':  allsourcesinfo['ypos'][choose_ind][order][i],
                          'ra':  allsourcesinfo['ra'][choose_ind][order][i],
                          'dec':  allsourcesinfo['dec'][choose_ind][order][i],
                          'mag': allsourcesinfo['mag'][choose_ind][order][i],
                          'mag_err': allsourcesinfo['emag'][choose_ind][order][i],
                          'filter': filt
                        })
    return refstars


def choose_ref(zquery, ra, dec):
    """ Choose a reference image to use, and download the file
    including the associated PSF cat

    Parameters
    ----------
    ra: position of source in decimal degrees
    dec: position of source in decimal degrees

    Returns
    -------
    the filename of what you just downloaded (IPAC link)
    the location of the file (local link)
    """
    zquery.load_metadata(kind="ref",radec=[ra, dec], size=0.0001)
    out = zquery.metatable
    # If no files are returned,
    if len(out) == 0:
        print("Error: couldn't find any reference at this position.")
    else:
        # choose the index of the file with the deepest maglimit
        ind = out['maglimit'].idxmax()
        ind = 1
        urls, dl_loc = zquery.download_data(nodl=True)
        imfile = dl_loc[ind]
        # Temp bug fix: check to make sure the URL is correct
        imfile_array = imfile.split("/")
        imfile_array[3] = imfile_array[4][5:8]
        imfile = '/'.join(imfile_array)
        # default is refimg
        download_single_url(urls[ind], dl_loc[ind], cookies=None) 
        # download the associated PSFcat
        urls, dl_loc = zquery.download_data(nodl=True, suffix='refpsfcat.fits')
        catfile = dl_loc[ind]
        download_single_url(
                urls[ind], dl_loc[ind], cookies=None)
        return imfile, catfile


def choose_sci(zquery, out, name, ra, dec):
    """ 
    Choose a science image to use, and download the file 
    """
    # for science object:
    if len(name)==12: 
        lc = get_lc(name)
        # Count the number of detections where limmag > 19.5
        # If 0, allow limmag > 19
        limmag = lc.limmag.values
        limmag_val = 19.5 # must be deeper than this value
        choose = limmag > limmag_val
        while sum(choose) == 0:
            limmag_val += 0.5
            choose = limmag > limmag_val

        if sum(choose) > 1:
            # Of all these images, choose the one where the transient is brightest
            ind = np.argmin(lc.mag.values[choose])
            jd_choose = lc['jdobs'][choose].values[ind] 
            # mag_choose = lc['mag'][choose].values[ind]
            # filt_choose = lc['filter'][choose].values[ind]
        elif sum(choose) == 1:
            # If there is only one choices...
            jd_choose = lc['jdobs'][choose].values[0]
            # mag_choose = lc['mag'][choose].values[0]
            # filt_choose = lc['filter'][choose].values[0]
            
        # Download the corresponding science image
        ind = np.argsort(np.abs(out.obsjd-jd_choose))[0]
        urls, dl_loc = zquery.download_data(nodl=True)
        imfile = dl_loc[ind]
        download_single_url(urls[ind], dl_loc[ind], cookies=None)
        urls, dl_loc = zquery.download_data(nodl=True, suffix='psfcat.fits')
        catfile = dl_loc[ind]
        download_single_url(urls[ind], dl_loc[ind], cookies=None)    
    # for host galaxy:
    elif len(name)==11:
        urls, dl_loc = zquery.download_data(nodl=True)
        urls = np.array(urls)
        dl_loc = np.array(dl_loc)
        assert len(urls)==len(out)
        
        # restrction on filter -- use r band
        ix = out['filtercode'].values=="zr"
        out = out[ix]
        urls = urls[ix]
        dl_loc = dl_loc[ix]
        
        # restrction on seeing
        seeing_cut = 2
        ix = out['seeing'] < seeing_cut
        if np.sum(ix)<10:
            seeing_cut = 2.5
            ix = out['seeing'].values < seeing_cut
        out = out[ix]
        urls = urls[ix]
        dl_loc = dl_loc[ix]
        
        # consider: host around the center of the image
        cinput = SkyCoord(ra, dec, frame='icrs', unit='deg')
        cimages = SkyCoord(out['ra'].values, out['dec'].values, frame='icrs', unit='deg')
        seps = cinput.separation(cimages).deg
        ix = seps < 0.45
        out = out[ix]
        urls = urls[ix]
        dl_loc = dl_loc[ix]
        
        # Download the corresponding science image
        ind = len(dl_loc)-1
        imfile = dl_loc[ind]
        imfile_url = urls[ind]
        download_single_url(imfile_url, imfile, cookies=None)
        catfile = imfile[:-11]+"psfcat.fits"
        catfile_url = imfile_url[:-11]+"psfcat.fits"
        download_single_url(catfile_url, catfile, cookies=None)    
    
    return imfile, catfile


def get_finder(ra, dec, name, rad=0.01, 
               target_mag = np.nan,
               starlist=None, print_starlist=True, 
               telescope="P200", minmag=15, maxmag=18.5):
    """ 
    Aim: Generate finder chart 
    
    Parameters:
    ra: float or string. 
        eg: ra="18h24m25.36s", 210.437583
        
    dec: float or string, type must be consistent with ra. 
        eg: dec="+44d07m50.0s", 46.215583
        
    name: ZTFname or host name. 
        eg: name="18aaslhxt_h", "ZTF18abclfee"
        
    rad: search radius in the unit of degree
    
    target_mag: magnitude of the target at time of observation, 
        default shoudl be r band
    
    starlist: name of the starlist
        eg: starlist="/Users/yuhanyao/Desktop/observation/20190428_LRIS/starlist"
        
        
    References:
    Keck starlist format: https://www2.keck.hawaii.edu/observing/starlist.html
        First field is the target name in columns 1-16 (tabs and spaces allowed). Maximum length is 15 characters.
        Next three space-separated tokens (beginning in column 17) are RA in the form HH MM SS.SS (including an arbitrary number of decimal places).
        Next three space-separated tokens are Dec in the form (-)DD MM SS.S (again, to an arbitrary number of decimal places). Note: if the Dec is between 0 and -1, the DD field MUST be -00).
        Next token is equinox in the form YYYY.YY (no parentheses; arbitrary number of decimal places), with <=1950 meaning FK4 and 2000 meaning FK5 and APP meaning apparent.
        
    """
    print ('Using search radius of %.1f arcsec.'%(rad*3600))
    
    name = str(name)
    assert len(name) in [12, 11]
    
    assert type(ra)==type(dec)
    if type(ra)==str:
        c1 = SkyCoord(ra, dec, frame='icrs')
        ra = c1.ra.degree
        dec = c1.dec.degree
    ra = float(ra)
    dec = float(dec)
    
    # prepare to print starlist
    if telescope == "Keck":
        name_length = 16 # maximum length of name is 15
        commentchar = "#"
        separator = ""
    elif telescope == "P200":
        name_length = 20
        commentchar = "!"
        separator = "!"
    
    #Write to the starlist if the name of the starlist was provided.
    rah, dech = deg2hour(ra, dec, sep=" ")
    if (not starlist is None) and (telescope =="Keck"):
        with open(starlist, "a") as f:
            f.write( "{:s}{:s} {:s} 2000.0 {:s} {:f} \n".format(name.ljust(name_length), rah, dech, commentchar, target_mag) ) 
            f.close()

    # Get metadata of all images at this location
    print("Querying for metadata...")
    zquery = query.ZTFQuery()
    zquery.load_metadata(radec=[ra,dec], size=rad)
    out = zquery.metatable

    # Do you need to use a reference image?
    need_ref = len(out) == 0
    if need_ref:
        print("Using a reference image")
        imfile, catfile = choose_ref(zquery, ra, dec)
    else:
        print("Using a science image")
        imfile, catfile = choose_sci(zquery, out, name, ra, dec)

    # get the cutout
    inputf = pyfits.open(imfile)
    im = inputf[0].data
    inputf.close()
    head = fits.getheader(imfile)

    # Get the x and y position of the target,
    # as per the IPAC catalog
    wcs = astropy.wcs.WCS(head)
    world =  np.array([[ra, dec]], np.float_)
    target_pix = wcs.wcs_world2pix(world, 0)[0]
    xpos = target_pix[0]
    ypos = target_pix[1]

    # adjust counts
    im[np.isnan(im)] = 0
    im[im > 30000] = 30000
    
    # extract 600x600 region around the position of the target
    width = 600
    height = 600
    xmax = xpos + width/2
    xmin = xpos - width/2
    ymax = ypos + height/2
    ymin = ypos - height/2

    plt.figure(figsize=(8,6))
    plt.set_cmap('gray_r')
    smoothedimage = gaussian_filter(im, 1.3)
    # pad the image
    im_padded = np.pad(smoothedimage, 300, mode='constant', constant_values=0)

    # If it's a reference image, you have to flip it up/down and left/right
    if need_ref:
        croppedimage = np.fliplr(np.flipud(
            im_padded[int(ymin)+300:int(ymax)+300,
                int(xmin)+300:int(xmax)+300]))

    # If it's a science image, you just flip it up/down
    else:
        croppedimage = np.flipud(
            im_padded[int(ymin)+300:int(ymax)+300,
                int(xmin)+300:int(xmax)+300])

    plt.imshow(
            croppedimage, origin='lower', # convention for IPAC images
            vmin=np.percentile(im.flatten(), 10),
            vmax=np.percentile(im.flatten(), 99.0))

    # Mark target: should just be the center of the image, now
    # horizontal line
    plt.plot([300+5,300+20],[300,300], 'g-', lw=2)
    # vertical line
    plt.plot([300,300],[300+5,300+20], 'g-', lw=2)

    # and the offset of the original coordinate system with the new coordinates
    offset_x = xpos-300
    offset_y = ypos-300

    # Choose offset stars
    cat = pyfits.open(catfile)[1].data
    zp = pyfits.open(catfile)[0].header['MAGZP']
    sep_pix = np.sqrt(
            (xpos-cat['xpos'])**2 + \
            (ypos-cat['ypos'])**2)

    # should be separated by at least 10 pixels
    crit_a = np.logical_and(sep_pix > 10, cat['flags']==0)
    crit_b = np.logical_and(
            cat['chi'] < 2, cat['snr'] > 10)
    crit_c = cat['sharp'] < 0.3
    crit_ab = np.logical_and(crit_a, crit_b)
    crit = np.logical_and(crit_ab, crit_c)

    # should be bright
    mag_crit = np.logical_and(cat['mag']+zp >= minmag, cat['mag']+zp <= maxmag)
    choose_ind = np.where(np.logical_and(crit, mag_crit))

    # mark the closest three stars
    nref = 3
    order = np.argsort(sep_pix[choose_ind])
    cols = ['orange', 'purple', 'red']

    for ii in np.arange(nref):
        ref_xpos_original = cat['xpos'][choose_ind][order][ii] - offset_x
        ref_ypos_original = cat['ypos'][choose_ind][order][ii] - offset_y

        # transform to flipped plot
        if need_ref:
            ref_xpos = 600-ref_xpos_original
            ref_ypos = 600-ref_ypos_original
        else:
            ref_xpos = ref_xpos_original
            ref_ypos = 600-ref_ypos_original

        plt.plot(
                [ref_xpos+5,ref_xpos+20],[ref_ypos,ref_ypos], 
                c=cols[ii], ls='-', lw=2)
        plt.plot(
                [ref_xpos,ref_xpos],[ref_ypos+5,ref_ypos+20], 
                c = cols[ii], ls='-', lw=2) 
        refra = cat['ra'][choose_ind][order][ii]
        refdec = cat['dec'][choose_ind][order][ii]
        if telescope == 'Keck':
            refrah, refdech = deg2hour(refra, refdec,sep=" ")
        elif telescope == 'P200':
            refrah, refdech = deg2hour(refra, refdec,sep=":")
        else:
            print("I don't recognize this telescope")
        refmag = cat['mag'][choose_ind][order][ii]+zp
        dra, ddec = get_offset(refra, refdec, ra, dec)

        offsetnum = 0.2
        plt.text(
                1.02, 0.60-offsetnum*ii, 
                'Ref S%s, mag %s' %((ii+1), np.round(refmag,1)), 
                transform=plt.axes().transAxes, fontweight='bold', color=cols[ii])
        plt.text(
                1.02, 0.55-offsetnum*ii, 
                '%s %s' %(refrah,refdech),
                color=cols[ii], transform=plt.axes().transAxes)
        plt.text(
                1.02, 0.50-offsetnum*ii, 
                str(np.round(ddec,2)) + "'' N, " + str(np.round(dra,2)) + "'' E", 
                color=cols[ii], transform=plt.axes().transAxes)

        # Print starlist for telescope
        if telescope == 'Keck':
            # Target name is columns 1-16
            # RA must begin in 17, separated by spaces
            print ("{:s}{:s} {:s} 2000.0 {:s} raoffset={:.2f} decoffset={:.2f} {:s} r={:.1f} ".format((name+"_S%s" %(ii+1)).ljust(name_length), refrah, refdech, separator, dra, ddec, commentchar, refmag))
        elif telescope == 'P200':
            print ("{:s} {:s} {:s}  2000.0 {:s} raoffset={:.2f} decoffset={:.2f} r={:.1f} {:s} ".format((name+"_S%s" %(ii+1)).ljust(name_length), refrah, refdech, separator, dra, ddec, refmag, commentchar))

        
        # and save to file if starlist name is provided
        if (not starlist is None) and (telescope =="Keck"):
            with open(starlist, "a") as f:
                f.write ("{:s}{:s} {:s} 2000.0 {:s} raoffset={:.2f} decoffset={:.2f} {:s} r={:.1f} \n".format((name+"_S%s" %(ii+1)).ljust(name_length), refrah, refdech, separator, dra, ddec, commentchar, refmag))
                f.close()
        elif (not starlist is None) and (telescope =="P200"):
            with open(starlist, "a") as f:
                f.write("{:s} {:s} {:s}  2000.0 {:s} raoffset={:.2f} decoffset={:.2f} r={:.1f} {:s} \n".format((name+"_S%s" %(ii+1)).ljust(name_length), refrah, refdech, separator, dra, ddec, refmag, commentchar))
                f.close()

    # Plot compass
    plt.plot(
            [width-10,height-40], [10,10], 'k-', lw=2)
    plt.plot(
            [width-10,height-10], [10,40], 'k-', lw=2)
    plt.annotate(
        "N", xy=(width-20, 40), xycoords='data',
        xytext=(-4,5), textcoords='offset points')
    plt.annotate(
        "E", xy=(height-40, 20), xycoords='data',
        xytext=(-12,-5), textcoords='offset points')

    # Get rid of axis labels
    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)

    # Set size of window (leaving space to right for ref star coords)
    plt.subplots_adjust(right=0.65,left=0.05, top=0.99, bottom=0.05)

    # List name, coords, mag of the target
    plt.text(1.02, 0.85, name, transform=plt.axes().transAxes, fontweight='bold')
    # Can't print mag, because we don't know how bright the target is
    #plt.text(1.02, 0.80, "%s"%mag, transform=plt.axes().transAxes, fontweight='bold')
    plt.text(1.02, 0.80, "%.5f %.5f"%(ra, dec),transform=plt.axes().transAxes)
    rah, dech = deg2hour(ra, dec)
    plt.text(1.02, 0.75,rah+"  "+dech, transform=plt.axes().transAxes)
    plt.savefig("finder_chart_%s.png" %name)
    plt.close()

 
    #get_finder(ra, dec, name, rad, telescope=telescope, debug=False, minmag=7, maxmag=18)
