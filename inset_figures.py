import numpy as np
from astropy import units as u
from astropy import coordinates
import pylab as pl
from astropy.io import fits
from astropy import wcs
import astropy.visualization
from mpl_plot_templates import asinh_norm
import matplotlib
import warnings
from visualization import make_scalebar
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

# this is just to make the script quieter; it's not necessary to run
warnings.filterwarnings('ignore', category=wcs.FITSFixedWarning, append=True)

distance = 415.*u.pc

pl.rcParams['image.interpolation'] = 'nearest'
pl.rcParams['image.origin'] = 'lower'
if int(matplotlib.__version__[0]) >= 2:
    pl.rcParams['figure.dpi'] = 75.
    pl.rcParams['savefig.dpi'] = 300.
    pl.rcParams['axes.labelsize'] = 8
    pl.rcParams['axes.titlesize'] = 8
    pl.rcParams['xtick.labelsize'] = 6
    pl.rcParams['ytick.labelsize'] = 6
    pl.rcParams['axes.linewidth'] = 0.15
    tick_fontsize = 6
    markersize = 3
else:
    pl.rcParams['savefig.dpi'] = 300.
    tick_fontsize = 10
    markersize = 3

# We predefine a set of zoom regions here that get used below
# The keywords should be relatively self-explanatory;
# the keys in the dictionary are just region names
# (this allows you to make different figures for different regions;
# for example, if you wanted a figure focused on the ISF but with
# zoom-ins of different regions around the ISF, then another figure
# with L1641 and zoom ins around it)
zoomregions = {'kelvinhelmholtzwaves': # this is a sub-region name
               {'bottomleft': coordinates.SkyCoord("5:34:12.841",
                                                   "-5:27:56.510",
                                                   unit=(u.h, u.deg),
                                                   frame='fk5'),
                'topright': coordinates.SkyCoord("5:33:45.184",
                                                 "-5:25:46.816",
                                                 unit=(u.h, u.deg),
                                                 frame='fk5'),
                'inregion': 'fullfield', # this tells you what parent region the sub-region is in
                'bbox':[0.55,0.9], # location of the inset axes defined in figure fraction units
                'loc': 1, # location of the inset axes within the bounding box (probably best to leave as 1?)
                'l1':3, # location of the line connecting the box corners (values 1,2,3,4)
                'l2':1,
                'min': 30, # min/max displayed value
                'max': 65,
                'zoom': 5, # zoom factor
               },
              }
# sometimes order matters; dictionaries don't preserve order
zoomregions_order = ['kelvinhelmholtzwaves']

# define legend locations
leglocs = {'fullfield': 'upper right'}
legbboxes = {'fullfield': [1,1]}

# define display limits
vlimits = {'fullfield': {'12CO': (-0.1, 200),}}

# a dictionary of filenames, in case you want to run the zooms on a set of different images
filenames = {'12CO': '12co_pix_2_Tmb.max.fits',
            }

# define scalebar location
scalebarpos = {'fullfield': coordinates.SkyCoord("5:34:33.892", "-6:04:36.568", unit=(u.h, u.deg), frame='fk5')}
scalebarlength = {'fullfield': 1*u.pc}
scalebarlabel = {'fullfield': "1 pc"}

# define field zoom region for the large-scale field
field_corners = {'fullfield':
                 {'bottomleft': coordinates.SkyCoord("5:38:48.102",
                                                     "-7:10:34.015",
                                                     unit=(u.h, u.deg),
                                                     frame='fk5'),
                  'topright': coordinates.SkyCoord("5:33:12.020",
                                                   "-4:49:59.662",
                                                   unit=(u.h, u.deg),
                                                   frame='fk5'),
                 }
                }

# define colormaps & norms
cmaps = {'fullfield': pl.cm.gray_r}
norms = {'fullfield': asinh_norm.AsinhNorm()}

# define colorbar labels (i.e., units)
cb_labels = {'12CO': "$T_B$ [K]"}
             # or something like: "$S_{3 mm}$ [mJy beam$^{-1}$]"}


# Loop over the large-scale regions, one figure per region
for regionname in ('fullfield', ):

    legloc = leglocs[regionname]
    leg_bbox = legbboxes[regionname]

    for line in ("12CO",):

        hdu_line = fits.open(filenames[line])[0]

        wcsaxes = mywcs = wcs.WCS(hdu_line.header).celestial

        fig3 = pl.figure(3)
        fig3.clf()
        ax = fig3.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcsaxes)

        ra = ax.coords['ra']
        ra.set_major_formatter('hh:mm:ss.s')
        dec = ax.coords['dec']
        ra.set_axislabel("RA (J2000)", fontsize=pl.rcParams['axes.labelsize'])
        dec.set_axislabel("Dec (J2000)", fontsize=pl.rcParams['axes.labelsize'], minpad=0.0)
        ra.ticklabels.set_fontsize(tick_fontsize)
        ra.set_ticks(exclude_overlapping=True)
        dec.ticklabels.set_fontsize(tick_fontsize)
        dec.set_ticks(exclude_overlapping=True)

        vmin, vmax = vlimits[regionname][line]

        bottomleft, topright = field_corners[regionname]['bottomleft'], field_corners[regionname]['topright'],

        im = ax.imshow(hdu_line.data.squeeze(),
                       transform=ax.get_transform(mywcs),
                       vmin=vmin, vmax=vmax, cmap=cmaps[regionname],
                       interpolation='nearest',
                       origin='lower', norm=norms[regionname])
        tr_fk5 = ax.get_transform("fk5")
        #(x1,y1),(x2,y2) = (1200,434),(2142,1743)
        # wrong (x1,y1),(x2,y2) = tr_fk5.transform_point([bottomleft.ra.deg, bottomleft.dec.deg]),tr_fk5.transform_point([topright.ra.deg, topright.dec.deg])
        (x1,y1),(x2,y2) = (mywcs.wcs_world2pix([[bottomleft.ra.deg,
                                                 bottomleft.dec.deg]],0)[0],
                           mywcs.wcs_world2pix([[topright.ra.deg,
                                                 topright.dec.deg]],0)[0]
                          )

        make_scalebar(ax, scalebarpos[regionname],
                      length=(scalebarlength[regionname] /
                              distance).to(u.arcsec, u.dimensionless_angles()),
                      color='k',
                      label=scalebarlabel[regionname],
                      text_offset=1.0*u.arcmin,
                     )


        ax.set_aspect(1)
        ax.axis([x1,x2,y1,y2])


        for zoomregion in zoomregions_order:

            ZR = zoomregions[zoomregion]
            if ZR['inregion'] != regionname:
                continue

            parent_ax = zoomregions[ZR['inset_axes']]['axins'] if 'inset_axes' in ZR else ax

            bl, tr = ZR['bottomleft'],ZR['topright'],
            (zx1,zy1),(zx2,zy2) = (mywcs.wcs_world2pix([[bl.ra.deg,
                                                         bl.dec.deg]],0)[0],
                                   mywcs.wcs_world2pix([[tr.ra.deg,
                                                         tr.dec.deg]],0)[0]
                                  )

            axins = zoomed_inset_axes(parent_ax, zoom=ZR['zoom'], loc=ZR['loc'],
                                      bbox_to_anchor=ZR['bbox'],
                                      bbox_transform=fig3.transFigure,
                                      axes_class=astropy.visualization.wcsaxes.core.WCSAxes,
                                      axes_kwargs=dict(wcs=wcsaxes))
            ZR['axins'] = axins
            imz = axins.imshow(hdu_line.data.squeeze(),
                               transform=parent_ax.get_transform(mywcs),
                               vmin=ZR['min'], vmax=ZR['max'],
                               cmap=cmaps[regionname],
                               interpolation='nearest',
                               origin='lower', norm=norms[regionname])


            ax.axis([x1,x2,y1,y2])
            axins.axis([zx1,zx2,zy1,zy2])

            axins.set_xticklabels([])
            axins.set_yticklabels([])

            lon = axins.coords['ra']
            lat = axins.coords['dec']
            lon.set_ticklabel_visible(False)
            lat.set_ticklabel_visible(False)

            # draw a bbox of the region of the inset axes in the parent axes and
            # connecting lines between the bbox and the inset axes area
            mark_inset(parent_axes=parent_ax, inset_axes=axins,
                       loc1=ZR['l1'], loc2=ZR['l2'], fc="none", ec="0.5")


            fig3.canvas.draw()
            assert np.abs(ax.bbox._bbox.x1 - 0.95) > 1e-4

            cax = fig3.add_axes([ax.bbox._bbox.x1+0.01, ax.bbox._bbox.y0, 0.02,
                                 ax.bbox._bbox.y1-ax.bbox._bbox.y0])
            cb = fig3.colorbar(mappable=im, cax=cax)
            cb.set_label(cb_labels[line])

        fig3.savefig("inset_figure_{line}_{regionname}.png"
                     .format(line=line,regionname=regionname),
                     bbox_inches='tight')
        fig3.savefig("inset_figure_{line}_{regionname}.pdf"
                     .format(line=line,regionname=regionname),
                     bbox_inches='tight')

        # if you want to overplot symbols etc.
        #leg = ax.legend(handles=coredots, labels=[cd.get_label() for cd in
        #                                          coredots],
        #                bbox_to_anchor=leg_bbox,
        #                loc=legloc,
        #                fontsize=8)

        #fig3.savefig("inset_figure_{line}_{regionname}_withlegend.png"
        #             .format(line=line,regionname=regionname),
        #             bbox_inches='tight')
        #fig3.savefig("inset_figure_{line}_{regionname}_withlegend.pdf"
        #             .format(line=line,regionname=regionname),
        #             bbox_inches='tight')
