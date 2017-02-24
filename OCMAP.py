# from astropy.table import QTable, hstack
import numpy as np
import time
import datetime
import glob
import itertools
import matplotlib.pyplot as plt
import os

class OpenClusters:
    """
    Container for open clusters.
    """

    __author__ = ['Locke Patton', 'Ellis Avelone','Katie Crotts']

    def __init__(self, cluster, cluster_title, filters_images,
                 path_in_cluster, path_in_standards, path_out, x_center, y_center,
                 t=None, verbose=True, verbose_absolute=True):
        """
        Parameters:
       	cluster: 	       string corresponding to the name of a cluster inside files
        cluster_title:     string corresponding to the desired plot titles.
        path_in_cluster:   string corresponding to the location of cluster .als files
        path_in_standards: string corresponding to the location of standard .mag files
        filters_images:    2D array-like. contains strings corresponding to the filters used and their corresponding image base names.
        path_out: 	       string corresponding to the desired saving location for output files and plots.
        verbose: 	       Boolean that lets OCMAP know whether or not to print extra information.
        verbose_absolute:  Boolean that lets OCMAP know whether or not to print important information.
        """

        #details about cluster under consideration
        self.cluster = cluster
        self.cluster_title = cluster_title

        #paths for reading in open cluster magnitude files, standard star magnitude files and output files
        self.path_in_cluster = path_in_cluster
        self.path_in_standards = path_in_standards
        self.path_out = path_out

        #defining global verbose and verbose_absolute variables for printing
        self.verbose = verbose
        self.verbose_absolute = verbose_absolute
        self.x_center = x_center
        self.y_center = y_center

        #checking to see if user specified an output time, t, for use in naming output files
        if t == None:
            self.t = datetime.datetime.now()
        elif t != None:
            self.t = t

        #printing output time in files
        if verbose_absolute:
            print 'output time seen in file names:', self.t

        #dictionary to contain magnitudes and images names, etc
        self.Centers = {}

        #reading in magnitudes from filters and image array
        self.filters, self.image_names = filters_images
        self.nstars = []

        for filter_, image_name_ in zip(self.filters,self.image_names):
            # building readin command
            iraf_als_file = glob.glob(self.path_in_cluster + image_name_)

            # opening iraf als photometry file
            with open(iraf_als_file[0]) as f_in:
                # intertools.islice slices file to only obtain mag line
                iraf_als = np.genfromtxt(itertools.islice(f_in, 0, None, 2),
                                         dtype=[('ID', '<f8'), ('XCENTER', '<f8'), ('YCENTER', '<f8'), ('MAG', '<f8'),
                                                ('MERR', '<f8'), ('MSKY', '<f8'), ('NITER', '<f8')])

            n_stars = len(iraf_als['ID'])
            self.nstars.append(n_stars)

            if self.verbose:
                print 'iraf_als',iraf_als

            self.Centers[filter_] = {}
            self.Centers[filter_]['ID'] = np.array(range(n_stars))
            self.Centers[filter_]['MAG'] = iraf_als['MAG']
            self.Centers[filter_]['MERR'] = iraf_als['MERR']
            self.Centers[filter_]['XCENTER'] = iraf_als['XCENTER']
            self.Centers[filter_]['YCENTER'] = iraf_als['YCENTER']
            self.Centers[filter_]['NSTARS'] = n_stars



            # printing filters, image name, iraf_als_file if verbose
            if self.verbose_absolute:
                print n_stars, 'stars in', filter_, 'from', image_name_
            if self.verbose:
                print iraf_als_file

path = r'/Cygwin/home/Katie/clusters/'
flname = glob.glob(path + '*turner11xy.txt')

    def plotXYCenter(self):
        plt.plot(self.x_center, self.y_center, linestyle='None', marker='o', markersize=10, alpha=.4)
        plt.xlabel('X PIX')
        plt.ylabel('Y PIX')
        plt.title(self.cluster_title + 'X-Center vs. Y-Center', fontsize='16')
        plt.legend(["b filter", "v filter", "y filter"], loc='upper left', fancybox=True, numpoints=1)

    @classmethod
    def from_txt(cls, flname):
        data = np.loadtxt(flname)
        print data
        x_center = data[:, 0]
        y_center = data[:, 1]
        return cls(x_center=x_center, y_center=y_center)

xycenters = [OpenClusters.from_txt(path) for path in flname]

for openclusters in xycenters:
    print(openclusters.x_center.mean())

    for openclusters in xycenters:
        openclusters.plot()
    fig = plt.gcf()
    fig.set_size_inches(8, 7)
    plt.xlim(0, 200)
    plt.ylim(0, 200)

        pass

    def PositionMatch(self, tol, n_iterations=None, shifts=None):
        """
        Matches the positions of data points across multiple images and filters.

        Parameters:
        self:   an OpenClusters object
        tol:    tolerance in matching files
        n_iterations:   user can run matching for any number of iterations less than total length of base filter stars
        shifts: an array containing horizontal and vertical shift of filters
            self: an OpenClusters object
            tol: tolerance in matching stars on cmb (should be in units of magnitudes inputted)
            shifts: an array containing horizontal and vertical position shifts between separate images (filters)
            image_names: names of images (in IRAF als format).
        """
        if self.verbose:
            print "\nRunning PositionMatch"
        #Determining the shifts for each filter in each filter in x and y direction

        #default shifts are 0s
        if shifts == None:
            self.shifts = []
            for filter_ in self.filters:
                self.shifts.append([0., 0.])
        else:
            self.shifts = shifts

        self.StarMatch = {}
        self.StarMatch_extra = {}

        for it_,(filter_, image_name_) in enumerate(zip(self.filters, self.image_names)):

            self.Centers[filter_]['XCENTER_SHIFTED'] = self.Centers[filter_]['XCENTER'] + self.shifts[it_][0]
            self.Centers[filter_]['YCENTER_SHIFTED'] = self.Centers[filter_]['YCENTER'] + self.shifts[it_][1]

            #printing filters, image name, iraf_als_file if verbose
            if self.verbose:
                print "filter, image:", filter_, image_name_
                print 'x, y shift:', self.shifts[it_]

            #building matching framework
            self.StarMatch[filter_]={}
            self.StarMatch[filter_]['ID'] = []
            self.StarMatch[filter_]['MAG'] = []
            self.StarMatch[filter_]['MERR'] = []
            self.StarMatch[filter_]['XCENTER'] = []
            self.StarMatch[filter_]['YCENTER'] = []
            self.StarMatch[filter_]['NSTARS'] = self.Centers[filter_]['NSTARS']
            self.StarMatch_extra[filter_] = {}

        #defining base filter as one with most stars
        magBase_index = np.argmax(self.nstars)
        magBase = self.filters[magBase_index]

        #defining filters 1 and 2 as first inputed into filters list (apart from the filter with max # of stars)
        mag1,mag2 = np.delete(self.filters, magBase_index)[:2]

        self.StarMatch_extra[mag1+'_'+magBase+'_radius'] = []
        self.StarMatch_extra[mag2+'_'+magBase+'_radius'] = []

        if self.verbose:
            print 'base filter:',magBase
            print 'matched filters:',mag1, mag2

        id_base = range(self.StarMatch[magBase]['NSTARS'])
        # id_mag1 = range(self.StarMatch[mag1]['NSTARS']) # v_id
        # id_mag2 = range(self.StarMatch[mag2]['NSTARS']) # b_id

        x_cen_magBase = self.Centers[magBase]['XCENTER_SHIFTED'] # y_xcen
        x_cen_mag1 = self.Centers[mag1]['XCENTER_SHIFTED'] # v_xcen
        x_cen_mag2 = self.Centers[mag2]['XCENTER_SHIFTED'] # b_xcen

        y_cen_magBase = self.Centers[magBase]['YCENTER_SHIFTED'] # y_ycen
        y_cen_mag1 = self.Centers[mag1]['YCENTER_SHIFTED'] # v_ycen
        y_cen_mag2 = self.Centers[mag2]['YCENTER_SHIFTED'] # b_ycen

        mag_magBase = self.Centers[magBase]['MAG']
        mag_mag1 = self.Centers[mag1]['MAG']
        mag_mag2 = self.Centers[mag2]['MAG']

        merr_magBase = self.Centers[magBase]['MERR']
        merr_mag1 = self.Centers[mag1]['MERR']
        merr_mag2 = self.Centers[mag2]['MERR']

        radius2 = lambda x,y : (x**2+y**2)

        #TODO Go through **2 mistake and re-run all open cluster membership without the mistake

        ratio = tol

        if n_iterations == None:
            n_iterations = len(id_base)

        for base_i in id_base[:n_iterations]:
            x_cen_magBase_i = x_cen_magBase[base_i]
            y_cen_magBase_i = y_cen_magBase[base_i]

            mag1_magBase_radius2 = radius2((x_cen_mag1 - x_cen_magBase_i),(y_cen_mag1 - y_cen_magBase_i))
            mag1_i = np.argmin(mag1_magBase_radius2)

            if mag1_magBase_radius2.min() <= ratio*radius2(self.Centers[mag1]['MERR'][mag1_i],
                                                           self.Centers[magBase]['MERR'][base_i]):

                mag2_magBase_radius2 = radius2((x_cen_mag2 - x_cen_magBase_i),(y_cen_mag2 - y_cen_magBase_i))
                mag2_i = np.argmin(mag2_magBase_radius2)

                if mag2_magBase_radius2.min() <= ratio * radius2(self.Centers[mag2]['MERR'][mag2_i],
                                                                 self.Centers[magBase]['MERR'][base_i]):

                    self.StarMatch[magBase]['ID'].append(base_i)
                    self.StarMatch[magBase]['MAG'].append(mag_magBase[base_i])
                    self.StarMatch[magBase]['MERR'].append(merr_magBase[base_i])
                    self.StarMatch[magBase]['XCENTER'].append(x_cen_magBase_i)
                    self.StarMatch[magBase]['YCENTER'].append(y_cen_magBase_i)

                    self.StarMatch[mag1]['ID'].append(mag1_i)
                    self.StarMatch[mag1]['MAG'].append(mag_mag1[mag1_i])
                    self.StarMatch[mag1]['MERR'].append(merr_mag1[mag1_i])
                    self.StarMatch[mag1]['XCENTER'].append(x_cen_mag1[mag1_i])
                    self.StarMatch[mag1]['YCENTER'].append(y_cen_mag1[mag1_i])

                    self.StarMatch[mag2]['ID'].append(mag2_i)
                    self.StarMatch[mag2]['MAG'].append(mag_mag2[mag2_i])
                    self.StarMatch[mag2]['MERR'].append(merr_mag2[mag2_i])
                    self.StarMatch[mag2]['XCENTER'].append(x_cen_mag2[mag2_i])
                    self.StarMatch[mag2]['YCENTER'].append(y_cen_mag2[mag2_i])

                    self.StarMatch_extra[mag1 + '_' + magBase + '_radius'].append(mag1_magBase_radius2.min())
                    self.StarMatch_extra[mag2 + '_' + magBase + '_radius'].append(mag2_magBase_radius2.min())

        if self.verbose_absolute:
            print '# Stars across filters',len(self.StarMatch[magBase]['ID']), '/', len(self.Centers[magBase]['ID'])

        """
        #Matching stars across filters

        StarMatch = {}
        StarMatch_extra = {}

        StarMatch['v_id']=[]
        StarMatch['v_MAG'] = []
        StarMatch['v_MERR'] = []
        StarMatch['v_XCENTER'] = []
        StarMatch['v_YCENTER'] = []
        StarMatch['b_id']=[]
        StarMatch['b_MAG'] = []
        StarMatch['b_MERR'] = []
        StarMatch['b_XCENTER'] = []
        StarMatch['b_YCENTER'] = []
        StarMatch['y_id']=[]
        StarMatch['y_MAG'] = []
        StarMatch['y_MERR'] = []
        StarMatch['y_XCENTER'] = []
        StarMatch['y_YCENTER'] = []

        StarMatch_extra['by_rad'] = []
        StarMatch_extra['vy_rad'] = []

        # by_err = (als_b['MERR']**2 + als_y['MERR']**2)**.5
        # vy_err = (als_v['MERR']**2 + als_y['MERR']**2)**.5

        ratio = 90
        for y in range(len(y_xcen)):
        # for y in range(10):
            y_xcen_i = y_xcen[y]
            y_ycen_i = y_ycen[y]
            vy_rad = (v_xcen - y_xcen_i)**2 + (v_ycen - y_ycen_i)**2
            v = np.where(vy_rad == vy_rad.min())[0][0]
            if vy_rad.min() <= ratio*(als_v['MERR'][v]**2 + als_y['MERR'][y]**2)**.5:
                by_rad = (b_xcen - y_xcen_i)**2 + (b_ycen - y_ycen_i)**2
                b = np.where(by_rad == by_rad.min())[0][0]
                if by_rad.min() <= ratio*(als_b['MERR'][b]**2 + als_y['MERR'][y]**2)**.5:
                    if verbose == 'no':
                        print np.where(by_rad == by_rad.min())
                        print np.where(by_rad == by_rad.min())[0][0]
                        print by_rad.min()
                        print by_rad[b]
                        print np.where(vy_rad == vy_rad.min())
                        print np.where(vy_rad == vy_rad.min())[0][0]
                        print vy_rad.min()
                        print vy_rad[v]

                    StarMatch['y_id'].append(y)
                    StarMatch['y_MAG'].append(als_y['MAG'][y])
                    StarMatch['y_MERR'].append(als_y['MERR'][y])
                    StarMatch['y_XCENTER'].append(y_xcen_i)
                    StarMatch['y_YCENTER'].append(y_ycen_i)

                    StarMatch['b_id'].append(b)
                    StarMatch['b_MAG'].append(als_b['MAG'][b])
                    StarMatch['b_MERR'].append(als_b['MERR'][b])
                    StarMatch['b_XCENTER'].append(b_xcen[b])
                    StarMatch['b_YCENTER'].append(b_ycen[b])

                    StarMatch['v_id'].append(v)
                    StarMatch['v_MAG'].append(als_v['MAG'][v])
                    StarMatch['v_MERR'].append(als_v['MERR'][v])
                    StarMatch['v_XCENTER'].append(v_xcen[v])
                    StarMatch['v_YCENTER'].append(v_ycen[v])

                    StarMatch_extra['by_rad'].append(by_rad.min())
                    StarMatch_extra['vy_rad'].append(vy_rad.min())

        if verbose_absolute == 'yes':
            print '# Stars across filters',len(StarMatch['y_id']), '/', len(als_v)
        """

    def Standardize(self):
        """
        TODO
        Parameters:
        self: an OpenClusters object
        """
        # TODO : make these standard names files - y_standards_names.txt

        #non-object oriented version from jupyter notebook
        # !cd $path_in_standards"v"; ls * fits.mag * > v_standards_names.txt
        # !cd $path_in_standards"b"; ls * fits.mag * > b_standards_names.txt
        # !cd $path_in_standards"y"; ls * fits.mag * > y_standards_names.txt
        #
        # if verbose == 'yes':
        #     !cd $path_in_standards"v"; cat v_standards_names.txt
        #     !cd $path_in_standards"b"; cat b_standards_names.txt
        #     !cd $path_in_standards "y"; cat y_standards_names.txt

        #example working with directory codes - maybe helpful?
        # directory = self.path + '/plots/'
        # if not os.path.exists(directory):
        #     os.makedirs(directory)
        #
        # filenameend = '.jpg'
        # if x != None:
        #     filenameend = '_x' + str(x[0]) + ':' + str(x[1]) + '.jpg'
        #
        # fig.savefig(directory + self.image + '_apertureplot' + filenameend, bbox_inches='tight')

        self.Standards = {}
        for filter_ in self.filters:

            stands = np.genfromtxt(self.path_in_standards + filter_ + '/' + filter_ + '_standards_names.txt', dtype=None)

            if self.verbose:
                print stands

            self.Standards[filter_] = {}
            self.Standards[filter_]['file'] = []
            self.Standards[filter_]['mag'] = []
            self.Standards[filter_]['merr'] = []
            self.Standards[filter_]['airmass'] = []
            self.Standards[filter_]['mean'] = []
        #
        #     # file name, magnitude and magnitude error found.
        #     # (This code can be adapted to include anything in mag files - see below for correct dtype and intertoools.islice)
        #     for file_ in enumerate(stands):
        #         Standards[filter_]['file'].append(file_[1])
        #         mag = glob.glob(path_in_standards + filter_ + '/' + file_[1])
        #         with open(mag[0]) as f_in:
        #             mag_line5 = np.genfromtxt(itertools.islice(f_in, 79, None, 5), autostrip=True,  # skip_header=75,
        #                                       dtype=[('RAPERT', '<f8'), ('SUM', '<f8'), ('AREA', '<f8'),
        #                                              ('FLUX', '<f8'),
        #                                              ('MAG', '<f8'), ('MERR', '<f8'), ('PIER', '<f8'),
        #                                              ('PERROR', '<f8')])
        #             Standards[filter_]['mag'].append(float(mag_line5['MAG']))
        #             Standards[filter_]['merr'].append(float(mag_line5['MERR']))
        #
        #     # adding airmass
        #     Standards_airmass = np.genfromtxt(path_in_standards + filter_ + '/' + filter_ + '_airmass.txt', dtype=None)
        #
        #     # adding mean value imstat
        #     Standards_mean = np.genfromtxt(path_in_standards + filter_ + '/' + filter_ + '_imstat_mean.txt', dtype=None)
        #     #     Standards_mean = QTable.read(path_in_standards+filter_+'/'+filter_+r'_imstat_mean.txt',format='ascii')
        #
        #     for i in range(len(Standards_airmass)):
        #         Standards[filter_]['airmass'].append(Standards_airmass[i])
        #         Standards[filter_]['mean'].append(Standards_mean[i])
        #
        #     SaveStandards = QTable()
        #     SaveStandards['file'] = Standards[filter_]['file']
        #     SaveStandards['mag'] = Standards[filter_]['mag']
        #     SaveStandards['merr'] = Standards[filter_]['merr']
        #     SaveStandards['airmass'] = Standards[filter_]['airmass']
        #     SaveStandards['mean'] = Standards[filter_]['mean']
        #     SaveStandards.write(path_out + cluster + '_' + str(t.date()) + '_Standard_Stars_' + filter_ + '.txt',
        #                         format='ascii')
        #     del SaveStandards
        #
        #     if verbose == 'yes':
        #         print 'Saved:', cluster + '_' + str(t.date()) + '_Standard_Stars_' + filter_ + '.txt'
        #         print ''
        #
        # # printing lengths
        # if verbose_absolute == 'yes':
        #     print 'b:', '#files', len(Standards['b']['file']), '# mags', len(
        #         Standards['b']['mag']), '# mag errors', len(Standards['b']['merr']), '# airmass', len(
        #         Standards['b']['airmass']), '# airmass', len(Standards['b']['mean'])
        #     print 'v:', '#files', len(Standards['v']['file']), '# mags', len(
        #         Standards['v']['mag']), '# mag errors', len(Standards['v']['merr']), '# airmass', len(
        #         Standards['v']['airmass']), '# airmass', len(Standards['v']['mean'])
        #     print 'y:', '#files', len(Standards['y']['file']), '# mags', len(
        #         Standards['y']['mag']), '# mag errors', len(Standards['y']['merr']), '# airmass', len(
        #         Standards['y']['airmass']), '# airmass', len(Standards['y']['mean'])
        #
        # if verbose == 'yes':
        #     print ''
        #     print Standards


            # All Information broadcast in mag files and #lines / dtype / intertools.isclice lines of code to read them in.
            # N IMAGE               XINIT     YINIT     ID    COORDS                 LID    \
            # U imagename           pixels    pixels    ##    filename               ##     \
            # 7 lines  (12448.0, 430.81, 465.01, 1.0, nan, 1.0, nan)
            # itertools.islice(f_in, 75, None, 5)   #5 is just to make sure it skips past the total # of lines. 75 is starting line. None means it doesn't end.
            #                               dtype=[('IMAGE', '<f8'),('XINIT', '<f8'), ('YINIT', '<f8'),('ID', '<f8'),('Nan1', '<f8'), ('COORDS', '<f8'),('Nan2', '<f8')])

            # N XCENTER    YCENTER    XSHIFT  YSHIFT  XERR    YERR            CIER CERROR   \
            # U pixels     pixels     pixels  pixels  pixels  pixels          ##   cerrors  \
            # 7 lines (431.386, 463.797, 0.576, -1.213, 0.001, 0.001, 107.0)
            # itertools.islice(f_in, 76, None, 5)
            #                               dtype=[('XCENTER', '<f8'), ('YCENTER', '<f8'),('XSHIFT', '<f8'),('YSHIFT', '<f8'),('XERR', '<f8'), ('YERR', '<f8'),('CIER', '<f8')])

            # N MSKY           STDEV          SSKEW          NSKY   NSREJ     SIER SERROR   \
            # U counts         counts         counts         npix   npix      ##   serrors  \
            # 6 lines (60.20601, 2.79533, 1.056994, 2795.0, 27.0, 0.0)
            # itertools.islice(f_in, 77, None, 5)
            #                               dtype=[('MSKY', '<f8'), ('STDEV', '<f8'),('SSKEW', '<f8'),('NSKY', '<f8'),('NSREJ', '<f8'), ('SIER', '<f8')])

            # N ITIME          XAIRMASS       IFILTER                OTIME                  \
            # U timeunit       number         name                   timeunit               \
            # 4 lines (1.0, nan, nan, nan)
            # itertools.islice(f_in, 78, None, 5)
            #                               dtype=[('ITIME', '<f8'), ('XAIRMASS', '<f8'),('IFILTER', '<f8'),('OTIME', '<f8')])

            # N RAPERT   SUM           AREA       FLUX          MAG    MERR   PIER PERROR   \
            # U scale    counts        pixels     counts        mag    mag    ##   perrors  \
            # 8 lines (14.6, 2385199.0, 669.9021, 2344867.0, 9.075, 0.001, 0.0, nan)
            # itertools.islice(f_in, 79, None, 5)
            #                               dtype=[('RAPERT', '<f8'), ('SUM', '<f8'),('AREA', '<f8'),('FLUX', '<f8'),('MAG', '<f8'), ('MERR', '<f8'),('PIER', '<f8'),('PERROR', '<f8')])

    def PlotCMD(self, col, mag):
        """
        Plots color-magnitude diagrams of both trilegal and standardized stars.

        Parameters:
        self: an OpenClusters object
        col:  a two-dimensional array containing the first magnitude m1 and second magnitude m2 within a color m1-m2 i.e. ['v','b']
        mag:  a string corresponding to a magnitude name, i.e. 'v'

        Returns:
        a color-magnitude diagram.
        """
        fig, ax = plt.subplots(1, 3)
        fig.set_size_inches(21, 6.5)
        fig.suptitle(self.cluster_title + ' Color Magnitude Diagrams', fontsize='16')

        def subplotCMD(n_, x_, y_, title, xlabel, ylabel, color, legend=None, setup=True):
            """
            Creates a subplot for the final CMD plot.

            Parameters:
            n_:     an integer corresponding to the subplot number
            x_:     a list of x values
            x_err_: a float corresponding to the error in x
            y_:     a list of y values
            y_err_: a float corresponding to the error in y
            title:  a string corresponding to the plot title
            xlabel: a string corresponding to the x-axis label
            ylabel: a string corresponding to the y-axis label
            color:  a single-letter string corresponding to the desired plot color
            legend: a string corresponding to the legend title (if a legend is desired)
            setup:  automatically True. If initial plot is already set up, and you wish to add more subplots, set to False

            Returns:
            a subplot.
            """
            if setup:
                ax[n_].set_title(self.cluster_title +' : '+ title)
                ax[n_].set_xlabel(xlabel)
                ax[n_].set_ylabel(ylabel)
                ax[n_].invert_yaxis()
            ax[n_].plot(x_, y_, marker='x', markersize='2', linestyle='', label=legend, c=color);
            if legend != None:
                ax[n_].legend(loc='center left', bbox_to_anchor=(1, .5), numpoints=1, markerscale=9, framealpha=1);

            #TODO x_err and y_err input

        color1,color2 = col

        color_func = lambda x,y : x - y                         # given two magnitudes, x and y, computes their combined color
        error_func = lambda x,y : np.sqrt(x**2 + y**2)          # computes the combined error of two values, x and y

        # x and y for standard CMD
        x = map(color_func, self.Standards[color1]['mag'], self.Standards[color2]['mag'])
        x_err = map(error_func, self.Standards[color1]['merr'], self.Standards[color2]['merr'])
        y = self.Standards[mag]['mag']
        y_err = self.Standards[mag]['merr']

        # x and y for TRILEGAL simulation CMD
        x_tri = map(color_func, self.Trilegal[color1]['mag'], self.Trilegal[color2]['mag'])
        x_tri_err = map(error_func, self.Trilegal[color1]['merr'], self.Trilegal[color2]['merr'])
        y_tri = self.Standards[mag]['mag']
        y_tri_err = self.Standards[mag]['merr']

        #color1 color2 mag
        x_label = color1 + ' - ' + color2
        y_label = mag

        # Final plot
        subplotCMD(0, x, y, 'Stars', x_label, y_label, 'g')
        subplotCMD(1, x_tri, y_tri, 'Model TRILEGAL Stars', x_label, y_label, 'b')
        subplotCMD(2, x, y, 'Stars & Model TRILEGAL Stars', x_label, y_label, 'g',legend=self.cluster_title)
        subplotCMD(2, x_tri, y_tri, '', x_label, y_label, 'b', legend='TRILEGAL', setup=False)

        fig.savefig(path_out + cluster + '_' + str(t.date()) + '_plot_colmag_tri_cluster.jpg', bbox_inches='tight')

    def PlotPositions(self, markers, colors, xlim=None,ylim=None):
        """
        Plots cluster spatial positions.

        Parameters:
        self: 	 an OpenClusters object
        markers: an array corresponding to the desired marker for each cluster. example: array containing . o
        colors:  an array corresponding to the desired color of each marker for each cluster
        xlim: 	 a two dimensional array corresponding to the x-axis limits (if limits are desired).
                Initially set to None. example: [0,400]
        ylim: 	 a two dimensional array corresponding to the y-axis limits (if limits are desired).
                Initially set to None. example: [0,400]
        """

        # Check that stars were matched correctly, by visual inspection.

        fig, ax = plt.subplots(2, 2)
        fig.set_size_inches(14, 6.5 * 2)
        fig.suptitle(self.cluster_title + ' Stars XCENTER vs. YCENTER', fontsize='16')

        def subplot_pos_setup(n_, m_, title, xlabel, ylabel, xlim, ylim):
            """
            subplot for the final position plot.

            Parameters:
            n_:     an integer corresponding to the x subplot number
            m_:     an integer corresponding to the y subplot number
            title:  a string corresponding to the plot title
            xlabel: a string corresponding to the x-axis label
            ylabel: a string corresponding to the y-axis label
            xlim: 	a two dimensional array corresponding to the x-axis limits, i.e. [0,400]
            ylim: 	a two dimensional array corresponding to the y-axis limits, i.e. [0,400]
            """

            ax[n_, m_].set_title(title)
            ax[n_, m_].set_xlabel(xlabel)
            ax[n_, m_].set_ylabel(ylabel)
            ax[n_, m_].set_xlim(xlim[0], xlim[1])
            ax[n_, m_].set_ylim(ylim[0], ylim[1])

        def subplot_pos(n_, m_, x_, y_, marker, color, label=None):
            """
            Creates a subplot for the final position plot.

            Parameters:
            n_:     an integer corresponding to the x subplot number
            m_:     an integer corresponding to the y subplot number
            x_:     a list of x values
            x_err_: a float corresponding to the error in x
            y_:     a list of y values
            y_err_: a float corresponding to the error in y
            legend: a string corresponding to the legend title (if a legend is desired)

            Returns:
            plot and legend commands for a subplot.
            """

            ax[n_, m_].plot(x_, y_, marker=marker, markersize='2', linestyle='', c=color, label=label);
            if label != None:
                ax[n_, m_].legend(loc='center left', bbox_to_anchor=(1, .5), numpoints=1, markerscale=9, framealpha=1);

        xlabel,ylabel = ['x pixel position','y pixel position']

        if xlim == None:
            xlim = [self.Phot['xcenter'].min()-20,self.Phot['xcenter'].max() + 20]
        if ylim == None:
            ylim = [self.Phot['ycenter'].min()-20,self.Phot['ycenter'].max() + 20]
        #TODO self.Phot has xcenter and ycenter


        subplot_pos_setup(0, 0, 'All Stars / Selection of Image', xlabel, ylabel, xlim, ylim)
        subplot_pos_setup(0, 1, 'Natched Stars / Selection of Image', xlabel, ylabel, xlim, ylim)
        subplot_pos_setup(1, 0, 'All Stars / Complete Image', xlabel, ylabel, xlim, ylim)
        subplot_pos_setup(1, 1, 'Matched Stars / Complete Image', xlabel, ylabel, xlim, ylim)

        plotvalue = [0, 1]

        # Check that the number of filters equals the number of markers
        # if len(self.filters) != len(self.markers):
        #raise ValueError
        # else:


        for plotval in plotvalue:
            # Iterate through each subplot and add data points

            ax[plotval, 0].set_xlabel('X PIX')
            ax[plotval, 0].set_ylabel('Y PIX')

            for filt_, marker_, color_ in zip(self.filters, markers, colors):
                #TODO x_ =
                #TODO y_ =
                subplot_pos(plotval, 0, x_, y_, marker_, color_, label= filt_ + ' filter')
                #TODO x_ =
                #TODO y_ =
                subplot_pos(plotval, 1, x_, y_, marker_, color_, label=filt_ + ' filter')

            # ax[plotval, 0].plot(als_b['XCENTER'], als_b['YCENTER'], marker='o', markersize='10', linestyle='',
            #                     alpha=.2, label='b filter');
            # ax[plotval, 0].plot(als_v['XCENTER'], als_v['YCENTER'], marker='x', markersize='10', linestyle='',
            #                     label='v filter');
            # ax[plotval, 0].plot(als_y['XCENTER'], als_y['YCENTER'], marker='+', markersize='10', linestyle='',
            #                     label='y filter');

            ax[plotval, 1].set_xlabel('X PIX')
            ax[plotval, 1].set_ylabel('Y PIX')
            ax[plotval, 1].plot(Phot['b_XCENTER'], Phot['b_YCENTER'], marker='o', markersize='10', linestyle='',
                                alpha=.2, label='b filter');
            ax[plotval, 1].plot(Phot['v_XCENTER'], Phot['v_YCENTER'], marker='x', markersize='10', linestyle='',
                                label='v filter');
            ax[plotval, 1].plot(Phot['y_XCENTER'], Phot['y_YCENTER'], marker='+', markersize='10', linestyle='',
                                label='y filter');
            ax[0, 1].legend(loc='center left', bbox_to_anchor=(1.1, .5), numpoints=1);
        fig.savefig(path_out + cluster + '_' + str(t.date()) + '_plot_selectstars.jpg', bbox_inches='tight')