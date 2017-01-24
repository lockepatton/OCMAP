#from astropy.table import QTable, hstack
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
    def __init__(self, cluster, cluster_title, filters,
                 path_in_cluster, path_in_standards, path_out,
                 t=None, verbose=1, verbose_absolute=1):
        """
        Parameters
        ----------
        cluster : string
            Inside file names
        cluster_title : string
            For plot titles
        path_in_cluster : sting
            location of cluster .als files
        path_in_standards : string
            location of standard .mag files
        filters : array-like
            example is ['v', 'b', 'y']
        path_out : basestring
            Files and plots saved here
        verbose : Boolean
            print extra info?
        verbose_absolute : Boolean
            print important info?
        """

        self.cluster = cluster
        self.cluster_title = cluster_title
        self.filters = filters

        self.path_in_cluster = path_in_cluster
        self.path_in_standards = path_in_standards
        self.path_out = path_out

        self.verbose = verbose
        self.verbose_absolute = verbose_absolute

        #checking to see if user specified an output time, t, for use in naming output files
        if t != None:
            self.t = datetime.datetime.now()
        else:
            self.t = t

        if verbose_absolute:
            print 'output time seen in file names:', t.date()


    def PositionMatch(self, tol, shifts=None):
        """

        :param shifts: array with horizontal and vertical shift of filters
        :param tol: tolerance in matching files
        :return:
        """
        # UNDER CONSTRUCTION
        # Reading in Open Cluster Field Data Files
        # Turner11_dao_b.fits.als.1

        if shifts == None:
            shifts = [np.zeros(len(filters)),np.zeros(len(filters))]

        for filter_ in filters:
            #building readin command
            iraf_als_file = glob.glob(path_in_cluster + cluster + '_' + str(filter_)+'.fits.als.1')
            #TODO - make it the last als file. Or make input image new variable for PositionMatch

            #printing file name if verbose
            if verbose:
                print iraf_als_file

            self.Centers = {}

            #opening iraf als photometry file
            with open(iraf_als_file[0]) as f_in:

                #CONTINUTE HERE
                # RENAMED als_file to iraf_als_file
                # RENAMED als_b to iraf_als

                # intertools.islice slices file to only obtain mag line
                iraf_als = np.genfromtxt(itertools.islice(f_in, 0, None, 2),
                                      dtype=[('ID', '<f8'), ('XCENTER', '<f8'), ('YCENTER', '<f8'), ('MAG', '<f8'),
                                             ('MERR', '<f8'), ('MSKY', '<f8'), ('NITER', '<f8')])

                n_stars = len(iraf_als_file['ID'])
                if verbose_absolute:
                    print '# '+str(filter_)+' stars', n_stars

                self.Centers[filter_] = {}
                self.Centers[filter_]['ID'] = range(n_stars)
                self.Centers[filter_]['XCENTER'] = iraf_als['XCENTER'] + shifts[0][filter_]
                self.Centers[filter_]['YCENTER'] = iraf_als['YCENTER'] + shifts[1][filter_]
    #
    #
    #         v_id = range(len(als_v['ID']))
    #         v_xcen = iraf_als['XCENTER'] + v_y_dx
    #         v_ycen = iraf_als['YCENTER'] + v_y_dy
    #         b_id = range(len(als_b['ID']))
    #         b_xcen = als_b['XCENTER'] + b_y_dx
    #         b_ycen = als_b['YCENTER'] + b_y_dy
    #         y_xcen = als_y['XCENTER']
    #         y_ycen = als_y['YCENTER']
    #
    #     #
    #     #
    #     # als_b = glob.glob(path_in_cluster + cluster + '_b.fits.als.1')
    #     # als_v = glob.glob(path_in_cluster + cluster + '_v.fits.als.1')
    #     # als_y = glob.glob(path_in_cluster + cluster + '_y.fits.als.1')
    #     #
    #     # if verbose == 'yes':
    #     #     print als_b
    #     #     print als_v
    #     #     print als_y
    #     #
    #     # with open(als_b[0]) as f_in:
    #     #     als_b = np.genfromtxt(itertools.islice(f_in, 0, None, 2),
    #     #                           dtype=[('ID', '<f8'), ('XCENTER', '<f8'), ('YCENTER', '<f8'), ('MAG', '<f8'),
    #     #                                  ('MERR', '<f8'), ('MSKY', '<f8'), ('NITER', '<f8')])
    #     # with open(als_v[0]) as f_in:
    #     #     als_v = np.genfromtxt(itertools.islice(f_in, 0, None, 2),
    #     #                           dtype=[('ID', '<f8'), ('XCENTER', '<f8'), ('YCENTER', '<f8'), ('MAG', '<f8'),
    #     #                                  ('MERR', '<f8'), ('MSKY', '<f8'), ('NITER', '<f8')])
    #     # with open(als_y[0]) as f_in:
    #     #     als_y = np.genfromtxt(itertools.islice(f_in, 0, None, 2),
    #     #                           dtype=[('ID', '<f8'), ('XCENTER', '<f8'), ('YCENTER', '<f8'), ('MAG', '<f8'),
    #     #                                  ('MERR', '<f8'), ('MSKY', '<f8'), ('NITER', '<f8')])
    #     #
    #     # if verbose_absolute == 'yes':
    #     #     print '# b stars', len(als_b['XCENTER'])
    #     #     print '# v stars', len(als_v['XCENTER'])
    #     #     print '# y stars', len(als_y['XCENTER'])
    #
    #     # Shifts between files - (imalign values).
    #     v_y_dx = 0  # average value of xcenter in v - average value xcenter in y
    #     v_y_dy = 0  # average value of ycenter in v - average value ycenter in y
    #     b_y_dx = 0  # average value of xcenter in b - average value xcenter in y
    #     b_y_dy = 0  # average value of ycenter in b - average value ycenter in y
    #
    #     v_id = range(len(als_v['ID']))
    #     v_xcen = als_v['XCENTER'] + v_y_dx
    #     v_ycen = als_v['YCENTER'] + v_y_dy
    #     b_id = range(len(als_b['ID']))
    #     b_xcen = als_b['XCENTER'] + b_y_dx
    #     b_ycen = als_b['YCENTER'] + b_y_dy
    #     y_xcen = als_y['XCENTER']
    #     y_ycen = als_y['YCENTER']
    #
    #
    def Standardize(self):
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
            Plots both trilegal and standardized stars.
            :param self:
            :param col:
            :param mag:
            :return:
            """
            fig, ax = plt.subplots(1, 3)
            fig.set_size_inches(21, 6.5)
            fig.suptitle(self.cluster_title + ' Color Magnitude Diagrams', fontsize='16')

            def subplotCMD(n_, x_, x_err_, y_, y_err_, title, xlabel, y_label, color, legend=None, setup=True):
                if setup:
                    ax[n_].set_title(self.cluster_title +' : '+ title)
                    ax[n_].set_xlabel(xlabel)
                    ax[n_].set_ylabel(y_label)
                    ax[n_].invert_yaxis()
                ax[n_].plot(x_, y_, marker='x', markersize='2', linestyle='', label=legend, c=color);
                if legend != None:
                    ax[n_].legend(loc='center left', bbox_to_anchor=(1, .5), numpoints=1, markerscale=9, framealpha=1);

            color1,color2 = col

            color_func = lambda x,y : x - y
            error_func = lambda x,y : np.sqrt(x**2 + y**2)

            #x and y color mag
            x = map(color_func, self.Standards[color1]['mag'], self.Standards[color2]['mag'])
            x_err = map(error_func, self.Standards[color1]['merr'], self.Standards[color2]['merr'])
            y = self.Standards[mag]['mag']
            y_err = self.Standards[mag]['merr']

            #x and y color mag TRILEGAL simulation
            x_tri = map(color_func, self.Trilegal[color1]['mag'], self.Trilegal[color2]['mag'])
            x_tri_err = map(error_func, self.Trilegal[color1]['merr'], self.Trilegal[color2]['merr'])
            y_tri = self.Standards[mag]['mag']
            y_tri_err = self.Standards[mag]['merr']

            color1 color2 mag

            x_label = color1 + ' - ' + color2
            y_label = mag

            subplotCMD(0, x, x_err, y, y_err, 'Stars', x_label, y_label, 'g')
            subplotCMD(1, x_tri, x_tri_err, y_tri, y_tri_err, 'Model TRILEGAL Stars', x_label, y_label, 'b')
            subplotCMD(2, x, x_err, y, y_err, 'Stars & Model TRILEGAL Stars', x_label, y_label, 'g',legend=self.cluster_title)
            subplotCMD(2, x_tri, x_tri_err, y_tri, y_tri_err, '', x_label, y_label, 'b', legend='TRILEGAL', setup=False)

            fig.savefig(path_out + cluster + '_' + str(t.date()) + '_plot_colmag_tri_cluster.jpg', bbox_inches='tight')