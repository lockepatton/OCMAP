class OpenClusters:
    """
    Container for open clusters.
    """

    __author__ = ['Locke Patton', 'Ellis Avelone', 'Katie Crotts']

    def __init__(self, cluster, cluster_title, filters_images,
                 path_in_cluster, path_in_standards, path_out,
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

    def plotXY(self, x=[None, None], y=[None, None], save_fig=False):
        if self.verbose:
            print "\nRunning plotXY"

        import matplotlib.pyplot as plt
        from astropy.visualization import astropy_mpl_style
        plt.style.use(astropy_mpl_style)

        fig, ax = plt.subplots(1, 1)
        fig.set_size_inches(10, 10)

        for filter_ in self.filters:
            x_center = self.Centers[filter_]['XCENTER']
            y_center = self.Centers[filter_]['YCENTER']

            ax.plot(x_center, y_center, linestyle='None', marker='o', markersize=10, alpha=.4, label=filter_)

            ax.set_xlabel('X PIX')
            ax.set_ylabel('Y PIX')

            ax.set_xlim(x[0], x[1])
            ax.set_ylim(y[0], y[1])

            # ax.set_title(self.cluster_title + 'X-Center vs. Y-Center', fontsize='16',legend=filter_)
            ax.legend(title='Filters', fancybox=True, loc="upper left", bbox_to_anchor=(1, 1))

        fig.show()

        if save_fig:
            import os
            directory = self.path_out + 'plots/'
            if not os.path.exists(directory):
                os.makedirs(directory)

            filenameend_x = ''
            filenameend_y = ''
            if x != [None, None]:
                filenameend_x = '_x' + str(x[0]) + ':' + str(x[1])
            if y != [None, None]:
                filenameend_y = '_y' + str(y[0]) + ':' + str(y[1])

            filenameend = filenameend_x + filenameend_y

            fig.savefig(directory + self.cluster + '_xy_centers' + filenameend + '.jpg', bbox_inches='tight')
            del os

        del plt
        del astropy_mpl_style

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
        # Determining the shifts for each filter in each filter in x and y direction

        # default shifts are 0s
        if shifts == None:
            self.shifts = []
            for filter_ in self.filters:
                self.shifts.append([0., 0.])
        else:
            self.shifts = shifts

        self.StarMatch = {}
        self.StarMatch_extra = {}

        for it_, (filter_, image_name_) in enumerate(zip(self.filters, self.image_names)):

            self.Centers[filter_]['XCENTER_SHIFTED'] = self.Centers[filter_]['XCENTER'] + self.shifts[it_][0]
            self.Centers[filter_]['YCENTER_SHIFTED'] = self.Centers[filter_]['YCENTER'] + self.shifts[it_][1]

            # printing filters, image name, iraf_als_file if verbose
            if self.verbose:
                print "filter, image:", filter_, image_name_
                print 'x, y shift:', self.shifts[it_]

            # building matching framework
            self.StarMatch[filter_] = {}
            self.StarMatch[filter_]['ID'] = []
            self.StarMatch[filter_]['MAG'] = []
            self.StarMatch[filter_]['MERR'] = []
            self.StarMatch[filter_]['XCENTER'] = []
            self.StarMatch[filter_]['YCENTER'] = []
            self.StarMatch[filter_]['NSTARS'] = self.Centers[filter_]['NSTARS']
            self.StarMatch_extra[filter_] = {}

        # defining base filter as one with most stars
        magBase_index = np.argmax(self.nstars)
        magBase = self.filters[magBase_index]

        # defining filters 1 and 2 as first inputed into filters list (apart from the filter with max # of stars)
        mag1, mag2 = np.delete(self.filters, magBase_index)[:2]

        self.StarMatch_extra[mag1 + '_' + magBase + '_radius'] = []
        self.StarMatch_extra[mag2 + '_' + magBase + '_radius'] = []

        if self.verbose:
            print 'base filter:', magBase
            print 'matched filters:', mag1, mag2

        id_base = range(self.StarMatch[magBase]['NSTARS'])
        # id_mag1 = range(self.StarMatch[mag1]['NSTARS']) # v_id
        # id_mag2 = range(self.StarMatch[mag2]['NSTARS']) # b_id

        x_cen_magBase = self.Centers[magBase]['XCENTER_SHIFTED']  # y_xcen
        x_cen_mag1 = self.Centers[mag1]['XCENTER_SHIFTED']  # v_xcen
        x_cen_mag2 = self.Centers[mag2]['XCENTER_SHIFTED']  # b_xcen

        y_cen_magBase = self.Centers[magBase]['YCENTER_SHIFTED']  # y_ycen
        y_cen_mag1 = self.Centers[mag1]['YCENTER_SHIFTED']  # v_ycen
        y_cen_mag2 = self.Centers[mag2]['YCENTER_SHIFTED']  # b_ycen

        mag_magBase = self.Centers[magBase]['MAG']
        mag_mag1 = self.Centers[mag1]['MAG']
        mag_mag2 = self.Centers[mag2]['MAG']

        merr_magBase = self.Centers[magBase]['MERR']
        merr_mag1 = self.Centers[mag1]['MERR']
        merr_mag2 = self.Centers[mag2]['MERR']

        radius2 = lambda x, y: (x ** 2 + y ** 2)

        # TODO Go through **2 mistake and re-run all open cluster membership without the mistake

        ratio = tol

        if n_iterations == None:
            n_iterations = len(id_base)

        for base_i in id_base[:n_iterations]:
            x_cen_magBase_i = x_cen_magBase[base_i]
            y_cen_magBase_i = y_cen_magBase[base_i]

            mag1_magBase_radius2 = radius2((x_cen_mag1 - x_cen_magBase_i), (y_cen_mag1 - y_cen_magBase_i))
            mag1_i = np.argmin(mag1_magBase_radius2)

            if mag1_magBase_radius2.min() <= ratio * radius2(self.Centers[mag1]['MERR'][mag1_i],
                                                             self.Centers[magBase]['MERR'][base_i]):

                mag2_magBase_radius2 = radius2((x_cen_mag2 - x_cen_magBase_i), (y_cen_mag2 - y_cen_magBase_i))
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
            print '# Stars across filters', len(self.StarMatch[magBase]['ID']), '/', len(self.Centers[magBase]['ID'])

    def plot_StarMatchXY(self, x=[None, None], y=[None, None], save_fig=False) :

        if self.verbose:
            print "\nRunning plot_StarMatchXY"

        import matplotlib.pyplot as plt
        from astropy.visualization import astropy_mpl_style
        plt.style.use(astropy_mpl_style)

        fig, ax = plt.subplots(1, 1)
        fig.set_size_inches(10, 10)

        for filter_ in self.filters:
            x_center = self.StarMatch[filter_]['XCENTER']
            y_center = self.StarMatch[filter_]['YCENTER']

            ax.plot(x_center, y_center, linestyle='None', marker='o', markersize=10, alpha=.4, label=filter_)

            ax.set_xlabel('X PIX')
            ax.set_ylabel('Y PIX')

            ax.set_xlim(x[0], x[1])
            ax.set_ylim(y[0], y[1])

            # ax.set_title(self.cluster_title + 'X-Center vs. Y-Center', fontsize='16',legend=filter_)
            ax.legend(title='Filters', fancybox=True, loc="upper left", bbox_to_anchor=(1, 1))

        fig.show()

        if save_fig:
            import os
            directory = self.path_out + 'plots/'
            if not os.path.exists(directory):
                os.makedirs(directory)

            filenameend_x = ''
            filenameend_y = ''
            if x != [None, None]:
                filenameend_x = '_x' + str(x[0]) + ':' + str(x[1])
            if y != [None, None]:
                filenameend_y = '_y' + str(y[0]) + ':' + str(y[1])

            filenameend = filenameend_x + filenameend_y

            fig.savefig(directory + self.cluster + 'StarMatch_xy_centers' + filenameend + '.jpg', bbox_inches='tight')
            del os

        del plt
        del astropy_mpl_style