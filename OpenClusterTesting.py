from OCMAP import OpenClusters as oc
import numpy as np
import glob

cluster = 'kron5'
cluster_title = 'Kronberger 5'
filters = [ 'b','v', 'y']
path_in_cluster = '/Users/lockepatton/Desktop/Research/Larson/OCMAP/OpenClusterKronberger5/'
path_in_standards = ''
path_out = '/Users/lockepatton/Desktop/Research/Larson/OCMAP/'

image_names_file = 'als_files.dat'
shifts = None
image_names = np.genfromtxt(path_in_cluster + image_names_file, dtype=None)

Kronberger5 = oc(cluster, cluster_title, [filters,image_names],
                 path_in_cluster, path_in_standards, path_out)

Kronberger5.PositionMatch(tol=1, image_names_file=image_names_file, shifts=shifts)


#
# image_names = np.genfromtxt(path_in_cluster + image_names_file, dtype=None)
#
# image_name = image_names[0]
# print image_name


#
#
# #testing inputs section
# print '----------------'
# #building readin command
# iraf_als_file = glob.glob( path_in_cluster + image_name)
#
# print "iraf_als_file:", iraf_als_file
#
# # opening iraf als photometry file
# with open(iraf_als_file[0]) as f_in:
#     # intertools.islice slices file to only obtain mag line
#     iraf_als = np.genfromtxt(itertools.islice(f_in, 0, None, 2),
#                              dtype=[('ID', '<f8'), ('XCENTER', '<f8'), ('YCENTER', '<f8'), ('MAG', '<f8'),
#                                     ('MERR', '<f8'), ('MSKY', '<f8'), ('NITER', '<f8')])