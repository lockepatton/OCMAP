from OCMAP import OpenClusters as oc
import numpy as np

cluster = 'kron5'
cluster_title = 'Kronberger 5'
filters = [ 'b', 'v', 'y']
path_in_cluster = '/Users/lockepatton/Desktop/Research/Larson/OCMAP/OpenClusterKronberger5/'
path_in_standards = ''
path_out = '/Users/lockepatton/Desktop/Research/Larson/OCMAP/'

image_names_file = 'als_files.dat'
shifts = None
image_names = np.genfromtxt(path_in_cluster + image_names_file, dtype=None)

Kronberger5 = oc(cluster, cluster_title, [filters,image_names],
                 path_in_cluster, path_in_standards, path_out,verbose=False)

Kronberger5.PositionMatch(tol=3, n_iterations=None, shifts=shifts)
Kronberger5.PlotXY(save_fig=True,x=[0,500],y=[0,500])
Kronberger5.PlotMatchedXY(save_fig=True,x=[None,None],y=[None,None])


# print Kronberger5.Centers
# print '\n', Kronberger5.StarMatch, '\n', Kronberger5.StarMatch_extra