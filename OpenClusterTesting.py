from OCMAP import OpenClusters as oc

cluster = 'kron5'
cluster_title = 'Kronberger 5'
filters = ['v', 'b', 'y']
path_in_cluster = '/Users/lockepatton/Desktop/Research/Larson/OCMAP/OpenClusterKronberger5'
path_in_standards = ''
path_out = '/Users/lockepatton/Desktop/Research/Larson/OCMAP/'

image_names = ['test.als.1','test.als.2']
shifts = [0,1]

Kronberger5 = oc(cluster, cluster_title, filters,
                 path_in_cluster, path_in_standards, path_out)

Kronberger5.PositionMatch(tol=1, image_names=image_names, shifts=shifts)