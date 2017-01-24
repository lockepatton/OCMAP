from OCMAP import OpenClusters as oc

cluster = 'kron5'
cluster_title = 'Kronberger 5'
filters = ['v', 'b', 'y']
path_in_cluster = ''
path_in_standards = ''
path_out = ''

Kronberger5 = oc(cluster, cluster_title, filters,
                 path_in_cluster, path_in_standards, path_out, verbose=1, verbose_absolute=1)

print Kronberger5

Kronberger5.Standardize()

Kronberger5.PositionMatch(tol=.2)