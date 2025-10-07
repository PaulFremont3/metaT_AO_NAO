from matplotlib.backends.backend_pdf import PdfPages
import phate
import pandas as pd
import numpy as np
import scprep
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import sys
import seaborn as sns

def load_matrix(name, sep):
    with open(name, 'r') as f:
        l = [[np.float64(num) for num in line.rstrip().split(sep)] for line in f]
    f.close()
    return(l)
def load_vector(name, sep):
    with open(name, 'r') as f:
        line =f.readline()
        line =line.rstrip()
        l = [np.float64(num) for num in line.split(sep)]
    f.close()
    return(l)
def write_matrix(mat,name, sep):
    with open(name, 'w') as f:
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                f.write(str(mat[i,j]))
                f.write(sep)
            f.write('\n')
    f.close()
def write_vector(mat,name, sep):
    with open(name, 'w') as f:
        for i in mat:
            to_w = str(i)
            f.write(to_w)
            f.write(sep)
    f.close()

if __name__ == '__main__':
	n_cores=int(sys.argv[1])
	frac=str(sys.argv[2])
	zeros=str(sys.argv[3])
	subs=str(sys.argv[4])
	norm=str(sys.argv[5])
	sub2=str(sys.argv[6])
	knn=int(sys.argv[7])
	if subs == '0':
		subs = ''
	data = load_matrix('phate_training_expression_'+frac+'_'+zeros+subs+'_'+norm+'_'+sub2+'.txt', sep=' ')
	stats = load_vector('stations_'+frac+'.txt', sep=' ')
	unigenes = load_vector('unigenes_phate_'+frac+'_'+zeros+subs+'_'+norm+'_'+sub2+'.txt', sep=' ')
	df = pd.DataFrame(data, columns = stats, index=unigenes)
	
	phate_operator = phate.PHATE(n_jobs=n_cores, knn=knn)
	tree_phate = phate_operator.fit_transform(df)
	write_matrix(tree_phate, 'phate_fit_expression_'+str(knn)+'_'+frac+'_'+zeros+subs+'_'+norm+'_'+sub2+'.txt', sep=' ')
	
	if sub2 != '0':
		sils = list()
		clusters_list = list()	
		for i in range(5, 15):
			clusters = phate.cluster.kmeans(phate_operator, n_clusters=i, random_state=1)
			sil = phate.cluster.silhouette_score(phate_operator, n_clusters=i, random_state=1)
			sils.append(sil)
			clusters_list.append(clusters)
			print(i)
		min_sil=np.min(sils)
		ind_min=np.array(sils).argmin()
	
		write_vector(clusters_list[ind_min], 'clusters_phate_fit_expression_'+str(knn)+'_'+frac+'_'+zeros+subs+'_'+norm+'_'+sub2+'.txt', sep=' ')
	
		write_vector(sils, 'silhouette_phate_fit_expression_'+str(knn)+'_'+frac+'_'+zeros+subs+'_'+norm+'_'+sub2+'.txt', sep=' ')
	#pp = PdfPages('silhouette_phate_fit_expression_'+str(knn)+'_'+frac+'_'+zeros+subs+'_'+norm+'_'+sub2+'.pdf')
	#plt.figure(figsize=(8,8), dpi=300)
	#plt.scatter(np.array(range(2, 20)), np.array(sils))
	#plt.xlabel('Number of clusters', fontsize=17)
	#plt.ylabel('Mean silhouette', fontsize=17)
	#pp.savefig()
	#plt.close()
		

