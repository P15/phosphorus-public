import numpy as np
import sys
import argparse
import multiprocessing

#author: Roman Shraga
#
#This script predicts ethnicity admixture for a file containing genotype information.
#
#use: python predict_ancestry.py -t <num_threads> -i <rsid_info_file> -g <genotype_file> -o <output file name>


def compute_likelihood(g,f,q):
	l = 0
	for i,g_v in enumerate(g):
		l_g = 0
		l_cg = 0
		l_g +=  np.log(np.dot(f[i], q))
		l_cg += np.log(np.dot(1.0 - f[i], q))
		
		l += float(g_v) * l_g
		l += float(2-g_v) * l_cg
	return l


def determine_likeliest(genotypes,num_regions,rsid_info,rsid_order,sample,result_queue):
	#q initial value is uniform across all geos
	q = [float(1)/float(num_regions)] * num_regions
	g = []
	f = []

	valid = set(['A','T','G','C'])

	#set up genotype and frequency vectors
	for ind,v in enumerate(genotypes):
		rsid = rsid_order[ind]
		ref_allele = rsid_info[rsid]["allele"]

		if v[0] in valid and v[1] in valid:
			matches = 0
			for i in v:
				if i == ref_allele:
					matches += 1

			g.append(matches)
			f.append(rsid_info[rsid]["freqs"])

	q = np.array(q)
	g = np.array(g)
	f = np.array(f)

	q_n_1 = q

	e = .01
	l_n = -1.0 * sys.maxint
	l_n_1 = compute_likelihood(g,f,q)

	c = 0
	while l_n_1 - l_n > e:
		c += 1
		q = q_n_1
		q_n_1 = [0] * len(q)

		for i,g_v in enumerate(g):
			a_denom = np.dot(q,f[i])
			b_denom = np.dot(q,1.0 - f[i])

			a = np.multiply(f[i],q) / a_denom
			b = np.multiply(1.0 - f[i],q) / b_denom

			q_n_1 += float(g_v) * a
			q_n_1 += float(2 - g_v) * b

		q_n_1 = (float(1)/float(2*len(g))) * q_n_1
		l_n = l_n_1
		l_n_1 = compute_likelihood(g,f,q_n_1)

	print "Sample: %s, Iterations: %d, Likelihood: %f" % (sample,c,l_n_1)


	result_string = [str(i) for i in q_n_1]

	result_queue.put("%s|%s\n" % (sample,"|".join(result_string)))
	
	return


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("-g",action="store", dest="g")
	parser.add_argument("-o",action="store", dest="o")
	parser.add_argument("-i",action="store", dest="i")
	parser.add_argument("-t",action="store", dest="t")


	args = parser.parse_args()
	rsid_info_file = args.i
	geno_file = args.g
	out_file = args.o
	num_threads = int(args.t)


	#set up rsid information dictionary using provided rsid info file
	f = open(rsid_info_file,"r")
	rsid_info = {}
	geo_regions = []
	for line in f:
		parts = line.strip().split("|")
		if parts[0] == "rsid":
			geo_regions = parts[2:]
			continue

		rsid = parts[0]
		rsid_info[rsid] = {}
		rsid_info[rsid]["allele"] = parts[1]
		rsid_info[rsid]["freqs"] = [float(i) for i in parts[2:]]

	f.close()


	print "\n\nPredicting Ancestry\n"
	jobs = []
	result_queue = multiprocessing.Queue()

	rsid_order = []
	f2 = open(geno_file,"r")
	o = open(out_file,"w")
	o.write("sample|%s\n" % ("|".join(geo_regions)))

	sample_num = 0
	for line in f2:
		parts = line.strip().split(",")
		if parts[0] == 'sample':
			rsid_order = parts[1:]
		else:
			sample_num += 1
			sample = parts[0]

			genotypes = np.array(parts[1:])

			p = multiprocessing.Process(target=determine_likeliest, args=(genotypes,len(geo_regions),rsid_info,rsid_order,sample,result_queue,))
			jobs.append(p)
			p.start()
			if sample_num % num_threads == 0:
				for job in jobs:
					job.join()
				for job in jobs:
					o.write(result_queue.get())
				jobs = []

	f2.close()

	if len(jobs) > 0:
		for job in jobs:
			job.join()

		for job in jobs:
			o.write(result_queue.get())
	o.close()


	print "Finished Predicting Ethnicity"

