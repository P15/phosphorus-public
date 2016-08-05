import sys
import numpy as np
from scipy.stats import rankdata

#author: Roman Shraga
#
#This script goes through a list of allele frequencies per population and picks a set of AIMs for ancestry estimation
#
#use: python pick_aims.py


def pairwise_fst(p1,p2):
	if abs(p1 - p2) < .0001:
		return 0.0
	else:
		h_1 = 2.0*p1*(1-p1)
		h_2 = 2.0*p2*(1-p2)
		h_s = (h_1 + h_2) / 2

		p_t = (p1+p2)/2.0
		h_t = 2.0 * p_t * (1-p_t)

		fst = 0.0
		if h_t > 0:
			fst = (h_t - h_s) / h_t

		return fst

#set up population to continent dictionary
pop_geo_region  = {}
f2 = open("country_to_continent","r")
for line in f2:
	pts = line.strip().split("|")
	pop_geo_region[pts[1]] = pts[0]
f2.close()

geo_regions = ["African","Central Asian","East Asian","European","Middle Eastern","Native American","Native Oceanian","South Asian"]

#array to keep track of fst values
fst_arr = []
rsid_pop_dict = {}
rsid_geo_dict = {}
rsid_allele = {}

pop_arr = []
geo_arr = []
ind_geo = {}

rsid_order = []
pair_order = []
f = open("alfred_allele_freqs","r")
for line in f:
	parts = line.strip().split("|")
	rsid = parts[0]
	rsid_allele[rsid] = parts[1]

	#if we are on the first line
	if rsid == "rsid":
		for ind,pop in enumerate(parts[2:]):
			#if the population is in our set of included populations
			if pop in pop_geo_region and pop_geo_region[pop] in geo_regions:
				#list to keep track of continent order
				geo_arr.append(pop_geo_region[pop])
				#list to keep track of population order
				pop_arr.append(pop)
				#index to population dict
				ind_geo[ind] = pop_geo_region[pop]

		geo_arr = np.array(geo_arr)
		continue
	
	geo_dict = {}
	pop_freq_array = []
	for ind,freq in enumerate(parts[2:]):
		if ind in ind_geo:
			geo = ind_geo[ind]
			if not geo in geo_dict:
				geo_dict[geo] = [0,0]

			geo_dict[geo][0] += float(freq)
			geo_dict[geo][1] += 1

			pop_freq_array.append(float(freq))

	rsid_pop_dict[rsid] = pop_freq_array

	geo_freq_arr = []

	#average population frequencies to get continet level frequency
	for geo in geo_regions:
		geo_freq_arr.append(geo_dict[geo][0] / float(geo_dict[geo][1]))

	rsid_geo_dict[rsid] = geo_freq_arr
	
	#compute all pairwise fsts
	avg_reg = {}
	for i in range(len(geo_freq_arr)):
		for j in range(i+1,len(geo_freq_arr)):
			pair_fst_v = pairwise_fst(geo_freq_arr[i],geo_freq_arr[j])
			reg1 = geo_regions[i]
			reg2 = geo_regions[j]
			if reg1 > reg2:
				t = reg1
				reg1 = reg2
				reg2 = t
			regs =  "%s,%s" % (reg1,reg2)
			avg_reg[regs] = pair_fst_v

	mean = np.mean(geo_freq_arr)
	std = np.std(geo_freq_arr)

	#compute global fst
	fst = 0
	if mean != 0 and std != 0:
		fst = (std ** 2) / (mean * (1.0 - mean))

	#get pairwise fsts in order
	pair_fst_vals = [avg_reg[k] for k in sorted(avg_reg.keys())]

	#add global fst
	pair_fst_vals.insert(0,fst)
	
	#add rsid to full list
	fst_arr.append(pair_fst_vals)
	rsid_order.append(rsid)
	pair_order = sorted(avg_reg.keys())

pair_order.insert(0,'Global')
fst_arr = np.array(fst_arr).T

#compute rank per pair per rsid
ranks = [rankdata(i) for i in fst_arr]

#set to keep track of final rsids
rsid_set = set([])

#number of rsids to take per fst type
num_per_pair = {'Global':75,'African,East Asian':20, 'African,European':20, 'African,Native American':20, 'African,Native Oceanian':20, 'African,South Asian':20, 'East Asian,European':50, 'East Asian,Native American':100, 'East Asian,Native Oceanian':50, 'East Asian,South Asian':50, 'European,Native American':50, 'European,Native Oceanian':50, 'European,South Asian':350, 'Native American,Native Oceanian':50, 'Native American,South Asian':50, 'Native Oceanian,South Asian':50, 'European,Middle Eastern':150,'Middle Eastern,South Asian':150}

for i,rsid in enumerate(rsid_order):
	for j,rd in enumerate(ranks):
		if pair_order[j] in num_per_pair:
			inclusion_thresh = len(rsid_order) - num_per_pair[pair_order[j]]

			#if rsid rank is within inclusion threshold
			if rd[i] >= inclusion_thresh:
				rsid_set.add(rsid)

#write final allele freqs per continent
o = open("aim_allele_freqs","w")
geos_included = []
for geo in geo_regions:
	if geo not in ("Middle Eastern","Central Asian"):
		geos_included.append(geo)

o.write("rsid|allele|%s\n" % ("|".join(geos_included)))

for rsid in rsid_set:
	final_freq_arr = []
	for i,geo in enumerate(geo_regions):
		if geo not in ("Middle Eastern","Central Asian"):
			final_freq_arr.append("%.3f" % max(min(rsid_geo_dict[rsid][i],.99),.01))

	o.write("%s|%s|%s\n" % (rsid,rsid_allele[rsid],"|".join(final_freq_arr)))
o.close()

#write final allele freqs per population
o = open("aim_allele_freqs_pop","w")
o.write("rsid|allele|%s\n" % ("|".join(pop_arr)))
for rsid in rsid_set:
	final_freq_arr = []
	for i,pop in enumerate(pop_arr):
		final_freq_arr.append("%.3f" % max(min(rsid_pop_dict[rsid][i],.99),.01))

	o.write("%s|%s|%s\n" % (rsid,rsid_allele[rsid],"|".join(final_freq_arr)))
