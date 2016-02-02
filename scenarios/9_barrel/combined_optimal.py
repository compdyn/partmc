
import sys, collections, operator

file_list = ["rmse_num_0908.dat", "rmse_num_0909.dat", "rmse_num_0925.dat"]
#file_list = ["rmse_num_0322.dat", "rmse_num_0908.dat", "rmse_num_0909.dat"]
d = collections.defaultdict(list)
for file in file_list:
	lines = open(file,'r').readlines()[5:]
	for line in lines:
		line_arr = line.strip().split()
		key = (float(line_arr[1]), float(line_arr[2]), float(line_arr[3]))
		value = float(line_arr[4])
		d[key].append(value)

res = {}
for key in d:
	li = d[key]
	if len(li) == len(file_list):
		res[key] = sum(li) / len(li)

sorted_res = sorted(res.items(), key=operator.itemgetter(1))
for key, value in sorted_res:
	print key, value
