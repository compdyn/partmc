
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
f_out_rmse_num = open('combined_rmse.txt', 'w')
f_out_rmse_num.write("# Colume 1: caseID\n")
f_out_rmse_num.write("# Colume 2: prefactor\n")
f_out_rmse_num.write("# Colume 3: exponent\n")
f_out_rmse_num.write("# Colume 4: fractal dimension\n")
f_out_rmse_num.write("# Colume 5: root mean square error of number distribution\n")
case = 0
for key, value in sorted_res:
	print key, value
	case += 1
	f_out_rmse_num.write("%04d   %.3f   %.2f   %.1f    %.4f\n" % (case, key[0], key[1], key[2], value))
f_out_rmse_num.close()
