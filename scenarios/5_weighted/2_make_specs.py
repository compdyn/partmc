#!/usr/bin/env python

for counter in ["10K", "100K"]:
    print "counter = ", counter
    filename_in = "weighted_template.spec"
    for run in ["flat", "wei-1", "wei-2", "wei-3"]: 
        filename_out = "spec/urban_plume_wc_%s_%s.spec" % (counter, run)
        print "filename_out ", filename_out
        f_in = open(filename_in, 'r')
        f_out = open(filename_out, 'w')

        for line in f_in:
            line = line.replace('%%OUTPUT_PREFIX%%', 'out/urban_plume_wc_%s_%s' % (counter, run))
            if run == "flat":
                line = line.replace('%%WEIGHTING_FUNC%%',  'none               ')
                if line.find('%%EXPONENT%%') > -1:
                    continue
                if line.find('%%REF_RADIUS%%') > -1:
                    continue
            else:
                line = line.replace('%%WEIGHTING_FUNC%%',  'power              ')
                line = line.replace('%%REF_RADIUS%%',  '1e-7              ')
                if run == "wei-1":
                    line = line.replace('%%EXPONENT%%',  '-1              ')
                elif run == "wei-2":
                    line = line.replace('%%EXPONENT%%',  '-2              ') 
                elif run == "wei-3":
                    line = line.replace('%%EXPONENT%%',  '-3              ')
            f_out.write(line)

        f_in.close()
        f_out.close()

#    filename_out_temp = "temp/temp_%02d.dat" % counter
#    print "filename_out_temp", filename_out_temp
#    f_out = open(filename_out_temp, 'w')

#    f_out.write("# time (s)\n")
#    f_out.write("# temp (K)\n")
#    f_out.write("time  %.2f %.2f\n" % ((counter-1) * 3600.0, (counter-1) * 3600.0 + 1200))
#    f_out.write("temp  290   280\n")   
#    f_out.close()

