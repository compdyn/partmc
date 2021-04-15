def find_between( s, start, end):
    return s[s.find(start)+len(start):s.rfind(end)]

outfile = 'cbmz_species.json'
f_out = open(outfile, 'w')
f_out.write('{\n')
f_out.write('"pmc-data" : [\n')
with open('cbmz_mosaic.spc', 'r') as fp:
   line = fp.readline()
   cnt = 1
   while line:
       temp = line.split()
       if (len(temp) > 0):
           print('line', temp)
           spc = temp[0]
           f_out.write('  {\n')
           f_out.write('    "name" : "%s",\n' %spc)
           f_out.write('    "type" : "CHEM_SPEC",\n')
           string = find_between( line[4:-1], "{", "}")
           print('string',string)
           f_out.write('    "description" : "%s"\n' %(string))
           f_out.write('  },\n') 
       line = fp.readline()
f_out.write(']}\n')
f_out.close()

outfile = 'cbmz_abs_tol.json'
f_out = open(outfile, 'w')
f_out.write('{\n')
f_out.write('  "pmc-data" : [\n')
f_out.write('    {\n')
f_out.write('      "type" : "RELATIVE_TOLERANCE",\n')
f_out.write('      "value" : 1.0E-04\n')
f_out.write('    },\n')

with open('cbmz_mosaic.spc', 'r') as fp:
   line = fp.readline()
   cnt = 1
   while line:
       temp = line.split()
       if (len(temp) > 0):
           print('line', temp)
           spc = temp[0]
           f_out.write('    {\n')
           f_out.write('      "name" : "%s",\n' %spc)
           f_out.write('      "type" : "CHEM_SPEC",\n')
           f_out.write('      "absolute integration tolerance" : 1.0E-03\n')
           f_out.write('    },\n')
       line = fp.readline()
f_out.write(']}\n')
f_out.close()
