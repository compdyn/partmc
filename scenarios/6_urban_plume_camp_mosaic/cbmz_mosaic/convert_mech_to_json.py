import re

def find_between( s, start, end):
    return s[s.find(start)+len(start):s.rfind(end)]

def convert_to_float(s):

    line = re.sub('_dp', '', s)
    line = re.sub('D', 'e', line)
    print(line)
    val = line 
    return (val)

def hasFactor(inputString):
    hasFactor = False
    if (len(inputString.split(" ")) > 1):
       hasFactor= True
    return (hasFactor) 

def getFactor(inputString):
    pos = 1
    s = inputString
    for i in range(len(s)):
       if (s[i] == "."):
          pos+=1
       elif (s[i].isalpha()):
          break
       else:
          pos+=1
    return(s[0:pos-1], s[pos-1:])

outfile = 'cbmz_mechanism.json'
f_out = open(outfile, 'w')

f_out.write('{ "pmc-data" : [\n')
f_out.write('  {\n')
f_out.write('    "name" : "CBMZ",\n')
f_out.write('    "type" : "MECHANISM",\n')

f_out.write('    "reactions" : [\n')

with open('cbmz_mosaic_short.eqn', 'r') as fp:
   line = fp.readline()
   cnt = 1
   while line:
       eq_number = line.find("}")
       equal = line.find('=')
       reactants = line[eq_number+1:equal]
       colon = line.find(':', eq_number+1)
       products = line[equal+1:colon]
       notes = line[colon+1:-1]
       reacts = reactants.split(sep="+")
       prods = products.split(sep="+")
       if (len(line) > 0):
           is_photolysis = False
           f_out.write('      {\n')
           f_out.write('        "rxn id" : "R%i",\n' %cnt)
           f_out.write('        "reactants" : {\n')
           for r, spec in enumerate(reacts):
               spec = spec.strip()
               if (hasFactor(spec)):
                   factor, species_name = getFactor(spec)
                   f_out.write('          "%s" : {%03f}' %(species_name, float(factor)))
                   if (r < len(reacts)-1):
                       f_out.write(',\n')
                   else:
                       f_out.write('\n')
               else:
                   if (spec != 'hv'):
                       f_out.write('          "%s" : {}' %(spec.strip()))
                       if (r < len(reacts)-1):
                           f_out.write(',\n')
                       else:
                           f_out.write('\n')
                   else:
                       is_photolysis = True
           f_out.write('      },\n')
           f_out.write('        "products" : {\n')
           for p, spec in enumerate(prods):
               spec = spec.strip()
               if (hasFactor(spec)):
                   factor, species_name = getFactor(spec)
                   f_out.write('          "%s" : {"yield" : %04f}' %(species_name.strip(), float((factor))))
                   if (p < len(prods)-1):
                       f_out.write(',\n')
                   else:
                       f_out.write('\n')
               else:
                   f_out.write('          "%s" : {}' %(spec))
                   if (p < len(prods)-1):
                       f_out.write(',\n')
                   else:
                       f_out.write('\n')
           f_out.write('        },\n')
           # 
           if ("ARR3(" in notes):
              print('arrhenius rection', notes)
              start = notes.find('ARR3(')
              end = notes.find(')')
              tmp = notes[start+5:end]
              print(tmp)
              A = convert_to_float(tmp.split(',')[0])
              C = convert_to_float(tmp.split(',')[1]) 
              f_out.write('         "type" : "ARRHENIUS",\n')
              f_out.write('         "A" : %s,\n' %A)
              f_out.write('         "B" : 0.0,\n')
              f_out.write('         "C" : %s\n' %C)
           elif ("TROEMS" in notes):
              print('TROE reaction', notes)
              start = notes.find('TROEEMS(')
              end = notes.find(')')
              tmp = notes[start+7:end]
              f_out.write('        "type" : "TROE",\n')
              f_out.write('        "k0_A" : 3.0E-31,\n')
              f_out.write('        "k0_B" : -3.3,\n')
              f_out.write('        "k0_C" : -0.0E+00,\n')
              f_out.write('        "kinf_A" : 1.5E-12,\n')
              f_out.write('        "kinf_B" : 0.0E+00,\n')
              f_out.write('        "kinf_C" : -0.0E+00,\n')
              f_out.write('        "Fc" : 0.6,\n')
              f_out.write('        "N" : 1.0\n')
           elif ("TROEEMS" in notes):
              print('TROE reaction', notes)
              start = notes.find('TROEEMS(')
              end = notes.find(')')
              tmp = notes[start+7:end]
              f_out.write('        "type" : "TROE",\n')
              f_out.write('        "k0_A" : 3.0E-31,\n')
              f_out.write('        "k0_B" : -3.3,\n')
              f_out.write('        "k0_C" : -0.0E+00,\n')
              f_out.write('        "kinf_A" : 1.5E-12,\n')
              f_out.write('        "kinf_B" : 0.0E+00,\n')
              f_out.write('        "kinf_C" : -0.0E+00,\n')
              f_out.write('        "Fc" : 0.6,\n')
              f_out.write('        "N" : 1.0\n')
           elif (is_photolysis):
              print('photolysis reaction', notes)
              photo_rate =  notes[:-2]
              f_out.write('        "base rate" : %.4E,\n' %float(photo_rate))
              f_out.write('        "Fast-J id" : "%s",\n' % reacts[0].strip())
              f_out.write('        "type" : "PHOTOLYSIS"\n')
           elif('peroxy' in notes):
              print('peroxy reaction', notes)
           elif('RK_2HO2' in notes):
              print('RK_2HO2 reaction', notes)
           else:
              print('base rate?', notes)
           f_out.write('      },\n')
       line = fp.readline()
       cnt += 1
       if (len(line) == 0):
           f_out.write('      }\n')
           f_out.write('    ]\n')
           f_out.write('}\n')
           f_out.write(']}')
           f_out.close()
           quit()
