#!/usr/bin/python

import os, cgi
import cgitb; cgitb.enable()

print "Content-Type: text/html"     # HTML is following
print                               # blank line, end of headers
print "<html>"
print "<head>"
print "<title>Equilibrium Result</title>"
print "</head>"
print "<body>"

import cgi
form = cgi.FieldStorage()
temp = form.getfirst("temp")
rh = form.getfirst("rh")
salt = form.getfirst("salt")
dust = form.getfirst("dust")

print "<pre>"
cmd = "./equilib %s %s %s %s" % (temp, rh, salt, dust)
f = os.popen(cmd)
for line in f:
    print line,
f.close()
print "</pre>"

print "<br><br><a href=\"equilib.html\">Back to input page</a>"

print "</body>"
print "</html>"

