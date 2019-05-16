# CMake generated Testfile for 
# Source directory: /tmp/tmp.aSqUIr3xL5/examples/arkode/C_serial
# Build directory: /tmp/tmp.aSqUIr3xL5/cmake-build-debug/examples/arkode/C_serial
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(ark_analytic "/usr/bin/python" "/tmp/tmp.aSqUIr3xL5/test/testRunner" "--verbose" "--testname=ark_analytic" "--executablename=/tmp/tmp.aSqUIr3xL5/cmake-build-debug/examples/arkode/C_serial/ark_analytic" "--outputdir=/tmp/tmp.aSqUIr3xL5/cmake-build-debug/Testing/output" "--nodiff" "--answerdir=/tmp/tmp.aSqUIr3xL5/examples/arkode/C_serial" "--answerfile=ark_analytic.out")
add_test(ark_robertson "/usr/bin/python" "/tmp/tmp.aSqUIr3xL5/test/testRunner" "--verbose" "--testname=ark_robertson" "--executablename=/tmp/tmp.aSqUIr3xL5/cmake-build-debug/examples/arkode/C_serial/ark_robertson" "--outputdir=/tmp/tmp.aSqUIr3xL5/cmake-build-debug/Testing/output" "--nodiff" "--answerdir=/tmp/tmp.aSqUIr3xL5/examples/arkode/C_serial" "--answerfile=ark_robertson.out")
