# Copyright (c) 2012,2013,2014 Tobias Brink
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

ifeq ($(wildcard ../Makefile.KIM_Config),)
  $(error ../Makefile.KIM_Config does not exist.  Something is wrong with your KIM API package setup)
endif
include ../Makefile.KIM_Config

MODEL_DRIVER_NAME := Tersoff_LAMMPS__MD_077075034781_001
MODEL_DRIVER_KIM_FILE_TEMPLATE := Tersoff_LAMMPS.kim.tpl
MODEL_DRIVER_INIT_FUNCTION_NAME := model_driver_init

LOCALOBJ = pair_tersoff.o model_driver_Tersoff.o

LOCALCLEAN =

include $(KIM_DIR)/$(builddir)/Makefile.ModelDriver
