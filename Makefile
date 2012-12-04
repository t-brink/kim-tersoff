# Copyright (c) 2012 Tobias Brink
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


MODEL_DRIVER_NAME := model_driver_Tersoff

LOCALOBJ = pair_tersoff.o

LOCALCLEAN =

# No changes should be needed below this line.

MODEL_DRIVER_BUILD_TARGET := $(strip $(MODEL_DRIVER_NAME)).a

ifndef KIM_API_DIR
   include $(KIM_DIR)/KIM_API/Include.mk
else
   include $(KIM_API_DIR)/Include.mk
endif


all: $(MODEL_DRIVER_BUILD_TARGET)

$(strip $(MODEL_DRIVER_NAME)).a: $(LOCALOBJ) $(strip $(MODEL_DRIVER_NAME)).o
	ar rcs $@ *.o

$(strip $(MODEL_DRIVER_NAME)).o: Makefile

Makefile: $(KIM_API_DIR)/Include.mk
	@touch Makefile

clean:
	rm -f *.o *.mod *.a *.so $(LOCALCLEAN)
