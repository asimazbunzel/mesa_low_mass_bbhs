# This is the makefile for the package template

#################################################################
#
# PREAMBLE

include $(MESA_DIR)/utils/makefile_header

#################################################################
#
# SOURCE FILES

SRCS = \
	cc_def.f90 \
	cc_data.f90 \
	cc_mechanism.f90 \
	cc_ctrls_io.f90 \
	cc_support.f90 \
   cc_lib.f90

#################################################################
#
# TARGETS

LIB = libcc.a
LIB_OBJS = $(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SRCS)))
$(LIB) : $(LIB_OBJS)
	$(LIB_TOOL) $(LIB) $(LIB_OBJS)

#################################################################
#
# COMPILATION RULES

COMPILE = $(COMPILE_TO_DEPLOY) $(FCfree)
#COMPILE = $(COMPILE_TO_TEST) $(FCfree)

COMPILE_FIXED = $(COMPILE_TO_DEPLOY) $(FCfixed)
#COMPILE_FIXED = $(COMPILE_TO_TEST) $(FCfixed)

COMPILE_CMD = $(COMPILE)

%.o : %.f
	$(COMPILE_CMD) $<
%.o : %.f90
	$(COMPILE_CMD) $<
%.mod : %.o
	@true
%.o : %.mod

clean:
	-@rm -f *.o *.mod *.modmic *genmod.f90 *.a .depend

#################################################################
#
# DEPENDENCIES

# Paths to search for source & include files

VPATH = $(MOD_PUBLIC_DIR):$(MOD_PRIVATE_DIR)

# Automatic dependency generation using makedepf90

nodeps : $(.DEFAULT_GOAL)

NODEPS = $(or $(filter nodeps,$(MAKECMDGOALS)),$(filter clean,$(MAKECMDGOALS)))

ifeq ($(NODEPS),)

  DEPEND_DIRS = $(MOD_PUBLIC_DIR):$(MOD_PRIVATE_DIR)

  .depend :
	makedepf90 -m %m.mod -I$(DEPEND_DIRS) $(SRCS) > .depend
  -include .depend

endif

### END ###
