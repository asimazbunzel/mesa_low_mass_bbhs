#################################################################

include $(MESA_DIR)/utils/makefile_header

# we do not use MESA default include and lib folders
# instead, we have local ones with our modules cc and ce
OWN_INCLUDE_DIR = ../include
OWN_LIB_DIR = ../lib

LOAD_OWN_MESA_BINARY = \
	-L$(MESA_LIB_DIR) -L$(OWN_LIB_DIR) -lcc -lce -lbinary $(STAR_LOAD_LIST)

LOAD_ALL = $(LOAD_EXTRAS1) $(LOAD_OWN_MESA_BINARY) $(LOAD_EXTRAS)

#################################################################

OWN_WORK_COMPILE = \
	$(FC) $(FCbasic) $(FCopenmp) $(FCchecks) $(FCdebug) $(FCfree) \
	-I$(MESA_INCLUDE_DIR) -I$(OWN_INCLUDE_DIR) $(INCLUDE_HDF5) -c

ifndef BINARY
	BINARY = bin2dco
endif

OBJS = $(BINARY_WORK_OBJS) bin2dco_utils.o bin2dco_io.o binary_run.o

WORK_DIR = ..
WORK_SRC_DIR = $(WORK_DIR)/src
STAR_JOB_DIR = $(MESA_DIR)/star/job
BINARY_JOB_DIR = $(MESA_DIR)/binary/job

$(BINARY) : $(OBJS)
	$(LOADER) $(FCopenmp) -o $(WORK_DIR)/$(BINARY) $(OBJS) $(LOAD_ALL)

#################################################################

# WORK_COMPILE defined in makefile_header

# STAR_WORK_OBJS

run_star_extras.o: $(WORK_SRC_DIR)/run_star_extras.f
	$(OWN_WORK_COMPILE) $<

run_star.o: $(STAR_JOB_DIR)/run_star.f
	$(OWN_WORK_COMPILE) $<


# BINARY_WORK_OBJS

binary_run.o: $(WORK_SRC_DIR)/binary_run.f
	$(OWN_WORK_COMPILE) $<

run_binary_extras.o: $(WORK_SRC_DIR)/run_binary_extras.f
	$(OWN_WORK_COMPILE) $<

run_binary.o: $(BINARY_JOB_DIR)/run_binary.f
	$(OWN_WORK_COMPILE) $<

# BIN2DCO_UTIL

bin2dco_utils.o: $(WORK_SRC_DIR)/bin2dco_utils.f
	$(OWN_WORK_COMPILE) $<

bin2dco_io.o: $(WORK_SRC_DIR)/bin2dco_io.f
	$(OWN_WORK_COMPILE) $<


clean:
	-@rm -f *.o *.mod $(WORK_DIR)/$(BINARY)

remk:
	-@rm -f run.o $(WORK_DIR)/$(BINARY)
