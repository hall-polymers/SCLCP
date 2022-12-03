# ------------------------------------------------------------
# soft link  fix_ext_alignment src files to LAMMPS_SRC directory
# ------------------------------------------------------------

# can set LAMMPS_SRC from command line using "make LAMMPS_SRC=MYPATH"

LAMMPS_SRC ?= '../lammps-stable/lammps-23Jun2022/src/'

files = fix_ext_alignment.cpp \
	fix_ext_alignment.h \
		
		
lammps_files := $(addprefix $(LAMMPS_SRC)/,$(files))
local_files := $(addprefix $(PWD)/,$(files))


help:
	@echo "Please use one of the following options:"
	@echo "  make link       to soft link the fix_ext_alignment C++ files into the LAMMPS src directory "
	@echo "  make unlink     to remove the soft links between the fix_ext_alignment and LAMMPS "
	@echo ""
	@echo "Note: Path to LAMMPS source is  LAMMPS_SRC=${LAMMPS_SRC}"
	@echo "This can be changed using 'make link LAMMPS_SRC=<new path>.'"


link:
	@#ln -sfr $(files) -t ${LAMMPS_SRC} 
	ln -sf $(local_files) -t ${LAMMPS_SRC}



unlink:
	rm ${lammps_files}

