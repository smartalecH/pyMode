
bin_PROGRAMS = wgms3d

wgms3d_SOURCES = \
	main.cc \
	wgms3d.cc \
	wgms3d.h \
	simulation_parameters.h \
	diffops.cc \
	diffops.h \
	stencil.h \
	sparse.h \
	mgp.cc \
	mgp.h \
	geometry.h \
	bezier.cc \
	bezier.h \
	complex_functions.h \
	fortran_interface.h \
	standarddiffop_2nd.h \
	standarddiffop_handmade.h \
	pml.h \
	solver.h \
	make_unique.h

if HAS_SLEPC
wgms3d_SOURCES += solver_slepc.cc
endif

if HAS_ARPACK
wgms3d_SOURCES += solver_superlu_arpack.cc
endif

wgms3d_CXXFLAGS = -Wall -pedantic

DIST_SUBDIRS = tests matlab doc

# The following is a Maple script used to generate the
# homogeneous-space differential operators (standarddiffop_2nd.h).
EXTRA_DIST = \
	standarddiffop.mpl \
	NOTES
