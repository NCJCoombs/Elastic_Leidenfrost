# Name of executable
noinst_PROGRAMS = dynamic_EL_gnl_structured dynamic_EL_gnl_structured_steady

dynamic_EL_gnl_structured_SOURCES = dynamic_EL_gnl_structured.cc
dynamic_EL_gnl_structured_steady_SOURCES = dynamic_EL_gnl_structured_steady.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources
dynamic_EL_gnl_structured_LDADD = -L@libdir@ -laxisym_cylindrical_solid -lgeneric -lconstitutive $(EXTERNAL_LIBS) $(FLIBS)
dynamic_EL_gnl_structured_steady_LDADD = -L@libdir@ -laxisym_cylindrical_solid -lgeneric -lconstitutive $(EXTERNAL_LIBS) $(FLIBS)
			     
AM_CPPFLAGS += -I@includedir@
