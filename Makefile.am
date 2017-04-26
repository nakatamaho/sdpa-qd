bin_PROGRAMS = sdpa_qd

if SYSTEM_MPACK
sdpa_qd_SOURCES = sdpa_chordal.cpp sdpa_dataset.cpp sdpa_io.cpp \
sdpa_jordan.cpp sdpa_linear.cpp sdpa_main.cpp sdpa_newton.cpp sdpa_parts.cpp \
sdpa_struct.cpp sdpa_tool.cpp sdpa_chordal.h sdpa_dataset.h \
sdpa_include.h sdpa_io.h sdpa_jordan.h sdpa_linear.h sdpa_newton.h \
sdpa_parts.h sdpa_right.h sdpa_struct.h sdpa_tool.h
else
sdpa_qd_SOURCES = sdpa_chordal.cpp sdpa_dataset.cpp sdpa_io.cpp \
sdpa_jordan.cpp sdpa_linear.cpp sdpa_main.cpp sdpa_newton.cpp sdpa_parts.cpp \
sdpa_struct.cpp sdpa_tool.cpp sdpa_chordal.h sdpa_dataset.h \
sdpa_include.h sdpa_io.h sdpa_jordan.h sdpa_linear.h sdpa_newton.h \
sdpa_parts.h sdpa_right.h sdpa_struct.h sdpa_tool.h \
mpack/Rdot.cpp  mpack/Rtrsm.cpp mpack/Rsyrk.cpp mpack/Raxpy.cpp \
mpack/Rgemm.cpp mpack/Rtrmm.cpp mpack/Rtrsv.cpp mpack/iMlaenv.cpp\
mpack/Rlamch.cpp mpack/Rlascl.cpp mpack/Rsytrd.cpp mpack/Rsterf.cpp\
mpack/Rorgtr.cpp mpack/Rlatrd.cpp mpack/Rsyr2k.cpp mpack/Rsytd2.cpp \
mpack/Rlanst.cpp mpack/Rlae2.cpp mpack/Rlapy2.cpp mpack/Rlasrt.cpp\
mpack/Rorgql.cpp mpack/Rorgqr.cpp mpack/Rsymv.cpp mpack/Rlarfg.cpp\
mpack/Rsyr2.cpp mpack/Rlassq.cpp mpack/Rorg2l.cpp mpack/Rlarft.cpp\
mpack/Rlarfb.cpp mpack/Rorg2r.cpp mpack/Rnrm2.cpp mpack/Rlarf.cpp \
mpack/Rger.cpp mpack/Rpotrf.cpp mpack/Mxerbla.cpp mpack/Rpotf2.cpp \
mpack/Mlsame.cpp mpack/Rscal.cpp mpack/Rcopy.cpp mpack/Rgemv.cpp \
mpack/Rtrmv.cpp mpack/Rsteqr.cpp mpack/Rlaset.cpp mpack/Rlaev2.cpp \
mpack/Rlasr.cpp mpack/Rlartg.cpp mpack/Rswap.cpp mpack/Rsyev.cpp\
mpack/Rlansy.cpp mpack/Mutils.cpp
endif

if ENABLE_SHARED
else
sdpa_qd_LDFLAGS = -static
endif

if SYSTEM_MPACK
sdpa_qd_LDADD = $(pthread_lib) -L$(SPOOLES_LIBDIR) -lspooles -L$(QD_LIBDIR) -lqd $(OPENMPLIB) -L$(MPACK_LIBDIR) -lmlapack_qd -lmblas_qd
sdpa_qd_CXXFLAGS= $(pthread_cflags) -I$(SPOOLES_INCLUDEDIR) -DUseMETIS=0 -I$(QD_INCLUDEDIR) $(OPENMPFLAGS) -I$(MPACK_INCLUDEDIR)
else
sdpa_qd_LDADD = $(pthread_lib) -L$(SPOOLES_LIBDIR) -lspooles -L$(QD_LIBDIR) -lqd $(OPENMPLIB)
sdpa_qd_CXXFLAGS= $(pthread_cflags) -I$(SPOOLES_INCLUDEDIR) -DUseMETIS=0 -I$(QD_INCLUDEDIR) $(OPENMPFLAGS) -I./mpack
endif

if BUILD_SPOOLES
SUBDIRS = spooles
endif
