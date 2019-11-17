# config to use the asdf files
HDF5=/opt/apps/intel18/impi18_0/phdf5/1.8.16/x86_64/lib
ASDF=/work/05880/tg851791/stampede2/asdf-library-1.0.0/asdf/lib/libasdf.a

./configure FC=ifort CC=icc CXX=icpc MPIFC=mpif90 VTK_LDFLAGS=-L$HDF5 ASDF_LIBS=$ASDF --with-asdf