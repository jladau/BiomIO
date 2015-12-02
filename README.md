# BiomIO
Lightweight, portable library for working with HDF5 BIOM files using Java.

The class BiomIO.java allows HDF5 BIOM files (http://biom-format.org) to be read using Java. The most recent NetCDF Java library (i.e., jar file; http://www.unidata.ucar.edu/software/thredds/current/netcdf-java) must be included to use the methods in this class, but no C libraries or other additional system libraries need to be installed. Hence, this class comprises a lightweight, portable lirbary for reading HDF5 BIOM files using Java. For additional documentation, see the 'doc/index.html' file. The data used for unit tests (tests are in the class 'BiomIOTest.java') is in the file 'data/rich_sparse_otu_table_hdf5.biom'.
