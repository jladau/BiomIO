# BiomIO
Lightweight, portable library for working with HDF5 BIOM files using Java.

The class BiomIO.java allows HDF5 BIOM files (http://biom-format.org) to be read using Java. The most recent NetCDF Java library (i.e., jar file; http://www.unidata.ucar.edu/software/thredds/current/netcdf-java) must be included to use the methods in this class, but no C libraries or other additional system libraries need to be installed. Hence, this class comprises a lightweight, portable lirbary for reading HDF5 BIOM files using Java. For additional documentation, see the 'doc/index.html' file. The data used for unit tests (tests are in the class 'BiomIOTest.java') is in the file 'data/rich_sparse_otu_table_hdf5.biom'.

Dependencies
------------

When developing applications using this library, the NetCDF-Java library version 4.6.3 or later must be included in the build path, which is available in the 'dependencies' directory or at https://github.com/Unidata/thredds/releases/tag/v4.6.3 (netcdfAll-4.6.3.jar). The NetCDF library can be bundled into compiled jar files (releases), so end-users do not need to include the NetCDF-Java library or install any additional dependencies to run applications that read HDF5 BIOM files using the BiomIO class.
