HDF5Handler::HDF5Handler()
{
    // define custom types for Proto::Point and Proto::Box
    char indices[] = "ijklmnopqrstuvwxyz";
    
    H5T_PROTO_POINT = H5Tcreate(H5T_COMPOUND, sizeof(Point));
    for (int ii = 0; ii < DIM; ii++)
    {
        std::string pvar = "intvect";
        pvar += indices[ii];
        H5Tinsert(H5T_PROTO_POINT, pvar.c_str(), HOFFSET(Point, m_tuple)+sizeof(int)*ii, H5T_NATIVE_INT); 
    }

    H5T_PROTO_BOX = H5Tcreate(H5T_COMPOUND, sizeof(Box));
    for (int ii = 0; ii < DIM; ii++)
    {
        std::string lowName = "lo_";
        lowName += indices[ii];
        H5Tinsert(H5T_PROTO_BOX, lowName.c_str(), HOFFSET(Box, m_low)+ii*sizeof(int), H5T_NATIVE_INT);
    }
    for (int ii = 0; ii < DIM; ii++)
    {
        std::string highName = "hi_";
        highName += indices[ii];
        H5Tinsert(H5T_PROTO_BOX, highName.c_str(), HOFFSET(Box, m_high)+ii*sizeof(int), H5T_NATIVE_INT);
    }

    H5T_PROTO_STRING = H5Tcopy(H5T_C_S1);
    H5Tset_size(H5T_PROTO_STRING, 100);
    //H5T_PROTO_STRING = H5Tcreate(H5T_STRING, H5T_VARIABLE);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR, typename... Args>
void HDF5Handler::readLevel(LevelBoxData<T, C, MEM, CTR>& a_data,
        std::string a_filename,
        Args... a_params)
{
    PR_TIME("HDF5Handler::readLevel");
    // PARSE THE FILE NAME
    char fname[100];
    if (a_filename.substr(a_filename.find_last_of(".") + 1) != "hdf5")
    {
        a_filename += ".hdf5";
    }
    
    // The following pragmas suppress some GCC warnings that are overly conservative
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat-security"
    sprintf(fname, a_filename.c_str(), a_params...);
#pragma GCC diagnostic pop

#ifdef PR_MPI
    auto mpi_comm = MPI_COMM_WORLD;
    auto mpi_info = MPI_INFO_NULL;

    auto p_access = H5Pcreate(H5P_FILE_ACCESS);
    assert(H5Pset_fapl_mpio(p_access, mpi_comm, mpi_info) >= 0);
    auto file = H5Fopen(fname, H5F_ACC_RDONLY, p_access);
    assert(H5Pclose(p_access) >= 0);
#else
    auto file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
#endif
    hid_t H5T_T;
    getH5DataType<T>(&H5T_T);
    unsigned int numPatches = a_data.layout().size();
    unsigned int patchSize = a_data.patchSize();
    unsigned int numData = numPatches*patchSize;

    unsigned int numPatches_local = a_data.layout().localSize();
    unsigned int numData_local = numPatches_local*patchSize;

    //T   rawData[numData_local];
    //Box boxData[numPatches_local];
    
    T* rawData = (T*)malloc(numData_local*sizeof(T));
    Box* boxData = (Box*)malloc(numPatches_local*sizeof(Box));

    auto ds_boxes = H5Dopen2(file, "/level_0/boxes", H5P_DEFAULT);
    auto ds_data = H5Dopen2(file, "/level_0/data:datatype=0", H5P_DEFAULT);
#ifdef PR_MPI
    hsize_t boxStart[] = {a_data.layout().offset()};
    hsize_t dataStart[] = {a_data.offset()};
    hsize_t stride[] = {1};
    hsize_t boxCount[] = {numPatches_local};
    hsize_t dataCount[] = {numData_local};
    
    auto slab_boxes = H5Dget_space(ds_boxes); 
    assert(H5Sselect_hyperslab(slab_boxes, H5S_SELECT_SET, boxStart,  stride, boxCount, NULL) >= 0);
    auto slab_data = H5Dget_space(ds_data); 
    assert(H5Sselect_hyperslab(slab_data,  H5S_SELECT_SET, dataStart, stride, dataCount, NULL) >= 0);
   
    auto s_boxes = H5Screate_simple(1, boxCount,  NULL);
    auto s_data =  H5Screate_simple(1, dataCount, NULL);
    
    assert(H5Dread(ds_boxes, H5T_PROTO_BOX, s_boxes, slab_boxes, H5P_DEFAULT, boxData) >= 0);
    assert(H5Dread(ds_data,  H5T_T,         s_data,  slab_data,  H5P_DEFAULT, rawData) >= 0);

    assert(H5Sclose(slab_boxes) >= 0);
    assert(H5Sclose(slab_data) >= 0);
    assert(H5Sclose(s_boxes) >= 0);
    assert(H5Sclose(s_data) >= 0);
#else
    assert(H5Dread(ds_boxes, H5T_PROTO_BOX, H5S_ALL, H5S_ALL, H5P_DEFAULT, boxData) >= 0);
    assert(H5Dread(ds_data, H5T_T, H5S_ALL, H5S_ALL, H5P_DEFAULT, rawData) >= 0);
#endif
    assert(H5Dclose(ds_boxes) >= 0);
    assert(H5Dclose(ds_data) >= 0);
    
    // READ GHOST SIZE
    Point ghost;
    auto g_level_meta = H5Gopen(file, "/level_0/data_attributes", H5P_DEFAULT);
    auto att_ghost = H5Aopen(g_level_meta, "ghost", H5P_DEFAULT);
    assert(H5Aread(att_ghost, H5T_PROTO_POINT, &ghost) >= 0);
    assert(H5Aclose(att_ghost) >= 0);
    assert(H5Gclose(g_level_meta) >= 0);

    // READ DOMAIN BOX
    Box domainBox;
    auto g_level = H5Gopen(file, "/level_0", H5P_DEFAULT);
    auto att_prob_domain = H5Aopen(g_level, "prob_domain", H5P_DEFAULT);
    assert(H5Aread(att_prob_domain, H5T_PROTO_BOX, &domainBox) >= 0);
    assert(H5Aclose(att_prob_domain) >= 0);

    // READ TIME VARIABLES
    auto att_time = H5Aopen(g_level, "time", H5P_DEFAULT);
    auto att_dt = H5Aopen(g_level, "dt", H5P_DEFAULT);
    assert(H5Aread(att_time, H5T_NATIVE_DOUBLE, &m_time) >= 0);
    assert(H5Aread(att_dt, H5T_NATIVE_DOUBLE, &m_dt) >= 0);
    assert(H5Aclose(att_time) >= 0);
    assert(H5Aclose(att_dt) >= 0);
    
    assert(H5Gclose(g_level) >= 0);
    
    // COMPUTE BOX SIZE (ASSUMING FIXED BOX SIZE)
    Point boxSize = boxData[0].high() - boxData[0].low() + Point::Ones();

    // BUILD DISJOINT BOX LAYOUT
    array<bool, DIM> periodicity;
    for (int ii = 0; ii < DIM; ii++) { periodicity[ii] = false; }
    ProblemDomain domain(domainBox, periodicity);
    DisjointBoxLayout layout(domain, boxSize);

    // DEFINE LEVEL DATA
    a_data.define(layout, ghost);
    a_data.linearIn(rawData);
   
    // EXCESS CLEANUP 
    free(rawData);
    free(boxData);

    barrier();
    assert(H5Fclose(file) >= 0);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR, typename... Args>
void HDF5Handler::writeLevel(
        const LevelBoxData<T, C, MEM, CTR>& a_data,
        std::string a_filename,
        Args... a_params)
{
    std::vector<std::string> varnames(C);
    for (int ii = 0; ii < C; ii++)
    {
        varnames[ii] = "var" + to_string(ii + 1);
    }
    std::array<double, DIM> dx;
    for (int ii = 0; ii < DIM; ii++) { dx[ii] = 1.0; }
    writeLevel(varnames, dx, a_data, a_filename, a_params...);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR, typename... Args>
void HDF5Handler::writeLevel(
        std::vector<std::string> a_varNames,
        const LevelBoxData<T, C, MEM, CTR>& a_data,
        std::string a_filename,
        Args... a_params)
{
    std::array<double, DIM> dx;
    for (int ii = 0; ii < DIM; ii++) { dx[ii] = 1.0; }
    writeLevel(a_varNames, dx, a_data, a_filename, a_params...);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR, typename... Args>
void HDF5Handler::writeLevel(
        std::array<double, DIM> a_dx,
        const LevelBoxData<T, C, MEM, CTR>& a_data,
        std::string a_filename,
        Args... a_params)
{
    std::vector<std::string> varnames(C);
    for (int ii = 0; ii < C; ii++)
    {
        varnames[ii] = "var" + to_string(ii + 1);
    }
    writeLevel(varnames, a_dx, a_data, a_filename, a_params...);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR, typename... Args>
void HDF5Handler::writeLevel(
        double a_dx,
        const LevelBoxData<T, C, MEM, CTR>& a_data,
        std::string a_filename,
        Args... a_params)
{
    std::vector<std::string> varnames(C);
    for (int ii = 0; ii < C; ii++)
    {
        varnames[ii] = "var" + to_string(ii + 1);
    }
    std::array<double, DIM> dx;
    for (int ii = 0; ii < DIM; ii++) {dx[ii] = a_dx;}
    writeLevel(varnames, dx, a_data, a_filename, a_params...);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR, typename... Args>
void HDF5Handler::writeLevel(
        std::vector<std::string> a_varNames,
        double a_dx,
        const LevelBoxData<T, C, MEM, CTR>& a_data,
        std::string a_filename,
        Args... a_params)
{
    std::array<double, DIM> dx;
    for (int ii = 0; ii < DIM; ii++) { dx[ii] = a_dx; }
    writeLevel(a_varNames, dx, a_data, a_filename, a_params...);
}

template<unsigned int C>
void HDF5Handler::writeHeader(
        hid_t* a_file,
        int a_numLevels,
        std::vector<std::string> a_varNames)
{
    auto s_scalar = H5Screate(H5S_SCALAR);
    auto g_root = H5Gopen(*a_file, "/", H5P_DEFAULT);
    
    // ATTRIBUTES
    auto att_filetype =       H5Acreate2(g_root, "filetype", H5T_PROTO_STRING, s_scalar, H5P_DEFAULT, H5P_DEFAULT);
    auto att_num_components = H5Acreate2(g_root, "num_components", H5T_NATIVE_INT, s_scalar, H5P_DEFAULT, H5P_DEFAULT);
    auto att_num_levels =     H5Acreate2(g_root, "num_levels", H5T_NATIVE_INT, s_scalar, H5P_DEFAULT, H5P_DEFAULT);

    std::string filetype = "VanillaAMRFileType";
    int numComponents = C;
     
    assert(H5Awrite(att_filetype, H5T_PROTO_STRING, filetype.c_str()) >= 0);
    assert(H5Awrite(att_num_components, H5T_NATIVE_INT, &numComponents) >= 0);
    assert(H5Awrite(att_num_levels, H5T_NATIVE_INT, &a_numLevels) >= 0);

    assert(H5Aclose(att_filetype) >= 0);
    assert(H5Aclose(att_num_components) >= 0);
    assert(H5Aclose(att_num_levels) >= 0);
    
    while (a_varNames.size() < C)
    {
        int index = a_varNames.size() + 1;
        a_varNames.push_back("var" + to_string(index));
    }
    
    for (int ii = 0; ii < C; ii++)
    {
        char compName[100];
        sprintf(compName, "component_%i", ii);

        auto att_component = H5Acreate2(g_root, compName, H5T_PROTO_STRING, s_scalar, H5P_DEFAULT, H5P_DEFAULT);
        std::string c = a_varNames[ii];
        assert(H5Awrite(att_component, H5T_PROTO_STRING, c.c_str()) >= 0);
        assert(H5Aclose(att_component) >= 0);
    }

    assert(H5Gclose(g_root) >= 0);
    
    // CHOMBO GROUP
    char g_chombo_name[] = "/Chombo_global";
    auto g_chombo = H5Gcreate(*a_file, g_chombo_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    auto att_SpaceDim = H5Acreate2(g_chombo, "SpaceDim", H5T_NATIVE_INT, s_scalar, H5P_DEFAULT, H5P_DEFAULT);
    auto att_testReal = H5Acreate2(g_chombo, "testReal", H5T_NATIVE_DOUBLE, s_scalar, H5P_DEFAULT, H5P_DEFAULT);

    int SpaceDim = DIM;
    assert(H5Awrite(att_SpaceDim, H5T_NATIVE_INT, &SpaceDim) >= 0);

    assert(H5Aclose(att_SpaceDim) >= 0);
    assert(H5Aclose(att_testReal) >= 0);
    assert(H5Gclose(g_chombo) >= 0);
    assert(H5Sclose(s_scalar) >= 0);
}

template<typename T, unsigned int C, Centering CTR, typename... Args>
void HDF5Handler::writeLevel(
        std::vector<std::string> a_varNames,
        array<double, DIM> a_dx,
        const LevelBoxData<T, C, DEVICE, CTR>& a_data,
        std::string a_filename,
        Args... a_params)
{
    LevelBoxData<T, C, HOST, CTR> tmp(a_data.layout(), a_data.ghost());
    // the loop is here because we don't want the automatic ghost cell
    // exchanging that LevelBoxData::copyTo does
    for (auto iter : a_data.layout())
    {
        auto& data_i = a_data[iter];
        auto& temp_i = tmp[iter];
        data_i.copyTo(temp_i);
    }
    //a_data.copyTo(tmp);
    writeLevel(a_varNames, a_dx, tmp, a_filename, a_params...);
}

template<typename T, unsigned int C, Centering CTR, typename... Args>
void HDF5Handler::writeLevel(
        std::vector<std::string> a_varNames,
        array<double, DIM> a_dx,
        const LevelBoxData<T, C, HOST, CTR>& a_data,
        std::string a_filename,
        Args... a_params)
{
    PR_TIME("HDF5Handler::writeLevel");
    char fname[100];
    if (a_filename.substr(a_filename.find_last_of(".") + 1) != "hdf5")
    {
        a_filename += ".hdf5";
    }
    
    // The following pragmas suppress some GCC warnings that are overly conservative
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat-security"
    sprintf(fname, a_filename.c_str(), a_params...);
#pragma GCC diagnostic pop
   
    // CREATE THE FILE
#ifdef PR_MPI
    auto mpi_comm = MPI_COMM_WORLD;
    auto mpi_info = MPI_INFO_NULL;

    auto p_access = H5Pcreate(H5P_FILE_ACCESS);
    assert(H5Pset_fapl_mpio(p_access, mpi_comm, mpi_info) >= 0);
    auto file = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, p_access);
    assert(H5Pclose(p_access) >= 0);
#else
    auto file = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
#endif
    writeHeader<C>(&file, 1, a_varNames);
    addLevel(&file, a_data, a_dx, Point::Ones(), 0);
    assert(H5Fclose(file) >= 0);
}

#ifdef PR_AMR
template<typename T, unsigned int C, MemType MEM, Centering CTR, typename... Args>
void HDF5Handler::writeAMRData(
        std::vector<std::string> a_varNames,
        double a_dx,
        const AMRData<T, C, MEM, CTR>& a_data,
        std::string a_filename,
        Args... a_params)
{
    std::array<double, DIM> dx;
    for (int ii = 0; ii < DIM; ii++) { dx[ii] = a_dx; }
    writeAMRData(a_varNames, dx, a_data, a_filename, a_params...);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR, typename... Args>
void HDF5Handler::writeAMRData(
        double a_dx,
        const AMRData<T, C, MEM, CTR>& a_data,
        std::string a_filename,
        Args... a_params)
{
    std::array<double, DIM> dx;
    for (int ii = 0; ii < DIM; ii++) { dx[ii] = a_dx; }
    std::vector<std::string> varnames(C);
    for (int ii = 0; ii < C; ii++)
    {
        varnames[ii] = "var" + to_string(ii + 1);
    }
    writeAMRData(varnames, dx, a_data, a_filename, a_params...);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR, typename... Args>
void HDF5Handler::writeAMRData(
        array<double, DIM> a_dx,
        const AMRData<T, C, MEM, CTR>& a_data,
        std::string a_filename,
        Args... a_params)
{
    std::vector<std::string> varnames(C);
    for (int ii = 0; ii < C; ii++)
    {
        varnames[ii] = "var" + to_string(ii + 1);
    }
    writeAMRData(varnames, a_dx, a_data, a_filename, a_params...);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR, typename... Args>
void HDF5Handler::writeAMRData(
        std::vector<std::string> a_varNames,
        array<double, DIM> a_dx,
        const AMRData<T, C, MEM, CTR>& a_data,
        std::string a_filename,
        Args... a_params)
{
    //pout() << "HDF5Handler::writeAMRData start" << std::endl;
    PR_TIME("HDF5Handler::writeLevel");
    char fname[100];
    if (a_filename.substr(a_filename.find_last_of(".") + 1) != "hdf5")
    {
        a_filename += ".hdf5";
    }
    
    // The following pragmas suppress some GCC warnings that are overly conservative
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat-security"
    sprintf(fname, a_filename.c_str(), a_params...);
#pragma GCC diagnostic pop
   
    // CREATE THE FILE
#ifdef PR_MPI
    auto mpi_comm = MPI_COMM_WORLD;
    auto mpi_info = MPI_INFO_NULL;

    auto p_access = H5Pcreate(H5P_FILE_ACCESS);
    assert(H5Pset_fapl_mpio(p_access, mpi_comm, mpi_info) >= 0);
    auto file = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, p_access);
    assert(H5Pclose(p_access) >= 0);
#else
    auto file = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
#endif
    writeHeader<C>(&file, a_data.numLevels(), a_varNames);
    auto dx = a_dx;
    for (int ii = 0; ii < a_data.numLevels(); ii++)
    {
        Point refRatio = Point::Ones();
        if (ii < a_data.numLevels()-1)
        {
            refRatio = a_data.grid().refRatio(ii);
        }
        //pout() << "HDF5Handler::writeAMRData -> addLevel " << ii << std::endl;
        addLevel(&file, a_data[ii], dx, refRatio, ii);
        //pout() << "HDF5Handler::writeAMRData <- addLevel " << ii << std::endl;
        for (int dir = 0; dir < DIM; dir++)
        {
            dx[dir] /= refRatio[dir];
        }
    }
    assert(H5Fclose(file) >= 0);
}
#endif

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void HDF5Handler::addLevel(hid_t* a_file,
        const LevelBoxData<T, C, MEM, CTR>& a_data,
        array<double, DIM>& a_dx0,              // dx maybe ought to be type T instead of double
        Point a_refRatio,
        int a_level)
{
    PR_TIMERS("HDF5Handler::addLevel");
    PR_TIMER("parallel_write", timer);
    hid_t H5T_T;
    getH5DataType<T>(&H5T_T);

    unsigned int numPatches = a_data.layout().size();
    unsigned int patchSize = a_data.patchSize();
    unsigned int numData = numPatches*patchSize;

    unsigned int numPatches_local = a_data.layout().localSize();
    unsigned int numData_local = numPatches_local*patchSize;

    auto s_scalar = H5Screate(H5S_SCALAR);
    hsize_t spaceDims[] = {DIM};
    auto s_dim = H5Screate_simple(1, spaceDims, NULL);
    
    // CREATE GROUP
    std::string g_level_name = "/level_" + to_string(a_level);
    auto g_level = H5Gcreate(*a_file, g_level_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // CREATE ATTRIBUTES
    auto att_dt   = H5Acreate2(g_level, "dt",   H5T_NATIVE_DOUBLE, s_scalar, H5P_DEFAULT, H5P_DEFAULT);
    auto att_time = H5Acreate2(g_level, "time", H5T_NATIVE_DOUBLE, s_scalar, H5P_DEFAULT, H5P_DEFAULT);
    auto att_dx   = H5Acreate2(g_level, "dx",   H5T_NATIVE_DOUBLE, s_dim, H5P_DEFAULT, H5P_DEFAULT);
    auto att_ref_ratio   = H5Acreate2(g_level, "ref_ratio",   H5T_NATIVE_INT, s_dim, H5P_DEFAULT, H5P_DEFAULT);
    auto att_prob_domain = H5Acreate2(g_level, "prob_domain", H5T_PROTO_BOX,  s_scalar, H5P_DEFAULT, H5P_DEFAULT); 

    double dx[DIM];
    int refRatio[DIM];
    for (int dir = 0; dir < DIM; dir++)
    {
        dx[dir] = a_dx0[dir];
        refRatio[dir] = a_refRatio[dir];
    }

    Box problemDomain = a_data.layout().domain().box();

    assert(H5Awrite(att_dt,   H5T_NATIVE_DOUBLE, &m_dt) >= 0);
    assert(H5Awrite(att_time, H5T_NATIVE_DOUBLE, &m_time) >= 0);
    assert(H5Awrite(att_dx,   H5T_NATIVE_DOUBLE, &dx) >= 0);
    assert(H5Awrite(att_ref_ratio,   H5T_NATIVE_INT, &refRatio) >= 0);
    assert(H5Awrite(att_prob_domain, H5T_PROTO_BOX,  &problemDomain) >= 0);

    assert(H5Aclose(att_dt) >= 0);
    assert(H5Aclose(att_time) >= 0);
    assert(H5Aclose(att_dx) >= 0);
    assert(H5Aclose(att_ref_ratio) >= 0);
    assert(H5Aclose(att_prob_domain) >= 0);
    
      
    // CREATE DATASPACES
    hsize_t procDims[] = {numProc()};
    hsize_t offsetDims[] = {numPatches + 1};
    hsize_t boxDims[] = {numPatches};
    hsize_t dataDims[] = {numData};

    hsize_t boxDims_local[] = {numPatches_local};
    hsize_t dataDims_local[] = {numData_local};

    auto s_processors = H5Screate_simple(1, procDims, NULL);
    auto s_offsets    = H5Screate_simple(1, offsetDims, NULL);
    auto s_boxes      = H5Screate_simple(1, boxDims, NULL);
    auto s_data       = H5Screate_simple(1, dataDims, NULL);
    
    auto s_boxes_local = H5Screate_simple(1, boxDims_local, NULL);
    auto s_data_local  = H5Screate_simple(1, dataDims_local, NULL);

    // CREATE DATASETS
    auto ds_processors = H5Dcreate2(*a_file, (g_level_name + "/processors").c_str(),     H5T_NATIVE_INT, s_processors,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    auto ds_offsets =    H5Dcreate2(*a_file, (g_level_name + "/data:offsets=0").c_str(), H5T_NATIVE_LONG, s_offsets,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    auto ds_boxes =      H5Dcreate2(*a_file, (g_level_name + "/boxes").c_str(),          H5T_PROTO_BOX, s_boxes,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    auto ds_data =    H5Dcreate2(*a_file, (g_level_name + "/data:datatype=0").c_str(), H5T_T, s_data,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    // GET RAW DATA
    long int* offsetData = (long int*)malloc((numPatches+1)*sizeof(long int));
    Box* boxData = (Box*)malloc(numPatches_local*sizeof(Box));
    T* rawData = (T*)malloc(numData_local*sizeof(T));
    
    offsetData[0] = 0;
    int offset = 0;
    for (int ii = 1; ii <= numPatches; ii++)
    {
        offset += patchSize;
        offsetData[ii] = offset;
    }
    int ii = 0;
    for (auto iter : a_data.layout())
    {
        boxData[ii] = a_data.layout()[iter];
        ii++;
    }
    
    a_data.linearOut(rawData);
    
    // WRITE DATA
    // FIXME: write processor data
    assert(H5Dwrite(ds_offsets, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, offsetData) >= 0);
    #ifdef PR_MPI
    hsize_t boxStart[] = {a_data.layout().offset()};
    hsize_t dataStart[] = {a_data.offset()};
    hsize_t stride[] = {1};
    
    auto slab_boxes = H5Dget_space(ds_boxes); 
    assert(H5Sselect_hyperslab(slab_boxes, H5S_SELECT_SET, boxStart,  stride, boxDims_local, NULL) >= 0);
    auto slab_data = H5Dget_space(ds_data); 
    assert(H5Sselect_hyperslab(slab_data,  H5S_SELECT_SET, dataStart, stride, dataDims_local, NULL) >= 0);

    PR_START(timer);
    assert(H5Dwrite(ds_boxes, H5T_PROTO_BOX, s_boxes_local, slab_boxes, H5P_DEFAULT, boxData) >= 0);
    assert(H5Dwrite(ds_data,  H5T_T,         s_data_local,  slab_data,  H5P_DEFAULT, rawData) >= 0);
    PR_STOP(timer);

    assert(H5Sclose(slab_boxes) >= 0);
    assert(H5Sclose(slab_data) >= 0);
    #else
    assert(H5Dwrite(ds_boxes, H5T_PROTO_BOX, H5S_ALL, H5S_ALL, H5P_DEFAULT, boxData) >= 0);
    assert(H5Dwrite(ds_data,  H5T_T,         H5S_ALL, H5S_ALL, H5P_DEFAULT, rawData) >= 0);
    #endif
    
    // CLEANUP
    assert(H5Dclose(ds_processors) >= 0);
    assert(H5Dclose(ds_boxes) >= 0);
    assert(H5Dclose(ds_offsets) >= 0);
    assert(H5Dclose(ds_data) >= 0);

    assert(H5Sclose(s_offsets) >= 0);
    assert(H5Sclose(s_processors) >= 0);
    assert(H5Sclose(s_boxes) >= 0);
    assert(H5Sclose(s_data) >= 0);
    assert(H5Sclose(s_boxes_local) >= 0);
    assert(H5Sclose(s_data_local) >= 0);

    // METADATA
    std::string g_level_meta_name = g_level_name + "/data_attributes";
    auto g_level_meta = H5Gcreate(*a_file, g_level_meta_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    auto att_comps = H5Acreate2(g_level_meta, "comps", H5T_NATIVE_INT,  s_scalar, H5P_DEFAULT, H5P_DEFAULT);
    auto att_ghost = H5Acreate2(g_level_meta, "ghost", H5T_PROTO_POINT, s_scalar, H5P_DEFAULT, H5P_DEFAULT);
    auto att_object_type = H5Acreate2(g_level_meta, "object_type", H5T_PROTO_STRING, s_scalar, H5P_DEFAULT, H5P_DEFAULT);
    auto att_outputGhost = H5Acreate2(g_level_meta, "outputGhost", H5T_PROTO_POINT,  s_scalar, H5P_DEFAULT, H5P_DEFAULT);
    
    int comps = C;
    Point ghost = a_data.ghost();
    Point outputGhost = a_data.ghost(); //FIXME No idea what this attribute is meant to be

    assert(H5Awrite(att_comps, H5T_NATIVE_INT, &comps) >= 0);
    assert(H5Awrite(att_ghost, H5T_PROTO_POINT, &ghost) >= 0);
    assert(H5Awrite(att_object_type, H5T_PROTO_STRING, "FArrayBox") >= 0);
    assert(H5Awrite(att_outputGhost, H5T_PROTO_POINT, &outputGhost) >= 0);

    assert(H5Aclose(att_comps) >= 0);
    assert(H5Aclose(att_ghost) >= 0);
    assert(H5Aclose(att_object_type) >= 0);
    assert(H5Aclose(att_outputGhost) >= 0);
    // FINAL CLEAN UP
    assert(H5Gclose(g_level_meta) >= 0);
    assert(H5Gclose(g_level) >= 0);
    assert(H5Sclose(s_scalar) >= 0);

    free(offsetData);
    free(boxData);
    free(rawData);
}

template<typename T, unsigned int C, MemType MEM>
void HDF5Handler::addPatch(hid_t* a_file,
        const BoxData<T, C, MEM>& a_data,
        array<double, DIM>& a_dx0,
        Point a_ghost)              // dx maybe ought to be type T instead of double
{
    PR_TIMERS("HDF5Handler::addPatch");
    hid_t H5T_T;
    getH5DataType<T>(&H5T_T);

    unsigned int numPatches = 1;
    unsigned int patchSize = a_data.size();
    unsigned int numData = numPatches*patchSize;

    unsigned int numPatches_local = numPatches;
    unsigned int numData_local = numData;

    auto s_scalar = H5Screate(H5S_SCALAR);
    hsize_t spaceDims[] = {DIM};
    auto s_dim = H5Screate_simple(1, spaceDims, NULL);
    
    // CREATE GROUP
    std::string g_level_name = "/level_0";
    auto g_level = H5Gcreate(*a_file, g_level_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // CREATE ATTRIBUTES
    auto att_dt   = H5Acreate2(g_level, "dt",   H5T_NATIVE_DOUBLE, s_scalar, H5P_DEFAULT, H5P_DEFAULT);
    auto att_time = H5Acreate2(g_level, "time", H5T_NATIVE_DOUBLE, s_scalar, H5P_DEFAULT, H5P_DEFAULT);
    auto att_dx   = H5Acreate2(g_level, "dx",   H5T_NATIVE_DOUBLE, s_dim, H5P_DEFAULT, H5P_DEFAULT);
    auto att_ref_ratio   = H5Acreate2(g_level, "ref_ratio",   H5T_NATIVE_INT, s_dim, H5P_DEFAULT, H5P_DEFAULT);
    auto att_prob_domain = H5Acreate2(g_level, "prob_domain", H5T_PROTO_BOX,  s_scalar, H5P_DEFAULT, H5P_DEFAULT); 

    double dx[DIM];
    int refRatio[DIM];
    for (int dir = 0; dir < DIM; dir++)
    {
        dx[dir] = a_dx0[dir];
        refRatio[dir] = 1;
    }

    Box problemDomain = a_data.box();

    assert(H5Awrite(att_dt,   H5T_NATIVE_DOUBLE, &m_dt) >= 0);
    assert(H5Awrite(att_time, H5T_NATIVE_DOUBLE, &m_time) >= 0);
    assert(H5Awrite(att_dx,   H5T_NATIVE_DOUBLE, &dx) >= 0);
    assert(H5Awrite(att_ref_ratio,   H5T_NATIVE_INT, &refRatio) >= 0);
    assert(H5Awrite(att_prob_domain, H5T_PROTO_BOX,  &problemDomain) >= 0);

    assert(H5Aclose(att_dt) >= 0);
    assert(H5Aclose(att_time) >= 0);
    assert(H5Aclose(att_dx) >= 0);
    assert(H5Aclose(att_ref_ratio) >= 0);
    assert(H5Aclose(att_prob_domain) >= 0);
    
      
    // CREATE DATASPACES
    hsize_t procDims[] = {numProc()};
    hsize_t offsetDims[] = {numPatches + 1};
    hsize_t boxDims[] = {numPatches};
    hsize_t dataDims[] = {numData};

    hsize_t boxDims_local[] = {numPatches_local};
    hsize_t dataDims_local[] = {numData_local};

    auto s_processors = H5Screate_simple(1, procDims, NULL);
    auto s_offsets    = H5Screate_simple(1, offsetDims, NULL);
    auto s_boxes      = H5Screate_simple(1, boxDims, NULL);
    auto s_data       = H5Screate_simple(1, dataDims, NULL);
    
    auto s_boxes_local = H5Screate_simple(1, boxDims_local, NULL);
    auto s_data_local  = H5Screate_simple(1, dataDims_local, NULL);

    // CREATE DATASETS
    auto ds_processors = H5Dcreate2(*a_file, (g_level_name + "/processors").c_str(),     H5T_NATIVE_INT, s_processors,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    auto ds_offsets =    H5Dcreate2(*a_file, (g_level_name + "/data:offsets=0").c_str(), H5T_NATIVE_LONG, s_offsets,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    auto ds_boxes =      H5Dcreate2(*a_file, (g_level_name + "/boxes").c_str(),          H5T_PROTO_BOX, s_boxes,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    auto ds_data =    H5Dcreate2(*a_file, (g_level_name + "/data:datatype=0").c_str(), H5T_T, s_data,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    // GET RAW DATA
    long int* offsetData = (long int*)malloc((numPatches+1)*sizeof(long int));
    Box* boxData = (Box*)malloc(numPatches_local*sizeof(Box));
    T* rawData = (T*)malloc(numData_local*sizeof(T));
    
    offsetData[0] = 0;
    int offset = 0;
    for (int ii = 1; ii <= numPatches; ii++)
    {
        offset += patchSize;
        offsetData[ii] = offset;
    }
    boxData[0] = a_data.box();
    
    a_data.linearOut(rawData, a_data.box(), {0,C-1});

    // WRITE DATA
    assert(H5Dwrite(ds_offsets, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, offsetData) >= 0);
    assert(H5Dwrite(ds_boxes, H5T_PROTO_BOX, H5S_ALL, H5S_ALL, H5P_DEFAULT, boxData) >= 0);
    assert(H5Dwrite(ds_data,  H5T_T,         H5S_ALL, H5S_ALL, H5P_DEFAULT, rawData) >= 0);
    
    // CLEANUP
    assert(H5Dclose(ds_processors) >= 0);
    assert(H5Dclose(ds_boxes) >= 0);
    assert(H5Dclose(ds_offsets) >= 0);
    assert(H5Dclose(ds_data) >= 0);

    assert(H5Sclose(s_offsets) >= 0);
    assert(H5Sclose(s_processors) >= 0);
    assert(H5Sclose(s_boxes) >= 0);
    assert(H5Sclose(s_data) >= 0);
    assert(H5Sclose(s_boxes_local) >= 0);
    assert(H5Sclose(s_data_local) >= 0);

    // METADATA
    std::string g_level_meta_name = g_level_name + "/data_attributes";
    auto g_level_meta = H5Gcreate(*a_file, g_level_meta_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    auto att_comps = H5Acreate2(g_level_meta, "comps", H5T_NATIVE_INT,  s_scalar, H5P_DEFAULT, H5P_DEFAULT);
    auto att_ghost = H5Acreate2(g_level_meta, "ghost", H5T_PROTO_POINT, s_scalar, H5P_DEFAULT, H5P_DEFAULT);
    auto att_object_type = H5Acreate2(g_level_meta, "object_type", H5T_PROTO_STRING, s_scalar, H5P_DEFAULT, H5P_DEFAULT);
    auto att_outputGhost = H5Acreate2(g_level_meta, "outputGhost", H5T_PROTO_POINT,  s_scalar, H5P_DEFAULT, H5P_DEFAULT);
    
    int comps = C;
    Point ghost = a_ghost;
    Point outputGhost = a_ghost; //FIXME No idea what this attribute is meant to be

    assert(H5Awrite(att_comps, H5T_NATIVE_INT, &comps) >= 0);
    assert(H5Awrite(att_ghost, H5T_PROTO_POINT, &ghost) >= 0);
    assert(H5Awrite(att_object_type, H5T_PROTO_STRING, "FArrayBox") >= 0);
    assert(H5Awrite(att_outputGhost, H5T_PROTO_POINT, &outputGhost) >= 0);

    assert(H5Aclose(att_comps) >= 0);
    assert(H5Aclose(att_ghost) >= 0);
    assert(H5Aclose(att_object_type) >= 0);
    assert(H5Aclose(att_outputGhost) >= 0);
    // FINAL CLEAN UP
    assert(H5Gclose(g_level_meta) >= 0);
    assert(H5Gclose(g_level) >= 0);
    assert(H5Sclose(s_scalar) >= 0);

    free(offsetData);
    free(boxData);
    free(rawData);
}

#ifdef PR_AMR
template<typename T, unsigned int C, MemType MEM, typename... Args>
void HDF5Handler::writeFluxRegister(
        LevelFluxRegister<T, C, MEM>& a_register,
        std::string a_filename,
        Args... a_params)
{
    std::vector<std::string> varnames(C);
    for (int ii = 0; ii < C; ii++)
    {
        varnames[ii] = "var" + to_string(ii + 1);
    }
    std::array<double, DIM> dx;
    for (int ii = 0; ii < DIM; ii++) { dx[ii] = 1.0; }
    writeFluxRegister(varnames, dx, a_register, a_filename, a_params...);
}

template<typename T, unsigned int C, MemType MEM, typename... Args>
void HDF5Handler::writeFluxRegister(
        array<double, DIM> a_dx,
        LevelFluxRegister<T, C, MEM>& a_register,
        std::string a_filename,
        Args... a_params)
{
    std::vector<std::string> varnames(C);
    for (int ii = 0; ii < C; ii++)
    {
        varnames[ii] = "var" + to_string(ii + 1);
    }
    writeFluxRegister(varnames, a_dx, a_register, a_filename, a_params...);
}

template<typename T, unsigned int C, MemType MEM, typename... Args>
void HDF5Handler::writeFluxRegister(
        std::vector<std::string> a_varnames,
        array<double, DIM> a_dx,
        LevelFluxRegister<T, C, MEM>& a_register,
        std::string a_filename,
        Args... a_params)
{
    LevelBoxData<T, C, MEM> data(a_register.crseLayout(), Point::Zeros());
    data.setToZero();
    a_register.reflux(data, 1.0);
    writeLevel(a_varnames, a_dx, data, a_filename, a_params...);
}
#endif

template<typename T, unsigned int C, MemType MEM, typename... Args>
void HDF5Handler::writePatch(
        array<double, DIM> a_dx,
        const BoxData<T, C, MEM>& a_data,
        std::string a_filename,
        Args... a_params)
{
    std::vector<std::string> varnames(C);
    for (int ii = 0; ii < C; ii++)
    {
        varnames[ii] = "var" + to_string(ii + 1);
    }
    writePatch(varnames, a_dx, a_data, a_filename, a_params...);
}

template<typename T, unsigned int C, MemType MEM, typename... Args>
void HDF5Handler::writePatch(
        double a_dx,
        const BoxData<T, C, MEM>& a_data,
        std::string a_filename,
        Args... a_params)
{
    std::vector<std::string> varnames(C);
    for (int ii = 0; ii < C; ii++)
    {
        varnames[ii] = "var" + to_string(ii + 1);
    }
    std::array<double, DIM> dx;
    for (int ii = 0; ii < DIM; ii++) { dx[ii] = a_dx; }
    writePatch(varnames, dx, a_data, a_filename, a_params...);
}

template<typename T, unsigned int C, MemType MEM, typename... Args>
void HDF5Handler::writePatch(
        std::vector<std::string> a_varnames,
        const BoxData<T, C, MEM>& a_data,
        std::string a_filename,
        Args... a_params)
{
    std::array<double, DIM> dx;
    for (int ii = 0; ii < DIM; ii++) { dx[ii] = 1.0; }
    writePatch(a_varnames, dx, a_data, a_filename, a_params...);
}

template<typename T, unsigned int C, MemType MEM, typename... Args>
void HDF5Handler::writePatch(
        std::vector<std::string> a_varnames,
        double a_dx,
        const BoxData<T, C, MEM>& a_data,
        std::string a_filename,
        Args... a_params)
{
    std::array<double, DIM> dx;
    for (int ii = 0; ii < DIM; ii++) { dx[ii] = a_dx; }
    writePatch(a_varnames, dx, a_data, a_filename, a_params...);
}

template<typename T, unsigned int C, typename... Args>
void HDF5Handler::writePatch(
        std::vector<std::string> a_varnames,
        array<double, DIM> a_dx,
        const BoxData<T, C, HOST>& a_data,
        std::string a_filename,
        Args... a_params)
{
    PR_TIME("HDF5Handler::writePatch");
    char fname[100];
    if (a_filename.substr(a_filename.find_last_of(".") + 1) != "hdf5")
    {
        a_filename += ".hdf5";
    }
    
    // The following pragmas suppress some GCC warnings that are overly conservative
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat-security"
    sprintf(fname, a_filename.c_str(), a_params...);
#pragma GCC diagnostic pop
   
    // CREATE THE FILE
    auto file = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    writeHeader<C>(&file, 1, a_varnames);
    addPatch(&file, a_data, a_dx, Point::Zeros());
    assert(H5Fclose(file) >= 0);
}

template<typename T, unsigned int C, typename... Args>
void HDF5Handler::writePatch(
        std::vector<std::string> a_varnames,
        array<double, DIM> a_dx,
        const BoxData<T, C, DEVICE>& a_data,
        std::string a_filename,
        Args... a_params)
{
    BoxData<T, C, HOST> tmp(a_data.box());
    a_data.copyTo(tmp);
    writePatch(a_varnames, a_dx, tmp, a_filename, a_params...);
}  

template<>
inline void
HDF5Handler::getH5DataType<double>(hid_t* a_type){(*a_type) = H5T_NATIVE_DOUBLE;}

template<>
inline void
HDF5Handler::getH5DataType<int>(hid_t* a_type){(*a_type) = H5T_NATIVE_INT;}

template<>
inline void
HDF5Handler::getH5DataType<short>(hid_t* a_type){(*a_type) = H5T_NATIVE_SHORT;}

template<>
inline void
HDF5Handler::getH5DataType<char>(hid_t* a_type){(*a_type) = H5T_NATIVE_CHAR;}
