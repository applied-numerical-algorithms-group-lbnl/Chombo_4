inline void test_agg_stencil_alloc(EBDataLoc*& a_ptr, unsigned int a_size, int  &compt_0, int& compt_1)
{
  a_ptr = new EBDataLoc[a_size];
  compt_1 = 0;
  compt_0 = 0;
  for(int i = 0 ; i < a_size ; i++)
  {
    a_ptr[i].m_dataID = ((unsigned int)(i*3.15)) % 2;
    if(a_ptr[i].m_dataID == 1)
    {
      a_ptr[i].m_offset = compt_1;
      compt_1++;
    }
   
    if(a_ptr[i].m_dataID == 0)
    {
      a_ptr[i].m_offset = compt_0;
      compt_0++;
    }
  }  
}

template<typename T>
inline void test_agg_stencil_alloc(Proto::pairPtr<T>& a_data, unsigned int a_size)
{
  a_data.ptr[1] = new T[a_size];
  a_data.ptr[0] = new T[a_size];
}

inline void test_agg_stencil_alloc(uint64_t*& a_data, unsigned int a_size)
{
  a_data = new uint64_t[a_size];
}

inline void test_agg_stencil_fill_sten_size_start(uint64_t* a_sizes, uint64_t* a_start, unsigned int a_size)
{
  unsigned int n = 0;
  for(int i = 0 ; i < a_size ; i++)
  {
    a_start[i] = n;
    a_sizes[i] = 3*i+2;
    n += 3*i+2;
  }
}

inline void test_agg_stencil_fill_eb_stencil(pair_t<double>*& a_eb_stencil, uint64_t* a_sizes, int a_size, int a_s0, int a_s1)
{
  int acc = 0;
  for(int i = 0; i < a_size ; i++)
    acc += a_sizes[i];
 
  a_eb_stencil = new pair_t<double>[acc];
  for(int i = 0; i < acc ; i++)
  {
    EBDataLoc tmp;
    tmp.m_dataID =  i % 2;//((unsigned int)(6.21*i)) % 2;
    if(tmp.m_dataID == 0) tmp.m_offset = acc % a_s0;
    if(tmp.m_dataID == 1) tmp.m_offset = acc % a_s1;
    
    a_eb_stencil[i].first = tmp;
    a_eb_stencil[i].second = 1;
  }
}

inline void test_agg_stencil_fill_value(pairPtr<double>& a_data, double a_val0, double a_val1, unsigned int a_size)
{
  for(int i = 0 ; i < a_size ; i++)
  {
    a_data.ptr[0][i] = a_val0;
    a_data.ptr[1][i] = a_val1;
  }
}

template<typename T>
inline void test_agg_stencil_copy_to_gpu(T* a_host, T*& a_devi, unsigned int size)
{
  unsigned int nbBytes = size * sizeof(T);
  T* tmp;
  protoMalloc(MEMTYPE_DEFAULT, tmp, nbBytes);
  a_devi = tmp;
  protoMemcpy(MEMTYPE_DEFAULT, a_devi,a_host,nbBytes,protoMemcpyHostToDevice);
}

template<typename T>
inline void test_agg_stencil_copy_to_gpu_eb_stencil(T* a_host, T*& a_devi, uint64_t* sizes, unsigned int size)
{
  int acc = 0;
  for(int i = 0; i<size ; i++) acc += sizes[i];
  unsigned int nbBytes = acc * sizeof(T);
  T* tmp;
  protoMalloc(MEMTYPE_DEFAULT, tmp, nbBytes);
  a_devi = tmp;
  protoMemcpy(MEMTYPE_DEFAULT, a_devi,a_host,nbBytes,protoMemcpyHostToDevice);
}

template<typename T>
inline void test_agg_stencil_copy_to_gpu(pairPtr<T>& a_host, pairPtr<T>& a_devi, unsigned int size)
{
  unsigned int nbBytes = size * sizeof(T);
  T* tmp0; 
  protoMalloc(MEMTYPE_DEFAULT, tmp0, nbBytes);
  protoMemcpy(MEMTYPE_DEFAULT, tmp0,a_host.ptr[0],nbBytes,protoMemcpyHostToDevice);
  a_devi.ptr[0] = tmp0;
  T* tmp1;
  protoMalloc(MEMTYPE_DEFAULT, tmp1, nbBytes);
  protoMemcpy(MEMTYPE_DEFAULT, tmp1,a_host.ptr[1],nbBytes,protoMemcpyHostToDevice);
  a_devi.ptr[1] = tmp1;
  
}


void test_agg_stencil_print(pairPtr<double>& a_data, unsigned int a_size)
{
  for(int i =0 ; i < a_size ; i++)
    std::cout << "(" << a_data.ptr[0][i] << "," << a_data.ptr[1][i] << ") ";
  std::cout << std::endl;
}

void test_agg_stencil_get_back_data(pairPtr<double>& a_host, pairPtr<double>& a_devi, unsigned int a_size)
{
  a_host.ptr[0] = new double[a_size];
  a_host.ptr[1] = new double[a_size];

  protoMemcpy(MEMTYPE_DEFAULT, a_host.ptr[0],a_devi.ptr[0],a_size*sizeof(double),protoMemcpyDeviceToHost);
  protoMemcpy(MEMTYPE_DEFAULT, a_host.ptr[1],a_devi.ptr[1],a_size*sizeof(double),protoMemcpyDeviceToHost);
  protoDeviceSynchronizeGPU();
}

bool test_answer_kernel_only_using(pairPtr<double>& a_res, double a_val0, double a_val1, unsigned int a_size)
{
  int compt_0 = 0;
  int compt_1 = 0;
  for(int i = 0 ; i < a_size ; i++)
  {
    int id = ((unsigned int)(i*3.15)) % 2;
    if(id == 1)
    {
      int dep = 2*i +3*(i*(i-1)/2);
      double res = (int)((3*i+2)/2) *(a_val1 + a_val0);

      if(i%2==1)
      {
        if(dep%2==1) res += a_val1;
        else res += a_val0;
      }

      if(a_res.ptr[1][compt_1] != res)
      {
	std::cout << " a_res.ptr[1]["<<compt_1<<"]: " << a_res.ptr[1][compt_1] << " diff of " << res << std::endl;
	return false;
      }
      compt_1++;
    }

    if( id == 0)
    {
      int dep = 2*i +3*(i*(i-1)/2);
      double res = 0;
      res =  (int)((3*i+2)/2) *(a_val1 + a_val0);
      if(i%2==1)
      {
        if(dep%2==1) res += a_val1;
        else res += a_val0;
      }
      if(a_res.ptr[0][compt_0] != res) 
      {
	std::cout << "i: " << i << " and a_res.ptr[0]["<<compt_0<<"]: " << a_res.ptr[0][compt_0] << " diff of " << res << std::endl;
	return false;
      }
      compt_0++;
    }
  }

  for(int i = compt_0 ; i < a_size ; i++)
    if(a_res.ptr[0][i] != 0) return false;


  for(int i = compt_1 ; i < a_size ; i++)
    if(a_res.ptr[1][i] != 0) return false;

  return true;
}

bool test_answer_increment_only(pairPtr<double> a_res, double a_val0, double a_val1, double a_val2, double a_val3, unsigned int a_size, bool increment_only)
{
  int compt_0 = 0;
  int compt_1 = 0;
  for(int i = 0 ; i < a_size ; i++)
  {
    int id = ((unsigned int)(i*3.15)) % 2;
    if(id == 1)
    {
      int dep = 2*i +3*(i*(i-1)/2);
      double res = (int)((3*i+2)/2) *(a_val1 + a_val0) + a_val3*increment_only;

      if(i%2==1)
      {
        if(dep%2==1) res += a_val1;
        else res += a_val0;
      }

      if(a_res.ptr[1][compt_1] != res)
      {
	std::cout << " a_res.ptr[1]["<<compt_1<<"]: " << a_res.ptr[1][compt_1] << " diff of " << res << std::endl;
	return false;
      }
      compt_1++;
    }

    if( id == 0)
    {
      int dep = 2*i +3*(i*(i-1)/2);
      double res = a_val2*increment_only;
      res +=  (int)((3*i+2)/2) *(a_val1 + a_val0);
      if(i%2==1)
      {
        if(dep%2==1) res += a_val1;
        else res += a_val0;
      }
      if(a_res.ptr[0][compt_0] != res) 
      {
	std::cout << "i: " << i << " and a_res.ptr[0]["<<compt_0<<"]: " << a_res.ptr[0][compt_0] << " diff of " << res << std::endl;
	return false;
      }
      compt_0++;
    }
  }

  for(int i = compt_0 ; i < a_size ; i++)
    if(a_res.ptr[0][i] != a_val2) return false;

  for(int i = compt_1 ; i < a_size ; i++)
    if(a_res.ptr[1][i] != a_val3) return false;

  return true;
}
 
bool test_agg_stencil_answer_scale(pairPtr<double> a_res, double a_val0, double a_val1, unsigned int a_size, double a_scale)
{
  int compt_0 = 0;
  int compt_1 = 0;
  for(int i = 0 ; i < a_size ; i++)
  {
    int id = ((unsigned int)(i*3.15)) % 2;
    if(id == 1)
    {
      int dep = 2*i +3*(i*(i-1)/2);
      double res = (int)((3*i+2)/2) *(a_val1 + a_val0)*a_scale;

      if(i%2==1)
      {
        if(dep%2==1) res += a_val1*a_scale;
        else res += a_val0*a_scale;
      }

      if(a_res.ptr[1][compt_1] != res)
      {
	std::cout << " a_res.ptr[1]["<<compt_1<<"]: " << a_res.ptr[1][compt_1] << " diff of " << res << std::endl;
	return false;
      }
      compt_1++;
    }

    if( id == 0)
    {
      int dep = 2*i +3*(i*(i-1)/2);
      double res = 0;
      res +=  (int)((3*i+2)/2) * (a_val1 + a_val0) * a_scale;
      if(i%2==1)
      {
        if(dep%2==1) res += a_val1 * a_scale;
        else res += a_val0 * a_scale;
      }
      if(a_res.ptr[0][compt_0] != res) 
      {
	std::cout << "i: " << i << " and a_res.ptr[0]["<<compt_0<<"]: " << a_res.ptr[0][compt_0] << " diff of " << res << std::endl;
	return false;
      }
      compt_0++;
    }
  }

  for(int i = compt_0 ; i < a_size ; i++)
    if(a_res.ptr[0][i] != 0) return false;

  for(int i = compt_1 ; i < a_size ; i++)
    if(a_res.ptr[1][i] != 0) return false;

  return true;
}

bool test_agg_stencil_is_equal(pairPtr<double> a_gpu, double* a_cpu[2], unsigned int a_size)
{
  for(int j = 0 ; j < 2 ; j++)
    for(int i = 0 ; i < a_size ; i ++)
      if(a_gpu.ptr[j][i] != a_cpu[j][i]) return false;

  return true;
}


