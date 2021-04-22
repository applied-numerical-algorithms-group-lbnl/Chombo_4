#pragma once
#include <chrono>
#include <ctime>
#include <cmath>
class test_timer
{

#ifdef PROTO_CUDA
	public:
	test_timer()
	{
		cudaEventCreate(&start);
		cudaEventCreate(&stop);
	}

	void begin()
	{
		cudaEventRecord(start);
	}

	void end()
	{
		cudaEventRecord(stop);
	}

	float duration()
	{
		float milliseconds;
		protoDeviceSynchronizeGPU();
		cudaEventElapsedTime(&milliseconds, start, stop);
		return milliseconds;
	}

	private:
	cudaEvent_t start; 
        cudaEvent_t stop;
#else
	public:
	test_timer()
	{
	}

	void begin()
	{
		start = std::chrono::system_clock::now();
	}

	void end()
	{
		stop = std::chrono::system_clock::now();
	}

	float duration()
	{
		return  std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();;
	}

	private:
    	std::chrono::time_point<std::chrono::system_clock> stop;
    	std::chrono::time_point<std::chrono::system_clock> start;
#endif
};
