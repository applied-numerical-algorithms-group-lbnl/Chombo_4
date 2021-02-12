class test_timer
{
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
		protoDeviceSynchronize();
		cudaEventElapsedTime(&milliseconds, start, stop);
		return milliseconds;
	}

	private:
	cudaEvent_t start; 
        cudaEvent_t stop;

};
