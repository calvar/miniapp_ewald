#include <stdio.h> 

int main() {
  int nDevices;

  cudaGetDeviceCount(&nDevices);
  printf("There are %d devices.\n", nDevices);
  for (int i = 0; i < nDevices; i++) {
    //printf("%d %d\n",i,nDevices);
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, i);
    printf("Device Number: %d\n", i);
    printf("  Device name: %s\n", prop.name);
    printf("  Compute capability: %d.%d\n", prop.major, prop.minor);
    printf("  Warp size: %d\n", prop.warpSize);
    printf("  Number of multiprocessors: %d\n", prop.multiProcessorCount);
    printf("  Max threads per multiprocessor: %d\n", prop.maxThreadsPerMultiProcessor);
    printf("  Max threads per block: %d\n", prop.maxThreadsPerBlock);
    printf("  Global memory size (bytes): %zu\n", prop.totalGlobalMem);
    printf("  Memory Clock Rate (KHz): %d\n",
           prop.memoryClockRate);
    printf("  Memory Bus Width (bits): %d\n",
           prop.memoryBusWidth);
    printf("  Peak Memory Bandwidth (GB/s): %f\n",
           2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
    printf("  L2 cache size (bytes): %d\n", prop.l2CacheSize);
    printf("  Shared memory per block (bytes): %zu\n\n", prop.sharedMemPerBlock);
  }

  return 0;
}