#ifndef THREAD_NUM_H
#define THREAD_NUM_H

int get_thread_num_final();

inline int get_thread_num_final() {
  int thread_num;
#ifdef _OPENMP
  // multithreaded OpenMP version of code
  thread_num = omp_get_thread_num();
#else
  // single-threaded version of code
  thread_num = 0;
#endif
  return thread_num;
}

#endif
