#include <iostream>

#if defined(USE_OPENMP) && defined(USE_MPI)
#error "Cannot define both USE_OPENMP and USE_MPI at the same time"
#endif

#if !(defined(USE_OPENMP) || defined(USE_MPI))
#error "Either USE_OPENMP or USE_MPI must be defined"
#endif

#if defined(USE_OPENMP)
#include <omp.h>
#endif

#if defined(USE_MPI)
#include <mpi.h>
#endif



#if defined(USE_OPENMP)

int main()
{
  int threads = 0;
  int sum = 0;

#pragma omp parallel reduction(+ : threads)
  {
    threads += 1;

#pragma omp for reduction(+ : sum)
    for (int i = 1; i <= 1000; ++i)
      sum += i;
  }

  std::cout << "using " << threads << " threads, sum(1:1000) = " << sum << std::endl;

  return 0;
}

#endif



#if defined(USE_MPI)

int main(int argc, char** argv)
{
  int id;
  int procs;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);

  int sum = 0;
  for (int i = id; i <= 1000; i += procs) {
    sum += i;
  }

  int allsum;
  MPI_Reduce(&sum, &allsum, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD);

  if (id == 0) {
    std::cout << "using " << procs << " processes, sum(1:1000) = " << allsum << std::endl;
  }

  MPI_Finalize();

  return 0;
}

#endif
