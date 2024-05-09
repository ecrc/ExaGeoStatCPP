#include <catch2/catch_all.hpp>

#ifdef USE_MPI
#include <mpi.h>
#endif

int main(int argc, char* argv[]) {

#ifdef USE_MPI
    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
#endif

    int result = Catch::Session().run( argc, argv );

#ifdef USE_MPI
    MPI_Finalize();
#endif

    return result;
}


