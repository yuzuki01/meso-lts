#ifndef CORE_MESOMPI_H
#define CORE_MESOMPI_H


namespace MESO::Math {
    /// defined in core/math.h
    class Vector;
}

namespace MESO::MPI {

    const Label mainRank = 0;
    extern Label processorNum;
    extern Label rank;

    void Initialize(Label *p_argc, char ***p_argv);

    void Finalize();

    bool isMainProcessor();

    void Barrier();

    /// Costume MPI DataType
    namespace DataType {
        extern MPI_Datatype MPI_VECTOR;

        void MPI_Vector_init();

        void MPI_Vector_free();
    }

    /// Send
    template<typename MesoType>
    void SendWithStatus(const MesoType &var, const Label &target, const Label &tag);

    /// Recv
    template<typename MesoType>
    void RecvWithStatus(MesoType &var, const Label &target, Label &tag);

    /// Bcast
    template<typename MesoType>
    void Bcast(MesoType &global, const Label &root);
}

#endif //CORE_MESOMPI_H
