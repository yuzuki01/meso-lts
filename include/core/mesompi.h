#ifndef CORE_MESOMPI_H
#define CORE_MESOMPI_H


namespace MESO::Math {
    /// defined in core/math.h
    class Vector;
}

namespace MESO {
    /// defined in mesh/object.h
    template<class FieldType>
    class Field;
}

namespace MESO::MPI {

    const int main_rank = 0;
    extern int processor_num;
    extern int rank;

    struct MPI_TaskObject {
        int start = 0;
        int size = 0;
    };
    typedef std::vector<MPI_TaskObject> MPI_Task;
    MPI_Task get_task_distribution(int total_num);

    void Initialize(int *p_argc, char*** p_argv);

    void Finalize();

    /// gather field (for partitioned mesh)
    // coded in solver/field.cpp
    void GatherFieldList(std::vector<Field<Scalar>> &local, ScalarList &global, ObjectId column_id);

    /// broadcast
    void Bcast(bool &global);
    void Bcast(Scalar &global);
    void Bcast(MESO::Vector &global);
    void Bcast(ObjectIdList &global);

    /// reduce
    void AllReduce(Scalar local, Scalar &global);
    void AllReduce(const MESO::Vector &local, MESO::Vector &global);
    // coded in solver/field.cpp
    void AllReduce(ScalarList &local, ScalarList &global);
    void AllReduce(VectorList &local, VectorList &global);
    void AllReduce(Field<Scalar> &local, Field<Scalar> &global);
    void AllReduce(Field<Vector> &local, Field<Vector> &global);

    /// UDF MPI DataType
    namespace UDF {
        extern MPI_Datatype MPI_Vector;
        extern MPI_Op MPI_VectorSum;
        void vector_sum(void* invec, void* inoutvec, const int* len, MPI_Datatype* datatype);
        void MPI_UDF_VectorReduce();
    }
}

#endif //CORE_MESOMPI_H
