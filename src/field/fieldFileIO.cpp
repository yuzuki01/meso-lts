#include "field/field.h"


namespaceMESO


template<>
void volField<Scalar>::output(const String &fieldName) {
    const auto values = fvm::processorCommAllData(*this);
    if (MPI::rank != MPI::mainRank) return;

    const int DATA_PRECISION = 18;
    const int LINE_DATA_NUM = 30;
    Utils::mkdir(time_.name());
    StringStream sPath;
    sPath << time_.name() << "/" << fieldName << ".plt";
    std::fstream fp;
    fp.open(sPath.str(), std::ios::out | std::ios::trunc);
    if (!fp.is_open()) {
        logger.warn << "Cannot open file: " << sPath.str() << std::endl;
        fp.close();
        FATAL_ERROR_THROW;
    }
    // write head
    fp << "TITLE = \"" << fieldName << "\"\n";
    fp << "FileType = \"SOLUTION\"\n"
          "VARIABLES = \"" << fieldName << "\"\n";
    fp << "ZONE N=" << mesh_.nodeNum() << ", E=" << mesh_.cellNum();
    fp << ", VARLOCATION=([1]=CELLCENTERED), DATAPACKING=BLOCK, ZONETYPE=";
    fp << (mesh_.dimension() == 2 ? "FEQUADRILATERAL" : "FEBRICK") << ", SOLUTIONTIME=" << time_.name() << std::endl;
    // write values
    int count;
    count = 0;
    fp << std::endl << "## " << fieldName << std::endl;
    for (auto value: values) {
        fp << "\t" << std::setprecision(DATA_PRECISION) << value;
        if (count++ >= LINE_DATA_NUM) {
            fp << std::endl;
            count = 0;
        }
    }
    // close
    fp.close();
    logger.note << "Write volScalarField data <" << fieldName << "> -> " << sPath.str() << std::endl;
}


template<>
void volField<Vector>::output(const String &fieldName) {
    const auto values = fvm::processorCommAllData(*this);
    if (MPI::rank != MPI::mainRank) return;

    const int DATA_PRECISION = 18;
    const int LINE_DATA_NUM = 30;
    Utils::mkdir(time_.name());
    StringStream sPath;
    sPath << time_.name() << "/" << fieldName << ".plt";
    std::fstream fp;
    fp.open(sPath.str(), std::ios::out | std::ios::trunc);
    if (!fp.is_open()) {
        logger.warn << "Cannot open file: " << sPath.str() << std::endl;
        fp.close();
        FATAL_ERROR_THROW;
    }
    // write head
    fp << "TITLE = \"" << fieldName << "\"\n";
    fp << "FileType = \"SOLUTION\"\n"
          "VARIABLES = \"" << fieldName << "-X\", \"" << fieldName << "-Y\", \" " << fieldName << "-Z\"\n";
    fp << "ZONE N=" << mesh_.nodeNum() << ", E=" << mesh_.cellNum();
    fp << ", VARLOCATION=([1-3]=CELLCENTERED), DATAPACKING=BLOCK, ZONETYPE=";
    fp << (mesh_.dimension() == 2 ? "FEQUADRILATERAL" : "FEBRICK") << ", SOLUTIONTIME=" << time_.name() << std::endl;
    // write values
    int count;
    for (int i = 0; i < 3; ++i) {
        count = 0;
        fp << std::endl << "## " << fieldName << "-Component: " << i + 1 << std::endl;
        for (auto value: values) {
            fp << "\t" << std::setprecision(DATA_PRECISION) << value[i];
            if (count++ >= LINE_DATA_NUM) {
                fp << std::endl;
                count = 0;
            }
        }
    }
    // close
    fp.close();
    logger.note << "Write volVectorField data <" << fieldName << "> -> " << sPath.str() << std::endl;
}
