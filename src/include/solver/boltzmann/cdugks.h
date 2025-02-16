#ifndef MESO_SOLVER_CDUGKS_H
#define MESO_SOLVER_CDUGKS_H

class CDUGKS : public BasicSolver {
private:
    volScalarField rho;
    volVectorField U;

    GeomMesh DV_;
    List<volScalarField> gVol_;
    List<surfScalarField> gSurf_;

public:
    /// Constructor
    explicit CDUGKS(const String &filePath);

    ~CDUGKS() = default;

    /// Initialize
    void initialize();

    /// Solution
    bool solution();

    /// File IO
    void output();
};

#endif //MESO_SOLVER_CDUGKS_H
