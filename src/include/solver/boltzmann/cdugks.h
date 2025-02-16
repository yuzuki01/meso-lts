#ifndef MESO_SOLVER_CDUGKS_H
#define MESO_SOLVER_CDUGKS_H

class CDUGKS : public BasicSolver {
private:
    const Scalar R;
    const Scalar TRef;
    const Scalar muRef;

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

private:
    /// CDUGKS
    void gEqMaxwell(volScalarField &g_, const volScalarField &rho_, const volVectorField &U_, const Mesh::Cell &dv_);
};

#endif //MESO_SOLVER_CDUGKS_H
