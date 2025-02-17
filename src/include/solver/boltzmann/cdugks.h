#ifndef MESO_SOLVER_CDUGKS_H
#define MESO_SOLVER_CDUGKS_H

class CDUGKS : public BasicSolver {
private:
    const Scalar R;
    const Scalar TRef;
    const Scalar muRef;

    volScalarField rhoVol, rhoVolOld;
    volVectorField UVol, UVolOld;

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
    template<Label PatchType>
    BasicField<Scalar, PatchType> gMaxwell(const BasicField<Scalar, PatchType> &rho_,
                                           const BasicField<Vector, PatchType> &U_,
                                           const Mesh::Cell &dv_);

    inline Scalar tau_f();
};

#endif //MESO_SOLVER_CDUGKS_H
