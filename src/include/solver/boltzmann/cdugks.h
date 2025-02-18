#ifndef MESO_SOLVER_CDUGKS_H
#define MESO_SOLVER_CDUGKS_H

class CDUGKS : public BasicSolver {
private:
    const Scalar R;
    const Scalar TRef;
    const Scalar muRef;

    Scalar dt_{};

    volScalarField rhoVol, rhoVolOld;
    volVectorField UVol, UVolOld;

    surfScalarField rhoSurf;
    surfVectorField USurf;

    GeomMesh DV_;
    List<volScalarField> gVol_;
    List<surfScalarField> gSurf_;
    List<volScalarField> fluxGVol_;

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
    Scalar gMaxwell(const Scalar &rho_, const Vector &U_, const Mesh::Cell &dv_);

    template<Label PatchType>
    BasicField<Scalar, PatchType> gMaxwell(const BasicField<Scalar, PatchType> &rho_,
                                           const BasicField<Vector, PatchType> &U_,
                                           const Mesh::Cell &dv_);

    template<Label PatchType>
    BasicField<Scalar, PatchType> tau_f(const BasicField<Scalar, PatchType> &rho_);

    void updateGbpSurf();

    void updateGSurf();

    void updateBC();

    void updateFVM();
};

#endif //MESO_SOLVER_CDUGKS_H
