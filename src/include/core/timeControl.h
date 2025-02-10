#ifndef MESO_TIMECONTROL_H
#define MESO_TIMECONTROL_H

namespace MESO {
    class Time {
    private:
        FileIO::ParamReader config_;
        Label step_;
        Scalar deltaT_;
        Scalar time_;
        Scalar endTime_;

        void operator=(const Time&);

    public:
        explicit Time(FileIO::ParamReader config);

        ~Time() = default;

        void update();

        bool isStoppable() const;

        /// Interfaces
        void runStep(Label stepNum=1);

        String name() const;

        const Scalar& deltaT() const;

        Scalar time() const;

        Label step() const;
    };
}

#endif //MESO_TIMECONTROL_H
