#ifndef MESO_TIMECONTROL_H
#define MESO_TIMECONTROL_H

namespace MESO {
    class Time {
    private:
        FileIO::ParamReader config_;
        Label step;
        Scalar deltaT_;
        Scalar time_;

        void operator=(const Time&);

    public:
        explicit Time(FileIO::ParamReader  config);

        ~Time() = default;

        void update();

        /// Interfaces

        String name() const;

        const Scalar& deltaT() const;
    };
}

#endif //MESO_TIMECONTROL_H
