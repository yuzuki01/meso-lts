#include "core/core.h"


namespaceMESO


Time::Time(FileIO::ParamReader config)
        : config_(std::move(config)), step(0), time_(0.0) {
    deltaT_ = config_.get("solver", "deltaT", 1.0);
}

String Time::name() const {
    if (time_ == 0.0) return "0";
    StringStream ss;
    ss << std::fixed << std::setprecision(15) << time_;
    return ss.str();
}

const Scalar &Time::deltaT() const {
    return deltaT_;
}

void Time::update() {
    config_.update();
    deltaT_ = config_.get("solver", "deltaT", 1.0);
}
