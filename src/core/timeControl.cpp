#include "core/core.h"


namespaceMESO


Time::Time(FileIO::ParamReader config)
        : config_(std::move(config)),
          step_(0),
          time_(config_.get("solver", "startTime", 0.0)) {
    deltaT_ = config_.get("solver", "deltaT", 1.0);
    endTime_ = config_.get("solver", "endTime", 10.0);
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

Scalar Time::time() const {
    return time_;
}

Label Time::step() const {
    return step_;
}

void Time::update() {
    config_.update();
    deltaT_ = config_.get("solver", "deltaT", 1.0);
    endTime_ = config_.get("solver", "endTime", 10.0);
}

void Time::runStep(Label stepNum) {
    if (stepNum <= 0) {
        logger.error << "MESO::Time::runStep() caught invalid value" << std::endl;
        FATAL_ERROR_THROW;
    }
    int count = 0;
    while (isStoppable() and (count < stepNum)) {
        if ((time_ + deltaT_) > endTime_) {
            deltaT_ = endTime_ - time_;
        }
        time_ += deltaT_;
        step_++;
        count++;
    }
}

bool Time::isStoppable() const {
    return (time_ < endTime_);
}
