#include "solver/solver.h"

using namespace MESO::Solver;

void Patch::set_value(MESO::StringList &string_list) {
    /**
     * Patch format:
     * <name>   <patch-type>    <data-type>     <data-value>
     **/
    if (string_list.size() < 2) return;
    auto name = string_list[0];
    auto type_s = string_list[1];
    if (type_s == PatchType::calculated_str) {
        type[name] = PatchType::calculated;
    } else if (type_s == PatchType::fixedValue_str) {
        type[name] = PatchType::fixedValue;
        auto vts = string_list[2];
        if (vts == "objectType") {
            ints[name] = stoi(string_list[3]);
        } else if (vts == "scalar") {
            scalars[name] = stod(string_list[3]);
        } else if (vts == "vector") {
            vectors[name] = {stod(string_list[3]), stod(string_list[4]), stod(string_list[5])};
        } else {
            std::stringstream err;
            err << "Patch::set_value caught unsupported type: " << name << " - " << type_s;
            logger.warn << err.str() << std::endl;
            throw std::invalid_argument(err.str());
        }
    } else if (type_s == PatchType::zeroGradient_str) {
        type[name] = PatchType::zeroGradient;
    } else if (type_s == PatchType::fromFile_str) {
        type[name] = PatchType::fromFile;
        auto vts = string_list[2];
        auto file = string_list[3];
        if (vts == "objectType") {
            ints_list[name] = Utils::read_np_file<ObjectType>(file);
        } else if (vts == "scalar") {
            scalars_list[name] = Utils::read_np_file<Scalar>(file);
        } else if (vts == "vector") {
            vectors_list[name] = Utils::read_np_file<Vector>(file);
        } else {
            std::stringstream err;
            err << "Patch::set_value caught unsupported type: " << name << " - " << type_s;
            logger.warn << err.str() << std::endl;
            throw std::invalid_argument(err.str());
        }
    }
}


MESO::ObjectType Patch::get_type(const MESO::String &key) {
    return type[key];
}

MESO::ObjectType Patch::get_int(const MESO::String &key) {
    return ints[key];
}

MESO::ObjectType Patch::get_file_int(const MESO::String &key, MESO::ObjectId id) {
    return ints_list[key][id];
}

MESO::Scalar Patch::get_scalar(const MESO::String &key) {
    return scalars[key];
}

MESO::Scalar Patch::get_file_scalar(const MESO::String &key, MESO::ObjectId id) {
    return scalars_list[key][id];
}

MESO::Vector Patch::get_vector(const MESO::String &key) {
    return vectors[key];
}

MESO::Vector Patch::get_file_vector(const MESO::String &key, MESO::ObjectId id) {
    return vectors_list[key][id];
}
