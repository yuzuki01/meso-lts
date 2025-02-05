#ifndef MESO_FILEREADER_H
#define MESO_FILEREADER_H

namespace MESO::FileIO {

    const String fileAnnotation = "#";

    class BasicReader {
    protected:
        const String file_;
    private:
        void operator=(const BasicReader&);
    public:
        explicit BasicReader(const String& filePath);

        ~BasicReader() = default;

        [[nodiscard]] StringList read_lines() const;
        [[nodiscard]] const String& name() const;
    };

    class ParamReader : private BasicReader {
    private:
        Dict<Dict<String>> data_;

        bool isVarExisted(const String& varRegion, const String& varName);

    public:

        explicit ParamReader(const String &filePath);

        void clear();

        template<typename MesoType> void set(const String& varRegion,
                                             const String& varName, const MesoType& varValue);

        void update();

        /// Interfaces
        template<typename MesoType> MesoType get(const String& varRegion, const String& varName,
                const MesoType& varDefault, bool throwNotFoundErr = true);
    };
}

#endif //MESO_FILEREADER_H
