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

    namespace Flag {
        enum {string, scalar, vector};
        const Dict<ObjectType> TypeName = {
                {"string", string},
                {"scalar", scalar},
                {"vector", vector}
        };
    }

    class ParamReader : public BasicReader {
    public:
        class PatchParam {
        protected:
            Dict<String> dataStr_;
            Dict<Scalar> dataScl_;
            Dict<Vector> dataVtr_;
        public:

            PatchParam() = default;

            String name();

            String type();

            void set(const String &line, Label lineNo=0);

            template<typename T>
            T get(const String &name);
        };

    protected:
        Dict<Dict<String>> data_;
        Dict<PatchParam> zones_;
        Dict<PatchParam> marks_;

        bool isVarExisted(const String& varRegion, const String& varName);

    public:

        explicit ParamReader(const String &filePath);

        void clear();

        template<typename MesoType> void set(const String& varRegion,
                                             const String& varName, const MesoType& varValue);

        void update();

        /// Interfaces
        template<typename MesoType> MesoType get(const String& varRegion, const String& varName,
                                                 const MesoType &varDefault, bool throwNotFoundErr = true);

        const PatchParam &zone(const String &name);

        const PatchParam &mark(const String &name);
    };
}

#endif //MESO_FILEREADER_H
