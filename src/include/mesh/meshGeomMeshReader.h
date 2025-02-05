private:
    void parseGambit(const List<String>& lines);
    void parseGambitHead(Label& i, const Label& size,
                         const StringList& lines);
    void parseGambitPoint(Label& i, const Label& size,
                          const StringList& lines);
    void parseGambitCell(Label& i, const Label& size,
                         const StringList& lines);
