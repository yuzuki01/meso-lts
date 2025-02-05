private:
    void parseGambitHead(Label& i, const Label& size,
                         const StringList& lines);
    void parseGambitNode(Label& i, const Label& size,
                          const StringList& lines);
    void parseGambitCell(Label& i, const Label& size,
                         const StringList& lines);
