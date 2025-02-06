if (parser.parse_switch("help")) {
    logger.note << "============================================" << std::endl;
    logger.info << logo
                << "--------------------------------------------\n"
                   " * http://github.com/yuzuki01/meso-lts\n"
                << std::endl;
    return 0;
}