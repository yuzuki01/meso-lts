MESO::MPI::Initialize(&argc, &argv);

MESO::ArgParser parser(argc, argv);

if (parser.parse_switch("debug")) {
    logger.level = -1;
    debug = true;
}
MESO::logger.debug << "Running in debug mode." << std::endl;
