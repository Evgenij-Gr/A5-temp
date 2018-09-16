#include "configuration.h"

Configuration::Configuration(int argc, char* argv[])
{
// the first argument is always the path to INI file
    char* pathToINIfile;
    pathToINIfile = argv[1];

    boost::property_tree::ptree config_file;
    boost::property_tree::ini_parser::read_ini(pathToINIfile, config_file);
    if (!checkConfigFileForConsistency(config_file))
    {
        std::cout<<"Config file is not consistent, terminating program"<<std::endl;
        isConfigurationConsistent = false;
    }
    else
    {
        isConfigurationConsistent = true;
    }
/////////////////////////////////////////////
    for (int i=2; i<argc; i++)
    {
        if (strcmp(argv[i], "--eigv")==0)
        {
            flags.isEigvalsNeeded = true;
        }
        if (strcmp(argv[i], "--lyap")==0)
        {
            flags.isLyapNeeded = true;
        }
        if (strcmp(argv[i], "--dist")==0)
        {
            flags.isDistNeeded = true;
        }
        if (strcmp(argv[i], "--plot")==0)
        {
            flags.isPlotNeeded = true;
        }
    }
/////////////////////////////////////////////
    std::string runMode = config_file.get<std::string>("mode");
/////////// ODE_System ///////////
    rInit = config_file.get<double>("OdeSystem.r");
    aInit = config_file.get<double>("OdeSystem.a");
    bInit = config_file.get<double>("OdeSystem.b");
/////////// Initial Fixed Point /////////
    initialFixedPoint = {config_file.get<double>("InitialFixedPoint.g1"),
                         config_file.get<double>("InitialFixedPoint.g2"),
                         config_file.get<double>("InitialFixedPoint.g3"),
                         config_file.get<double>("InitialFixedPoint.g4")};
/////////// Integrator Params ///////////

    rtol = config_file.get<double>("IntegratorParams.rtol");
    atol = config_file.get<double>("IntegratorParams.atol");
/////////// Event Params ////////////////
    eventTol = config_file.get<double>("EventParams.eventTol");
    approximateReturnTime = config_file.get<double>("EventParams.approximateReturnTime");
    skipTime = config_file.get<double>("EventParams.skipTime");
/////////// Fixed Point Params //////////
    fixPtTol = config_file.get<double>("FixedPointParams.fixPtTol");
    fixPtMaxIter = config_file.get<int>("FixedPointParams.fixPtMaxIter");
    jacobiNewtonEps = config_file.get<double>("FixedPointParams.jacobiNewtonEps");
    jacobiFixPtTypeEps = config_file.get<double>("FixedPointParams.jacobiFixPtTypeEps");
    fixPtEpsNewtonStep = config_file.get<double>("FixedPointParams.fixPtEpsNewtonStep");
/////////// Trajectory Params
    mapIterSkip = config_file.get<int>("TrajectoryParams.mapIterSkip");
    mapIterLast = config_file.get<int>("TrajectoryParams.mapIterLast");
/////////// Run Params ///////////
    if (runMode=="Single")
    {
        v2_N    = config_file.get<int>("Single.v_N")+1;
        v2_min  = config_file.get<double>("Single.v_min");
        v2_max  = config_file.get<double>("Single.v_max");
        v2_name = config_file.get<std::string>("Single.v_name");
        v2_att_stride = config_file.get<int>("Single.v_att_stride");

        v1_N    = 1;
        const std::array<std::string, 3> allParameterNames = {"r", "a", "b"};
        v1_name = *std::find_if(std::begin(allParameterNames), std::end(allParameterNames), [=](const std::string& s){return (v2_name != s);});
        v1_min  = config_file.get<double>("OdeSystem."+v1_name);
        v1_max  = config_file.get<double>("OdeSystem."+v1_name);
        v1_att_stride = 0;
    }
    else if (runMode=="Biparametric")
    {
        v1_N    = config_file.get<int>("Biparametric.v1_N")+1;
        v1_min  = config_file.get<double>("Biparametric.v1_min");
        v1_max  = config_file.get<double>("Biparametric.v1_max");
        v1_name = config_file.get<std::string>("Biparametric.v1_name");
        v1_att_stride = config_file.get<int>("Biparametric.v1_att_stride");

        v2_N    = config_file.get<int>("Biparametric.v2_N")+1;
        v2_min  = config_file.get<double>("Biparametric.v2_min");
        v2_max  = config_file.get<double>("Biparametric.v2_max");
        v2_name = config_file.get<std::string>("Biparametric.v2_name");
        v2_att_stride = config_file.get<int>("Biparametric.v2_att_stride");
    }
    else if (runMode=="Once")
    {
        v1_name = "r";
        v1_N = 1;
        v1_min = config_file.get<double>("OdeSystem.r");
        v1_max = config_file.get<double>("OdeSystem.r");
        v1_att_stride = 0;

        v2_name = "a";
        v2_N = 1;
        v2_min = config_file.get<double>("OdeSystem.a");
        v2_max = config_file.get<double>("OdeSystem.a");
        v2_att_stride = 0;
    }
    else if (runMode=="Flow")
    {
    /////////// Flow Params ///////////
        flowInitPoint = {config_file.get<double>("Flow.g1"),
                         config_file.get<double>("Flow.g2"),
                         config_file.get<double>("Flow.g3"),
                         config_file.get<double>("Flow.g4")};
        flowMaxTime = config_file.get<double>("Flow.maxTime");
        timeStep = config_file.get<double>("Flow.timeStep");
    }
// also the crossing direction can be defined beforehand
// a bad practice, really: assumes that we started from section, which isn't
// always true
    Ashwin5osc system(rInit, aInit, bInit);
    Ash_Section_Event evt;
    AshStateType gamma0 = initialFixedPoint;
    AshStateType dgammadt0;
    system(gamma0, dgammadt0, 0.0);
    AshStateType gradEvent0 = evt.gradient(gamma0);
    crossingDirection = dgammadt0[0]*gradEvent0[0]+dgammadt0[1]*gradEvent0[1]+dgammadt0[2]*gradEvent0[2]+dgammadt0[3]*gradEvent0[3];
    if (crossingDirection > 1e-15)
    {
        crossingDirection = 1.0;
    }
    else if (crossingDirection < -1e-15)
    {
        crossingDirection = -1.0;
    }
    else
    {
        std::cerr<<"Initial point has vector tangent to cross-section, terminating program"<<std::endl;
        isConfigurationConsistent = false;
    }
}

std::map<std::string, double> Configuration::getParameterValues(int i, int j) const
{
    std::map<std::string, double> params;
    // initial setup for parameters: taken from "default" section of INI-file
    params["r"] = rInit;
    params["a"] = aInit;
    params["b"] = bInit;
    if (v2_N != 1)
    {
        params[v2_name] = v2_min + i * (v2_max-v2_min) / (v2_N - 1);
    }
    if (v1_N !=1)
    {
        params[v1_name] = v1_min + j * (v1_max-v1_min) / (v1_N - 1);
    }
    return params;
}

bool Configuration::checkConfigFileForConsistency(boost::property_tree::ptree config) const
{
    bool hasNoErrors = true;
    if (config.get<std::string>("mode")=="Once")
    {
        ;
    }
    else if (config.get<std::string>("mode")=="Single")
    {
        if (config.get<double>("Single.v_min") != config.get<double>("OdeSystem."+config.get<std::string>("Single.v_name")))
        {
            std::cout<<"Config error: [Single.v_min] doesn't coincide with [OdeSystem."<<config.get<std::string>("Single.v_name")<<"]"<<std::endl;
            hasNoErrors = false;
        }
    }
    else if (config.get<std::string>("mode")=="Biparametric")
    {
        if (config.get<double>("Biparametric.v1_min") != config.get<double>("OdeSystem."+config.get<std::string>("Biparametric.v1_name")))
        {
            std::cout<<"Config error: [Biparametric.v1_min] doesn't coincide with [OdeSystem."<<config.get<std::string>("Biparametric.v1_name")<<"]"<<std::endl;
            hasNoErrors = false;
        }
        if (config.get<double>("Biparametric.v2_min") != config.get<double>("OdeSystem."+config.get<std::string>("Biparametric.v2_name")))
        {
            std::cout<<"Config error: [Biparametric.v2_min] doesn't coincide with [OdeSystem."<<config.get<std::string>("Biparametric.v2_name")<<"]"<<std::endl;
            hasNoErrors = false;
        }
    }
    return hasNoErrors;
}
