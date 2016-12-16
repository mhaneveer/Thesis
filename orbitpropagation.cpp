/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <iostream>
#include <cstdio>
#include <ctime>


boost::shared_ptr< tudat::propagators::PropagationTerminationSettings > getTerminationSettings( const int testType, const int days, const int start_epoch)
{
    // Define stopping conditions, depending on test case.
    boost::shared_ptr< tudat::propagators::PropagationTerminationSettings > terminationSettings;
    switch( testType )
    {
    // Stop at given time.
    case 0:
        terminationSettings = boost::make_shared< tudat::propagators::PropagationTimeTerminationSettings >( start_epoch + tudat::physical_constants::JULIAN_DAY * days  );
        break;

        // Stop at given altitude
    case 2:
        terminationSettings = boost::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                    boost::make_shared< tudat::propagators::SingleDependentVariableSaveSettings >(
                        tudat::propagators::altitude_dependent_variable, "Satellite" ), 80.0E3, 1 );


        break;

        // Stop when a single of the conditions 0-3 is fulfilled.
    case 4:
    {
        std::vector< boost::shared_ptr< tudat::propagators::PropagationTerminationSettings > > constituentSettings;
        constituentSettings.push_back( getTerminationSettings( 0 , days, start_epoch) );
        constituentSettings.push_back( getTerminationSettings( 2 , days, start_epoch) );

        terminationSettings = boost::make_shared< tudat::propagators::PropagationHybridTerminationSettings >(
                    constituentSettings, 1 );
        break;
    }
    }
    return terminationSettings;
}

//! Execute propagation of orbit of Satellite around the Earth.
int main()
{
    const int testType = 4;
    std::clock_t start;
    double duration;

    start = std::clock();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat;
    using namespace simulation_setup;
    using namespace propagators;
    using namespace numerical_integrators;
    using namespace orbital_element_conversions;
    using namespace basic_mathematics;
    using namespace gravitation;
    using namespace numerical_integrators;

    // Reading and Writing Data Files
    using tudat::input_output::DoubleKeyTypeVectorXdValueTypeMap;
    using tudat::input_output::writeDataMapToTextFile;
    using tudat::input_output::readMatrixFromFile;

    // Convert degrees to radians.
    using tudat::unit_conversions::convertDegreesToRadians;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     READ INPUT AND SETTING FILES         //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    const std::string& homePath = "D:/Documents/TU_Delft/MSc/Thesis/02_Program/TUDAT/tudatBundle/tudatApplications/Application2/TemplateApplication/";

    const std::string& filePath          = homePath + "data.txt";
    const std::string& file2Path         = homePath + "IC.txt";
    const std::string& filePathSetting   = homePath + "settings.txt";

    const std::string& WGS84path   = homePath + "WGS_84/settingsWGS84.txt";
    const std::string& cosinePath  = homePath + "WGS_84/cosine_coefficients.txt";
    const std::string& sinePath    = homePath + "WGS_84/sine_coefficients.txt";

    const std::string& separators = "\t";
    const std::string& skipLinesCharacter = "%";

    Eigen::MatrixXd FileData         =  readMatrixFromFile(filePath,        separators, skipLinesCharacter );
    Eigen::MatrixXd rvData           =  readMatrixFromFile(file2Path,       separators, skipLinesCharacter );
    Eigen::MatrixXd FileDataSetting  =  readMatrixFromFile(filePathSetting, separators, skipLinesCharacter );

    Eigen::MatrixXd WGS84_parameters    =  readMatrixFromFile(WGS84path,  ",", skipLinesCharacter );
    Eigen::MatrixXd cosineCoefficients  =  readMatrixFromFile(cosinePath, ",", skipLinesCharacter );
    Eigen::MatrixXd sineCoefficients    =  readMatrixFromFile(sinePath,   ",", skipLinesCharacter );

    const double gravitationalParameter  = WGS84_parameters(0);
    const double referenceRadius         = WGS84_parameters(1);

    //int rows = sizeof rvData - sizeof rvData[0];
    int rows = rvData(0,0);;

    //std::cout << cosineCoefficients   << std::endl << std::endl ;
    //std::cout << sineCoefficients     << std::endl << std::endl ;
    std::cout << rows <<std::endl ;
    //std::cout << rvData <<std::endl ;

    for (int ii = 1; ii < rows + 1; ii = ii + 1){

    const int days = FileDataSetting(0,0);

    std::cout << "Computing " << days << " days"  << std::endl;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels.
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );


    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back(  "Sun"   );
    bodiesToCreate.push_back(  "Earth" );
    bodiesToCreate.push_back(  "Moon"  );
    bodiesToCreate.push_back(  "Mars"  );
    bodiesToCreate.push_back(  "Venus" );



        // Set simulation time settings.
        const long double simulationStartEpoch = rvData(ii,10);
        const double start_epoch = simulationStartEpoch;
        const double simulationEndEpoch = tudat::physical_constants::JULIAN_DAY * days + simulationStartEpoch;

        // Create body objects.
        std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
                //getDefaultBodySettings( bodiesToCreate,  - 3000.0,  3000.0 );
               getDefaultBodySettings( bodiesToCreate, simulationStartEpoch, simulationEndEpoch);

        bodySettings[ "Earth" ]->atmosphereSettings = boost::make_shared< AtmosphereSettings >( nrlmsise00 );

        bodySettings[ "Earth" ]->gravityFieldSettings = boost::make_shared< GravityFieldSettings >( spherical_harmonic );
        bodySettings[ "Earth" ]->gravityFieldSettings = boost::make_shared< SphericalHarmonicsGravityFieldSettings >(
          gravitationalParameter,
          referenceRadius,
          cosineCoefficients,
          sineCoefficients,
          "IAU_Earth" );

        for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
        {
            bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
            bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
        }

        NamedBodyMap bodyMap = createBodies( bodySettings );

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Create spacecraft object.
        double referenceArea = rvData(ii,7);
        double scMass = rvData(ii,8);
        double aerodynamicCoefficient = rvData(ii,9);

        bodyMap[ "Satellite" ] = boost::make_shared< simulation_setup::Body >( );
        bodyMap[ "Satellite" ]->setConstantBodyMass( scMass );

        boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
                boost::make_shared< ConstantAerodynamicCoefficientSettings >(
                    referenceArea, aerodynamicCoefficient * Eigen::Vector3d::UnitX( ), 1, 1 );
        bodyMap[ "Satellite" ]->setAerodynamicCoefficientInterface(
                    createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Satellite" ) );

        boost::shared_ptr< RadiationPressureInterfaceSettings > SatelliteRadiationPressureSettings =
                boost::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                    "Sun", referenceArea, aerodynamicCoefficient, boost::assign::list_of( "Earth" ) );

        bodyMap[ "Satellite" ]->setRadiationPressureInterface(
                    "Sun", createRadiationPressureInterface(
                     SatelliteRadiationPressureSettings, "Satellite", bodyMap ) );

        // Finalize body creation.
        setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Define propagator settings variables.
        SelectedAccelerationMap    accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;

        // Define propagation settings.
        std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfSatellite;
        accelerationsOfSatellite[ "Earth" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 5, 5 ) );

        accelerationsOfSatellite[ "Sun"   ].push_back( boost::make_shared< AccelerationSettings >(
                                                         basic_astrodynamics::central_gravity ) );
        accelerationsOfSatellite[ "Moon"  ].push_back( boost::make_shared< AccelerationSettings >(
                                                         basic_astrodynamics::central_gravity ) );
        accelerationsOfSatellite[ "Mars"  ].push_back( boost::make_shared< AccelerationSettings >(
                                                         basic_astrodynamics::central_gravity ) );
        accelerationsOfSatellite[ "Venus" ].push_back( boost::make_shared< AccelerationSettings >(
                                                         basic_astrodynamics::central_gravity ) );
        accelerationsOfSatellite[ "Sun"   ].push_back( boost::make_shared< AccelerationSettings >(
                                                         basic_astrodynamics::cannon_ball_radiation_pressure ) );
        accelerationsOfSatellite[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >(
                                                         basic_astrodynamics::aerodynamic ) );

        accelerationMap[  "Satellite" ] = accelerationsOfSatellite;
        bodiesToPropagate.push_back( "Satellite" );
        centralBodies.push_back( "Earth" );
 std::cout << "Computing " << days << " 11"  << std::endl;
        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodyMap, accelerationMap, bodiesToPropagate, centralBodies );
 std::cout << "Computing " << days << " 12"  << std::endl;
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Define list of dependent variables to save.
        std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;

        dependentVariables.push_back(
                    boost::make_shared< SingleDependentVariableSaveSettings >( altitude_dependent_variable,
                                                                               "Satellite", "Earth" ) );

        // --------------------------------------------------------------------------------------------
        // READ INITIAL STATE FROM CARTESIAN ELEMENTS
        // --------------------------------------------------------------------------------------------

        Vector6d SatelliteInitialState;
        SatelliteInitialState( 0 ) = rvData(ii,1);
        SatelliteInitialState( 1 ) = rvData(ii,2);
        SatelliteInitialState( 2 ) = rvData(ii,3);
        SatelliteInitialState( 3 ) = rvData(ii,4);
        SatelliteInitialState( 4 ) = rvData(ii,5);
        SatelliteInitialState( 5 ) = rvData(ii,6);

        boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                boost::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, SatelliteInitialState,getTerminationSettings( testType, days, start_epoch ),  cowell,
                  boost::make_shared< DependentVariableSaveSettings >( dependentVariables ) );

        boost::shared_ptr< IntegratorSettings< > > integratorSettings =
                boost::make_shared< RungeKuttaVariableStepSizeSettings < > >
                ( numerical_integrators::rungeKuttaVariableStepSize, simulationStartEpoch, 60.0 ,
                  numerical_integrators::RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg78, 1, 300,
                  1.0E-12,
                  1.0E-12,
                  1,
                  0.8,
                  4.0,
                  0.1 );

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator(
                    bodyMap, integratorSettings, propagatorSettings, true, false, false );
        std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

        std::map< double, Eigen::VectorXd > dependentVariableSoution = dynamicsSimulator.getDependentVariableHistory( );

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        std::cout << "Writing Data" << std::endl;
        // Write Satellite propagation history to file.
        writeDataMapToTextFile( integrationResult,
                                "RK78_" + boost::lexical_cast< std::string >( rvData(ii,0) ) +  ".dat",
                                homePath + "/RK78",
                                "",
                                std::numeric_limits< double >::digits10,
                                std::numeric_limits< double >::digits10,
                                "," );

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.

    std::cout << std::endl << "Initial State (r,v) in m and m/s:  " << std::endl << SatelliteInitialState << std::endl << std::endl;
    std::cout << "Reference area set to " <<  referenceArea << std::endl;
    std::cout << "Spacecraft mass set to " <<  scMass << std::endl;
    std::cout << "Initial drag coefficient set to " <<  aerodynamicCoefficient << std::endl << std::endl;
    std::cout << "Finished satellite " << ii <<" out of " << rows  <<std::endl << std::endl;

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;


    std::cout<<"Total computation time: "<< duration <<'\n' << std::endl;
    }
    return EXIT_SUCCESS;
}

