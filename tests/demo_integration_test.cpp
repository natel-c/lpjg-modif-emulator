#if __cplusplus >= 201703L
#include "catch.hpp"

#include "commandlinearguments.h"
#include "framework.h"
#include "guess.h"
#include <fstream>
#include <iterator>
#include <iostream>

// for some reason, some people do not want to upgrade to newer C++ versions, so we need to add this check here...

#include <filesystem>

static const char *const CPOOL_OUTPUT_FILEPATH = "../tests/test_outputs/cpool.out";
static const double TOLERANCE_kgC = 5.0;

static std::string getLastLineOfFile(const char *filepath);

TEST_CASE("Simple demo run works for one grid cell and total carbon is in a reasonable range", "[integrationtest][demoinput]"){

    // remove cpool.out from test outputs, to make sure that we are really writing a cpool.out file and not checking a file that is lying around there from previous runs.
    std::error_code error_code;
    bool was_deleted = std::filesystem::remove(CPOOL_OUTPUT_FILEPATH, error_code);
    if(was_deleted){
        std::cout << "Output file " << CPOOL_OUTPUT_FILEPATH << " deleted. Ready to start the test." << std::endl;
    } else {
        if (error_code.value() != 0) {
            std::cout << "Output file " << CPOOL_OUTPUT_FILEPATH << " could not be deleted." << std::endl;
        } else{
            std::cout << "Output file " << CPOOL_OUTPUT_FILEPATH << " does not exist. This is not an error." << std::endl;
        }
    }

    std::cout << "Starting LPJ-GUESS" << std::endl;
    // this is a global variable. We should remove all pfts because it could interfere with other tests.
    pftlist.killall();
    char* castedArgs[4] = { const_cast<char*>("guess"), const_cast<char*>("-input"), const_cast<char*>("demo"), const_cast<char*>("../tests/test_insfiles/europe_demo_for_test.ins") };
    framework(CommandLineArguments(4, castedArgs));

    std::cout << "LPJ-GUESS has finished. Checking outputs..." << std::endl;
    // this is a global variable. We should remove all pfts because it could interfere with other tests.
    pftlist.killall();

    // check value of total carbon
    std::string lastLine = getLastLineOfFile(CPOOL_OUTPUT_FILEPATH);
    std::istringstream iss(lastLine);
    std::vector<std::string> tokens(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());

    double cpoolTotalAtSimulationEnd = std::stof(tokens.back());
    double expectedTotalCarbon = 25.0;
    REQUIRE(cpoolTotalAtSimulationEnd == Approx(expectedTotalCarbon).margin(TOLERANCE_kgC));
}

static std::string getLastLineOfFile(const char *filepath) {
    std::ifstream file(filepath);
    std::string lastLine;

    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
            lastLine = line;
        }
        file.close();
    }
    return lastLine;
}

#endif