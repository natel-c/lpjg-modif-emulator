///////////////////////////////////////////////////////////////////////////////////////
/// \file main.cpp
/// \brief Main function for the unit tests
///
/// \author Joe Siltberg
/// $Date: 2024-06-14 16:07:56 +0200 (Fri, 14 Jun 2024) $
///
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at http://mozilla.org/MPL/2.0/.
///
///////////////////////////////////////////////////////////////////////////////////////

// CATCH_CONFIG_NO_POSIX_SIGNALS is needed because catch.hpp uses 'SIGSTKSZ' constant which in newer compilers is no longer available.
#define CATCH_CONFIG_NO_POSIX_SIGNALS
#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include "shell.h"

int main(int argc, char** argv) {

	// Set a shell so we have working dprintf, fail etc.
	set_shell(new CommandLineShell("tests.log"));

	// Let CATCH do the rest
	int result = Catch::Session().run(argc, argv);

	return result;
}
