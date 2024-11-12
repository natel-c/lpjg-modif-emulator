///////////////////////////////////////////////////////////////////////////////////////
/// \file string_test.cpp
/// \brief Unit tests functionality in guessstring.h
///
/// \author Joe Siltberg
/// $Date: 2024-06-14 16:07:56 +0200 (Fri, 14 Jun 2024) $
///
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at http://mozilla.org/MPL/2.0/.
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "catch.hpp"

#include "guessstring.h"

TEST_CASE("trim", "[string]") {
	REQUIRE(trim("") == "");

	REQUIRE(trim(" ") == "");
	
	REQUIRE(trim("   trim \t\n  ") == "trim");

	REQUIRE(trim("a b") == "a b");
	
	REQUIRE(trim(" a b ") == "a b");
}

TEST_CASE("to_upper", "[string]") {
	REQUIRE(to_upper("abc") == "ABC");

	REQUIRE(to_upper("a b c") == "A B C");

	REQUIRE(to_upper("1") == "1");
}

TEST_CASE("to_lower", "[string]") {
	REQUIRE(to_lower("AbC") == "abc");

	REQUIRE(to_lower("A B C") == "a b c");

	REQUIRE(to_lower("1") == "1");
}
