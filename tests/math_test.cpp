///////////////////////////////////////////////////////////////////////////////////////
/// \file math_test.cpp
/// \brief Unit tests functionality in guessmath.h
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

#include "guessmath.h"

TEST_CASE("Historic/add", "[historic]") {
	Historic<double, 3> history;

	REQUIRE(history.size() == 0);

	history.add(1);

	REQUIRE(history.size() == 1);
	REQUIRE(history.sum() == Approx(1));
	REQUIRE(history[0] == 1);
	REQUIRE(history.mean() == Approx(1));

	history.add(2);

	REQUIRE(history.size() == 2);
	REQUIRE(history.sum() == Approx(3));
	REQUIRE(history[0] == 1);
	REQUIRE(history[1] == 2);
	REQUIRE(history.mean() == Approx(1.5));

	history.add(3);

	REQUIRE(history.size() == 3);
	REQUIRE(history.sum() == Approx(6));
	REQUIRE(history[0] == 1);
	REQUIRE(history[1] == 2);
	REQUIRE(history[2] == 3);
	REQUIRE(history.mean() == Approx(2));

	history.add(4);

	REQUIRE(history.size() == 3);
	REQUIRE(history.sum() == Approx(9));
	REQUIRE(history[0] == 2);
	REQUIRE(history[1] == 3);
	REQUIRE(history[2] == 4);
	REQUIRE(history.mean() == Approx(3));
}

TEST_CASE("variation_coefficient", "[variation]") {

	double single_value[] = { 7 };

	REQUIRE(variation_coefficient(single_value, 1) == Approx(-1));

	double two_values[] = { 7 , 7};

	REQUIRE(variation_coefficient(two_values, 2) == Approx(0));

	double values[] = { 2, 4, 4, 4, 5, 5, 7, 9 };

	REQUIRE(variation_coefficient(values, 8) == Approx(0.427618));
}
TEST_CASE("Historic/periodic", "[historic]") {
	Historic<double, 3> history;

	REQUIRE(history.size() == 0);

	history.add(1);
	history.add(2);
	history.add(3);

	REQUIRE(history.periodicsum(2) == Approx(5));
	REQUIRE(history.periodicmean(2) == Approx(2.5));
	REQUIRE(history.min() == Approx(1));
	REQUIRE(history.max() == Approx(3));
	history.add(4);
	REQUIRE(history.periodicsum(2) == Approx(7));
	REQUIRE(history.periodicmean(2) == Approx(3.5));
	REQUIRE(history.lastadd() == Approx(4));
	REQUIRE(history.min() == Approx(2));
	REQUIRE(history.max() == Approx(4));
}
