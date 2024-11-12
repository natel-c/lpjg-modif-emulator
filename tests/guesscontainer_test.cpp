///////////////////////////////////////////////////////////////////////////////////////
/// \file guesscontainer_test.cpp
/// \brief Unit tests GuessContainer
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

#include "guesscontainer.h"

TEST_CASE("indexing", "[guesscontainer]") {
	GuessContainer<int> container;

	container.push_back(new int(1));
	container.push_back(new int(2));

	REQUIRE(container[0] == 1);
	REQUIRE(container[1] == 2);
	REQUIRE(container.size() == 2);
}

class InstanceCountingClass {
public:
	static int living_instances;

	InstanceCountingClass() {
		++living_instances;
	}

	~InstanceCountingClass() {
		--living_instances;
	}
};

int InstanceCountingClass::living_instances = 0;

TEST_CASE("memory", "[guesscontainer]") {

	GuessContainer<InstanceCountingClass> container;

	container.push_back(new InstanceCountingClass());
	container.push_back(new InstanceCountingClass());
	container.push_back(new InstanceCountingClass());
	
	REQUIRE(InstanceCountingClass::living_instances == 3);

	container.erase(container.begin());

	REQUIRE(InstanceCountingClass::living_instances == 2);

	container.clear();

	REQUIRE(InstanceCountingClass::living_instances == 0);
}

TEST_CASE("iteration", "[guesscontainer]") {
	GuessContainer<int> container;

	container.push_back(new int(1));
	container.push_back(new int(2));

	GuessContainer<int>::iterator itr = container.begin();

	REQUIRE(itr == container.begin());
	REQUIRE(itr++ == container.begin());
	REQUIRE(*itr == 2);

	itr = container.begin();

	REQUIRE(++itr != container.begin());
	REQUIRE(*itr == 2);
	
	REQUIRE(++itr == container.end());
}

TEST_CASE("erase", "[guesscontainer]") {
	GuessContainer<int> container;

	container.push_back(new int(1));
	container.push_back(new int(2));
	container.push_back(new int(3));
	
	GuessContainer<int>::iterator itr = container.erase(container.begin());

	REQUIRE(itr == container.begin());

	++itr;

	itr = container.erase(itr);

	REQUIRE(itr == container.end());
	
	itr = container.begin();
	itr = container.erase(itr);
	
	REQUIRE(itr == container.end());
	REQUIRE(container.size() == 0);
}
