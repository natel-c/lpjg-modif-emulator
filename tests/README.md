LPJ-GUESS Unit Testing
======================

# Unit tests
Unit tests are small and fast tests to test the functionality of a single, small *unit* of the code.
This could be, for example, one function.

Testing serves multiple purposes:
1. Shows that a newly implemented function does what it's supposed to do
2. Makes others aware if their code changes break functionality
3. Serves as a great means of documentation for a piece of code
4. Aids at having understandable, short, concise and **testable** pieces of code!

# Writing tests
Ideally, a newly implemented function should come with at least one test.
A unit test typically has a few simple rules. It should

* be short
* be fast
* test one single thing, for example one function
* be understandable
* have an understandable name as to understand what it tests
* be able to run independenly of other tests, in particular tests must be able to run in any order

A test normally consists of
1. setup
2. call of the code under test
3. assertion(s) to verify the correctness of the code
4. (Clean-up)

Please refer to `ncompete_test.cpp` or `management_test.cpp` for examples of this.

Defining the "one single thing" you test is a matter of opinion, i.e., if you test a function with multiple inputs, you don't strictly have to put every combination of function call and assertion in its own test, but sometimes it is more natural to group a couple of function calls and assertions into one logical test. 

## Integration tests
These are tests that run more than just one component. For instance, `demo_integration_test.cpp` is an integration test that runs an entire LPJ-GUESS simulation.
Please note that when creating such tests, we should not add large files to the repository that the tests require to run.
It is currently under discussion how to deal with tests requiring large input files.


## A Few Additional Rules

#### Red to Green
Make sure you have seen a test fail before you commit it. That way you can make sure the test does not contain a bug and always passes (and thus does not test anything).

#### Test Independence, Set-Up, and Clean-Up
Tests must run completely independent of others. In particular, we must be able to run them in any order, because test frameworks will run them in random order.
One key to this that each tests properly sets up and cleans up the environment.
In LPJ-GUESS we need to specifically make sure that, e.g., global flags like `ifbvoc` are reset to their defaults if we choose to change them in a test.

In the Catch2 test framework, the differentiation between `TEST_CASE` and `SECTION` can be used for this.

## Catch2 Framework
Please refer to 
https://github.com/catchorg/Catch2/blob/devel/docs
for documentation on the test framework we use in LPJ-GUESS

### Tags
Catch2 provides the concept of tags. Using these we can tag tests with concepts of LPJ-GUESS they test, e.g. `harvest`. That way, if we change code in the management module, we can make use of that tag to run all tests that deal with harvesting.

## Benefit vs additional work
Writing tests along new code can actually help to develop this code and make it more clear and concise.
However since we work in a scientific environment we often write code that does not live very long, because we are trying many things out.
It is thus up to the developer to decide upon whether/how much testing is necessary.
Any code that is likely be used by others should be tested however.

# Running tests

## Running tests on build
Adding the flag `-DUNIT_TESTS=ON` to `cmake` (e.g. in your IDE or using `ccmake`) will result in all tests being run after every build.

## Running tests from within an IDE
This should work out of the box. If there are problems, you might need to restart the IDE and invalidate caches.
Running the tests from the IDE also allows for running it with a debugger.
In the IDE you can also select tags for which you want to write tests.