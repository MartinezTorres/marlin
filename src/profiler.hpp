/***********************************************************************

profiler: a simple profiler for the Marlin codec (not thread safe)

Usage:

 	// For each event to be profiled
 	Profiler::start("event_name");
 	...
 	Profiler::end("event_name");

 	Profiler::report(std::cout);

Notes:

  * Events can be nested, but the start and end calls must be consistent.

  * If the same event name is started and ended several times, total duration is accumulated.

  * To completely disable, define the NO_PROFILER macro.

MIT License

Copyright (c) 2018 Manuel Martinez Torres, Miguel Hernández-Cabronero

Marlin: A Fast Entropy Codec

MIT License

Copyright (c) 2018 Manuel Martinez Torres

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

***********************************************************************/

#ifndef PROFILER_HPP
#define PROFILER_HPP

#include <time.h>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>

namespace marlin {

/**
 * Class to allow seamless profiling.
 *
 * This is not a thread-safe utility.
 */
class Profiler {

public:
	// Only the provided static methods should be used
	Profiler(const Profiler& other) = delete;
	void operator=(const Profiler& other) = delete;

	/**
	 * Start a new event now.
	 */
	static void start(std::string event_name);

	/**
	 * End the last started event. If a string is provided,
	 * it is verified that the finishing event has a matching name.
	 */
	static void end(std::string event_name="");

	/**
	 * Report the Profiler results to out.
	 *
	 * If csv_format is false, a hierarchical print is shown
	 * If csv_format is true, flattened results are printed in CSV format.
	 */
	static void report(std::ostream& out, bool csv_format=false);

	/**
	 * Create a text file at output_path and print the profiler report there.
	 *
	 * If csv_format is false, a hierarchical print is shown
	 * If csv_format is true, flattened results are printed in CSV format.
	 */
	static void report(std::string output_path, bool csv_format=false);


protected:
	// Only the provided static methods should be used
	Profiler() {}
};

}

#endif /* PROFILER_HPP */
