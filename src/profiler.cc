/***********************************************************************

profiler: a simple profiler for the Marlin codec (not thread safe)

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

#include "profiler.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#ifdef NO_PROFILER

namespace marlin {

// Empty implementations that can be easily factored out by the compiler.

void Profiler::start(std::string event_name) {}

void Profiler::end(std::string event_name) {}

void Profiler::report(std::ostream& out, bool csv_format) {}

void Profiler::report(std::string output_path, bool csv_format) {}

Profiler::Profiler() {}

}

#else

namespace {

	/// Clock types to be used for all events
	const std::vector< std::pair<const clockid_t, const std::string> > clock_types_names = {
			{CLOCK_PROCESS_CPUTIME_ID, "cpu"},
			{CLOCK_REALTIME, "wall"}
	};


	/// Represent a named event to be tracked by the profiler
	class Event {

	public:
		const std::string name;
		/// Number of previous instances of this event
		uint32_t times;
		/// Pointer to the most recently opened event when this is created,
		/// or nullptr if this event is not nested in any other.
		Event *const parent;

		static Event* getRoot();

		/**
		 * Start a child event with the given name and
		 * return its reference.
		 */
		Event* start_child(std::string name_);

		/**
		 * Start the event.
		 */
		void setStart();

		/**
		 * End this run of the event.
		 *
		 * Increase accumulated duration and time counter.
		 */
		void setEnd();

		/**
		 * Is this a finished event?
		 */
		bool finished();

		/**
		 * Recursively report the time measurements of this event and all descendents
		 * in CSV format.
		 */
		void report_csv(std::ostream& out);

		/**
		 * Recursively report the time measurements of this event and all descendents
		 * in plain-text format.
		 *
		 * @param indentation_level depth of the event in the event tree
		 */
		void report_plain(std::ostream& out, uint32_t indentation_level=0);

		/**
		 * Get a map of durations indexed by clock id.
		 * This can be called even for an "open" event, that is,
		 * out for which the end time has not been set. In that case,
		 * time is calculated up to now.
		 */
		std::map<clockid_t, double> getDurations();

		/**
		 * Recursively release any used memory
		 */
		~Event();

	protected:
		/// Clocks with start times
		std::map<clockid_t, timespec> start_clocks;
		/// Clocks with end times
		std::map<clockid_t, timespec> end_clocks;
		/// Accumulated durations in seconds per clock, not including the current duration.
		std::map<clockid_t, double> durations;
		/// Map of children indexed by name
		std::map<std::string, Event*> children_by_name;
		/// List of children names sorted by addition order
		std::vector<std::string> sorted_child_names;

		/**
		 * Set to now all the timespecs in clock_map
		 */
		void setClocksNow(std::map<clockid_t, timespec> &clock_map);

	private:
		/**
		 * Create an empty Event with durations set to 0,
		 * starting times set to now and no end times.
		 */
		Event(std::string name_, Event* parent_=nullptr);
	};

	Event* Event::start_child(std::string name_) {
		auto it_found_child = children_by_name.find(name_);
		Event * child;
		if (it_found_child != children_by_name.end()) {
			child = it_found_child->second;
			if (! child->finished()) {
				std::stringstream ss;
				ss << "Cannot re-start running event '" << name_ << "'" << std::endl;
				throw std::runtime_error(ss.str());
			}
			child->end_clocks.clear();
			child->setStart();
		} else{
			child = new Event(name_, this);
			children_by_name.emplace(name_, child);
			sorted_child_names.push_back(name_);
		}
		return child;
	}

	/**
	 * Start the event.
	 */
	void Event::setStart() {
		setClocksNow(start_clocks);
	}

	/**
	 * End this run of the event.
	 *
	 * Increase accumulated duration and time counter.
	 */
	void Event::setEnd() {
		times++;

		// Stop timer and calculate new durations
		setClocksNow(end_clocks);
		durations = getDurations();

		// Set start time as end time so that further calls
		// to getDuration do not add false time
		for (auto clockid_name : clock_types_names) {
			start_clocks.at(clockid_name.first) = end_clocks.at(clockid_name.first);
		}
	}

	/**
	 * Is this a finished event?
	 */
	bool Event::finished() {
		return ! end_clocks.empty();
	}

	/**
	 * Get a map of durations indexed by clock id.
	 * This can be called even for an "open" event, that is,
	 * out for which the end time has not been set. In that case,
	 * time is calculated up to now.
	 */
	std::map<clockid_t, double> Event::getDurations() {
		std::map<clockid_t, double> duration_map;

		for (auto it : clock_types_names) {
			clockid_t clock_id = it.first;

			// Previous duration
			double duration = 1.0 * durations.at(clock_id); // Avoid modifying reference

			// Current start time
			const auto it_start_time = start_clocks.find(clock_id);
			if (it_start_time == start_clocks.end()) {
				throw std::runtime_error("Error: event without start time");
			}
			const timespec start_time = it_start_time->second;

			// End times or now
			const auto it_end_time = end_clocks.find(clock_id);
			if (it_end_time == end_clocks.end()) {
				// Still running
				timespec now;
				clock_gettime(clock_id, &now);
				duration += (now.tv_sec - start_time.tv_sec) + 1E-9 * (now.tv_nsec - start_time.tv_nsec);
			} else {
				duration += (it_end_time->second.tv_sec - start_time.tv_sec) \
								 + 1e-9 * (it_end_time->second.tv_nsec - start_time.tv_nsec);
			}

			duration_map.emplace(clock_id, duration);
		}

		return duration_map;
	}

	Event::~Event() {
		for (auto it : children_by_name) {
			delete it.second;
		}
	}

	/**
	 * Set to now all the timespecs in clock_map
	 */
	void Event::setClocksNow(std::map<clockid_t, timespec> &clock_map) {
		for (auto it : clock_types_names) {
			const clockid_t& clock_id = it.first;
			timespec now;
			clock_gettime(clock_id, &now);
			auto emplace_it = clock_map.emplace(clock_id, now);
			if (! emplace_it.second) {
				// the entry already existed: replace it
				emplace_it.first->second = now;
			}
		}
	}

	Event::Event(std::string name_, Event* parent_) : name(name_), parent(parent_), times(0) {
		for (auto it : clock_types_names) {
			durations.emplace(it.first, 0.0);
		}
		setStart();
	}

	Event* Event::getRoot() {
		static Event root("total");
		return &root;
	}

	void Event::report_csv(std::ostream& out) {
		static const std::string separator(",");

		if (this == Event::getRoot()) {
			// Write the CSV header only for the root event
			out << "event_name" << separator << "times";
			for (auto clockid_name : clock_types_names) {
				out << separator << clockid_name.second;
			}
			out << std::endl;
		}

		// Report this event
		out << name << separator << times;
		auto durations_ = getDurations();
		for (auto clockid_name : clock_types_names) {
			out << separator << durations_.at(clockid_name.first);
		}
		out << std::endl;

		// Report all children
		for (auto child_name : sorted_child_names) {
			children_by_name.at(child_name)->report_csv(out);
		}
	}

	void Event::report_plain(std::ostream& out, uint32_t indentation_level) {
		static const std::string indentation("  ");

		for (uint32_t i=0; i<indentation_level; i++) {
			out << indentation;
		}

		std::map<clockid_t, double> durations_this = getDurations();
		std::map<clockid_t, double> durations_children = getDurations();

		// Print this node
		out << "[" << times << "] " << name << ":";
		for (auto clockid_name : clock_types_names) {
			out << " " << clockid_name.second << "=" << durations_this.at(clockid_name.first);
			durations_children.at(clockid_name.first) = 0.0;
		}
		out << std::endl;

		// Print all children and accumulate their total times
		for (auto child_name : sorted_child_names) {
			Event* child = children_by_name.at(child_name);
			child->report_plain(out, indentation_level+1);
			for (auto clockid_name : clock_types_names) {
				durations_children.at(clockid_name.first) += child->getDurations().at(clockid_name.first);
			}
		}

		// Print information of time unaccounted for by the children
		if (! sorted_child_names.empty()) {
			for (uint32_t i = 0; i < indentation_level + 1; i++) {
				out << indentation;
			}
			out << "(remaining@" << name << ") :";
			for (auto clockid_name : clock_types_names) {
				out << " " << clockid_name.second << "="
				    << durations_this.at(clockid_name.first) - durations_children.at(clockid_name.first);
			}
			out << std::endl;
		}
	}



	Event* current_event = Event::getRoot();
}

namespace marlin {

void Profiler::start(std::string event_name) {
	current_event = current_event->start_child(event_name);
}

void Profiler::end(std::string event_name) {
	if (current_event == Event::getRoot()) {
		throw std::runtime_error("Cannot end root event");
	}
	if (! event_name.empty()) {
		if (current_event->name != event_name) {
			std::stringstream ss;
			ss << "End name '" << event_name << "' does not match currently open event '"
			   << current_event->name << "'" << std::endl;
			throw std::runtime_error(ss.str());
		}
	}

	current_event->setEnd();
	current_event = current_event->parent;
}

void Profiler::report(std::ostream& out, bool csv_format) {
	if (csv_format) {
		current_event->report_csv(out);
	} else {
		current_event->report_plain(out);
	}

}



void Profiler::report(std::string output_path, bool csv_format) {
	if (! output_path.empty()) {
		std::ofstream out(output_path);
		Profiler::report(out, csv_format);
	}
}

}

#endif